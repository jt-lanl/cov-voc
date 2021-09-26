import sys
import re
from collections import defaultdict,Counter
import matplotlib.pyplot as plt
import random
import datetime
import argparse

import covid
import colornames
#import readseq
import sequtil
import embersplot
#from embers import get_daterange,date_from_seqname

OTHER = 'other'

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--weekly",action="store_true",
        help="Make weekly average plots instead of cumulative")
    paa("--daily",type=int,default=7,
        help="daily=1 for daily, daily=7 (default) for weekly, daily=0 for cumulative")
    paa("--lineplot",action="store_true",
        help="Make log-linear line plot instead of linear stacked bar plot")
    paa("--onsets",action="store_true",
        help="plot onset dates for each mutant")
    #paa("--nolegend",action="store_true",help="avoid putting legend on plot")
    paa("--legend",type=int,default=0,choices=(0,1,2),
        help="0: no legend, 1: legend, 2: big legend (with seq patterns)")
    paa("--other",action="store_true",
        help="write out the lineages in the 'other' class")
    paa("--output","-o",help="write plot to file")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    if args.weekly:
        warnings.warn("'--weekly' deprecated: weekly is now default")
    if (args.daily and args.weekly):
        raise RuntimeError("'--weekly' deprecated; use '--daily 1' for daily, "
                           "default is '--daily 7' which is weekly")
    
    return args


def xdate_fromiso(s):
    return sequtil.date_fromiso(s)

def days_in_month(yyyy,mm):
    # https://stackoverflow.com/questions/42950/how-to-get-the-last-day-of-the-month
    next_month = datetime.date(int(yyyy),int(mm),28) + datetime.timedelta(days=4)
    return (next_month - datetime.timedelta(days=next_month.day)).day

def date_fromiso(s):
    if type(s) == datetime.date:
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError: 
        if s == ".":
            return None
        try:
            yyyy,mm = s.split("-")
            day = random.randint(1,days_in_month(yyyy,mm))
            dt = datetime.date(int(yyyy),int(mm),day)
            return dt
        except ValueError as e:
            return None
    return None #raise RuntimeError(f"Invalid Date {s}")


def date_from_seqname(sname):
    ## a not-very-robust way to get the date
    ## alternative would be to search for /d/d/d/d-/d/d-/d/d
    tokens = sname.split('.')
    try:
        date = date_fromiso(tokens[4])
    except IndexError:
        date = None
    return date

def get_daterange(datecounter,argsdates):
    ''' find range of dates, in ordinal numbers, based on:
    datecounter[mutant][date] = count of mutants in date, and
    argsdates is basically args.dates '''
    
    all_dateset=set()
    for p in datecounter:
        all_dateset.update(datecounter[p])

    vprint("Range of dates:",min(all_dateset),max(all_dateset))
    ordmin = min(all_dateset).toordinal()
    ordmax = max(all_dateset).toordinal()

    ordplotmin = ordmin
    ordplotmax = ordmax
    if argsdates and argsdates[0] and argsdates[0] != ".":
        ordplotmin = date_fromiso(argsdates[0]).toordinal()
    if argsdates and argsdates[1] and argsdates[1] != ".":
        ordplotmax = date_fromiso(argsdates[1]).toordinal()

    ordmin = min([ordmin,ordplotmin])
    ordmax = max([ordmax,ordplotmax])

    return ordmin, ordmax, ordplotmin, ordplotmax



def main(args):
    seqs = covid.read_seqfile(args)
    ## filter by country, but not by date
    seqs = covid.filter_seqs_by_pattern(seqs,args)
    seqs = covid.filter_seqs_by_date(seqs,args)
    seqs = covid.fix_seqs(seqs,args)
    seqs = sequtil.checkseqlengths(seqs)

    table = [
        ('other',   OTHER,                    'LightGrey'),
        ('Alpha',  'B.1.1.7',                 'Orange'),
        ('Beta',   'B.1.351',                 'Plum'),
        ('Gamma',  'P.1.*',                   'FireBrick'),
        ('Delta',  r'(B\.1\.617\.2)|(AY\.[0-9]+)', 'BlueViolet'),
        ('Lambda', 'C.37',                    'Green'),
        ('Mu',     'B.1.621(.1)?',            'Black'),
    ]

    patterns = []
    colors = dict()
    fullnames = dict()

    for name,patt,color in table:
        colors[patt] = colornames.tohex(color)
        patterns.append(patt)
        fullnames[patt] = name
    
    DG_datecounter = {m: Counter() for m in patterns}
    other_lineages = Counter()
    for s in seqs:

        seqdate = date_from_seqname(s.name)
        if not seqdate:
            continue

        if (seqdate.year,seqdate.month) < (2019,11):
            print("bad seqdate:",seqdate)
            continue

        vocmatch = [voc for voc in patterns[1:]
                    if re.search(r'\.'+voc+'$',s.name)]
        for voc in vocmatch:
            DG_datecounter[voc][seqdate] += 1
        if not vocmatch:
            if args.other:
                lineage = covid.get_lineage_from_name(s.name)
                other_lineages[lineage] += 1
                vprint(s.name)
            DG_datecounter[OTHER][seqdate] += 1

    if args.other:
        otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
        for otherlin in otherlist:
            print("%6d %s" % (other_lineages[otherlin],otherlin))
            
    nmatches = sum(sum(DG_datecounter[p].values()) for p in patterns)
    vprint("matched sequences:",nmatches)
    if nmatches==0:
        raise RuntimeError("No sequences for: " + " ".join(args.filterbyname))

    for p in DG_datecounter:
        vprint(p,sum(DG_datecounter[p].values()))

    ## Onset times for each pattern
    ## Don't include OTHER or patterns that don't appear in sequence set
    onset = {m: min(DG_datecounter[m]) for m in patterns if DG_datecounter[m]}
    for m in onset:
        vprint("onset",m,onset[m])

    ordmin, ordmax, ordplotmin, ordplotmax = get_daterange(DG_datecounter,args.dates)
    Ndays = ordmax+1-ordmin
    vprint("Days:",Ndays,ordmin,ordmax)
   
    DG_cum={m: list() for m in patterns} #defaultdict(list)
    for ord in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord)
        for m in DG_datecounter:
            DG_cum[m].append( sum(DG_datecounter[m][dt] for dt in DG_datecounter[m] if dt <= day ) )

    if args.daily: ## daily=7 (default) for weekly averages
        DAYSPERWEEK=args.daily
        DG_weekly=dict()
        for m in DG_cum:
            DG_weekly[m] = DG_cum[m][:]
            for n in range(DAYSPERWEEK,Ndays):
                DG_weekly[m][n] = DG_cum[m][n] - DG_cum[m][n-DAYSPERWEEK]
            DG_cum[m] = DG_weekly[m]
            
    ## Only keep data that is within the specified date range
    ## That way, automatic scaling on the y-axis will be based on available data
    if args.dates:
        for m in DG_cum:
            ztmp = []
            for i,cnt in enumerate(DG_cum[m]):
                if ordplotmin <= ordmin+i <= ordplotmax:
                    ztmp.append(cnt)
            DG_cum[m] = ztmp
        ordmin = ordplotmin
        Ndays = ordplotmax+1-ordmin

    for f in (False,True):
        embersplot.embersplot(DG_cum,fullnames,colors,ordmin,
                              ordplotrange = (ordplotmin,ordplotmax),
                              title = covid.get_title(args) + ": %d sequences" % (nmatches,),
                              legend=1,fraction=f,lineplot=args.lineplot)

        plt.ylabel("Fraction" if f else "Counts")
        plt.tight_layout() ## do it again to accomodate the ylabel

        linbar = "line" if args.lineplot else "bar"
        if args.output:
            fc = "f" if f else "c"
            wk = "wk" if args.daily==7 \
                else "dy" if args.daily == 1 \
                     else "cm"
            
            plt.savefig("-".join([fc,wk,linbar,args.output]))

    if not args.output:
        plt.show()
    
 

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

