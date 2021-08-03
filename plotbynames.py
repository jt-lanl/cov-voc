import sys
import re
from collections import defaultdict,Counter
import matplotlib.pyplot as plt
import datetime
import argparse

import covid
import colornames
#import readseq
#import sequtil
import embersplot
#from embers import get_daterange,date_from_seqname

OTHER = 'other'

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--weekly",action="store_true",
        help="Make weekly average plots instead of cumulative")
    paa("--daily",action="store_true",
        help="Make daily plots instead of cumulative")
    paa("--lineplot",action="store_true",
        help="Make log-linear line plot instead of linear stacked bar plot")
    paa("--onsets",action="store_true",
        help="plot onset dates for each mutant")
    #paa("--nolegend",action="store_true",help="avoid putting legend on plot")
    paa("--legend",type=int,default=0,choices=(0,1,2),
        help="0: no legend, 1: legend, 2: big legend (with seq patterns)")
    paa("--output","-o",help="write plot to file")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    if (args.daily and args.weekly):
        raise RuntimeError("Pick only one of daily or weekly")
    return args


def date_fromiso(s):
    if isinstance(s,datetime.date):
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError:
        #print("Invalid date:",s)
        return None

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
    seqs = covid.filter_seqs(seqs,args)
    seqs = covid.checkseqlengths(seqs)

    colors = {OTHER       : 'LightGrey',
              'B.1.1.7'   : 'Orange', #Alpha
              'B.1.351'   : 'Plum', #Beta
              'P.1'       : 'FireBrick', #Gamma
              'C.37'      : 'Green', #Lambda
              'B.1.617.2' : 'BlueViolet', #Delta
              'AY.1'      : 'Orchid',
              'AY.2'      : 'LightPink',
              'AY.3'      : 'Red',
    }
    patterns = list(colors)
    for p in patterns:
        colors[p] = colornames.tohex(colors[p])
              
    DG_datecounter = {m: Counter() for m in patterns} 
    for s in seqs:

        seqdate = date_from_seqname(s.name)
        if not seqdate:
            continue

        if (seqdate.year,seqdate.month) < (2019,11):
            print("bad seqdate:",seqdate)
            continue

        vocmatch = [voc for voc in patterns[1:] if re.search(r'\.'+voc+'$',s.name)]
        for voc in vocmatch:
            DG_datecounter[voc][seqdate] += 1
        if not vocmatch:
            DG_datecounter[OTHER][seqdate] += 1
            
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
   
    DG_cum=defaultdict(list)
    for ord in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord)
        for m in DG_datecounter:
            DG_cum[m].append( sum(DG_datecounter[m][dt] for dt in DG_datecounter[m] if dt <= day ) )

    if args.weekly or args.daily: ## Weekly averages:
        DAYSPERWEEK=7 if args.weekly else 1
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
        embersplot.embersplot(DG_cum,None,colors,ordmin,
                              ordplotrange = (ordplotmin,ordplotmax),
                              title = covid.get_title(args),
                              legend=1,fraction=f,linePlot=args.lineplot)

        plt.ylabel("Fraction" if f else "Counts")
        plt.tight_layout() ## do it again to accomodate the ylabel

        linbar = "line" if args.lineplot else "bar"
        if args.output:
            fc = "f" if f else "c"
            wk = "wk" if args.weekly else "dy"
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
    

