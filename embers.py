'''
Stacked barplots (also, optionally, line-plots) of variant counts over time
'''

import sys
import warnings

import re
from collections import Counter,defaultdict
import datetime

import argparse
import matplotlib.pyplot as plt

import sequtil
import intlist
from spikevariants import SpikeVariants
import embersplot
import covid

def _getargs():
    ''' read command line arguments using argparse package'''
    ap = argparse.ArgumentParser(description=__doc__,
                                 conflict_handler='resolve')
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
    paa("--legend",type=int,default=0,choices=(0,1,2),
        help="0: no legend, 1: legend, 2: big legend (with seq patterns)")
    paa("--colormut","-c",
        help="read SpikeVariants structure from color-mut file")
    paa("--ctable","-t",
        help="write a count table to this file")
    paa("--output","-o",help="write plot to file")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()

    if args.weekly:
        warnings.warn("'--weekly' deprecated: weekly is now default")
    if (args.daily and args.weekly):
        raise RuntimeError("'--weekly' deprecated; use '--daily 1' for daily")
    return args

def date_fromiso(s):
    ''' return datetime.date value based on yyyy-mm-dd string '''
    if isinstance(s,datetime.date):
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError:
        #print("Invalid date:",s)
        return None

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

def filename_prepend(pre,file):
    '''prepend a string to a file name; eg
       "pre","file" -> "prefile", but also
       "pre","dir/file" -> "dir/prefile"
    '''
    if not file:
        return file
    return re.sub(r"(.*/)?([^/]+)",r"\1"+pre+r"\2",file)

def relativepattern(master,mutant,dittochar='_',noditto='-'):
    '''ABC,ABD -> __D'''
    return "".join((dittochar if (a==b and a not in noditto) else b)
                   for a,b in zip(master,mutant))

def reltoabspattern(master,mutant,dittochar='_'):
    '''ABC,__D -> ABD'''
    return "".join((a if b==dittochar else b)
                   for a,b in zip(master,mutant))


def check_dups(xlist):
    '''see if there are any duplicates in the list xlist'''
    duplist=[]
    xset=set()
    for x in xlist:
        if x in xset:
            duplist.append(x)
        xset.add(x)
    return duplist

def lineage_counts(sitelist,master,voclist,cpatt,n_sequences):
    ''' yield print strings for lineage counts '''
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent Lineage"
    for voc in voclist[::-1]:
        patt = voc.flat_pattern
        fpatt = cpatt[patt]/n_sequences if n_sequences else 0
        rpatt = relativepattern(master,patt)
        yield "%s %6d %6.2f%% %s" % (rpatt,cpatt[patt],100*fpatt,voc.name)

def missing_patterns_with_nearby(sitelist,master,voclist,xpatt,n_sequences):
    ''' for each of the seq patterns in xpatt, find nearby patterns in voclist '''
    ## routine yields lines that are meant to be printed
    yield f"Missing: {len(xpatt)} patterns, for a total of {sum(xpatt.values())} sequences"
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent NearbyLineage(s)"
    for seq in sorted(xpatt,key=xpatt.get,reverse=True)[:50]: #hardcoded 15
        rseq = relativepattern(master,seq)
        dv = {str(v): sum(bool(x != y and y != ".")
                          for x,y in zip(seq,v.flat_pattern))
              for v in voclist}
        vnearby_names = [v.name.strip() for v in voclist if dv[str(v)]<2]
        fpatt = xpatt[seq]/n_sequences if n_sequences else 0
        yield "%s %6d %6.2f%% %s" % (rseq,xpatt[seq],100*fpatt,
                             ", ".join(vnearby_names))

def main(args):
    ''' embers main '''
    args_keepx = args.keepx
    args.keepx = False
    seqs = covid.read_filter_seqfile(args)
    args.keepx = args_keepx
    first,seqs = sequtil.get_first_item(seqs)

    if args.colormut:
        svar = SpikeVariants.from_colormut(args.colormut,refseq=first.seq)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar = SpikeVariants.default(refseq=first.seq)

    OTHERNAME  = svar.OTHERNAME
    OTHERCOLOR = svar.OTHERCOLOR

    ## try this?? makes very little difference!!
    # svar.less_exact()

    master = svar.master
    sitelist = svar.ssites()
    voclist = svar.vocs

    for v in voclist:
        v.flat_pattern = svar.flatpattern(v)
        vvprint(f"{v} {v.name}: {v.exact} {v.flat_pattern}")

    mutants = [v.flat_pattern for v in voclist]
    patterns = mutants + [OTHERNAME]
    dups = check_dups(patterns)
    if dups:
        raise RuntimeError(f"Duplicated patterns {dups}")

    colors = [v.color for v in voclist] + [OTHERCOLOR]
    mcolors = dict(zip(patterns,colors))
    dups = check_dups(colors)
    if dups:
        vprint("Duplicated colors:",dups)

    namelist = [v.name for v in voclist] + [OTHERNAME]
    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    for v in voclist:
        v.name = namefmt % v.name
    mnames =  {m: namefmt%n for m,n in zip(patterns,namelist)}
    dups = check_dups(namelist)
    if dups:
        vprint("Duplicated names:",dups)

    def relpattern(mut):
        if mut == OTHERNAME:
            return "." * len(master)
        return relativepattern(master,mut,dittochar='_')

    mrelpatt = {p: relpattern(p) for p in patterns}

    for p in patterns:
        vprint(mnames[p],mrelpatt[p],mcolors[p])

    if args.legend == 0:
        fullnames = None
        legendtitle = None
    if args.legend == 1:
        fullnames = mnames
        legendtitle = None
    if args.legend == 2:
        fullnames = {p: " ".join([mnames[p],mrelpatt[p]]) for p in patterns}
        legendtitle = " "*(maxnamelen+2) + master

    svar.checkmaster(first.seq) ## ensure master agrees with first seqlist

    seqlist = list(seqs)
    n_sequences = len(seqlist)-1  ## -1 not to count the reference sequence

    if not args.keepx:
        seqlist = [s for s in seqlist if "X" not in s.seq]
        vprint("Removed",n_sequences+1-len(seqlist),"sequences with X")
        n_sequences = len(seqlist)-1

    ## How many of each sequence
    c = Counter(s.seq for s in seqlist[1:])

    ## How many of each mutant
    cpatt = Counter()
    xpatt = Counter()
    for seq in c:

        vocmatch = svar.vocmatch(seq)

        for voc in vocmatch[:1]:
            cpatt[voc.flat_pattern] += c[seq]

        ## Ideally just one match, if zero or more than one, then...
        if len(vocmatch)==0:
            sseq = svar.shorten(seq)
            xpatt[sseq] = c[seq]
        elif len(vocmatch)>1:
            if any(voc.name != vocmatch[0].name for voc in vocmatch):
                warn_msg = f"\n{svar.shorten(seq)} matches\n"
                warn_msg += " and\n".join(f"{relpattern(v.flat_pattern)} {v.name} {v}"
                                          for v in vocmatch)
                warn_msg += f"\n{svar.master} Master"
                warnings.warn(warn_msg)

    vprint("Unmatched sequences:",sum(xpatt.values()))

    ## Write counts table to file
    if args.ctable:
        vprint(cpatt)
        with open(args.ctable,"w") as fout:
            for line in lineage_counts(sitelist,master,voclist,cpatt,n_sequences):
                print(line,file=fout)

    ## Write x-counts table to file (sequences that don't match patterns)
    ## Include the nearby c-pattern that is closest
    if args.ctable:
        xctable = filename_prepend("x-",args.ctable)
        vprint("Unmatched sequence patterns in file:",xctable)
        with open(xctable,"w") as fout:
            for line in missing_patterns_with_nearby(sitelist,master,voclist,xpatt,n_sequences):
                print(line,file=fout)

    ## now go through the sequences and tally dates

    DG_datecounter = {m: Counter() for m in patterns}
    for s in seqlist[1:]:
        if not args.keepx and "X" in s.seq:
            raise RuntimeError("X's should have already been filtered out")

        seqdate = covid.date_from_seqname(s)
        if not seqdate:
            continue

        vocmatch = svar.vocmatch(s.seq)
        for voc in vocmatch[:1]:
            ## if multiple matches, only count the first
            ## if that's an error, then you'll see a warning from earlier in the computation
            DG_datecounter[voc.flat_pattern][seqdate] += 1
        if not vocmatch:
            DG_datecounter[OTHERNAME][seqdate] += 1

    nmatches = sum(sum(DG_datecounter[p].values()) for p in patterns)
    vprint("matched sequences:",nmatches)
    if nmatches==0:
        raise RuntimeError("No sequences for: " + " ".join(args.filterbyname))

    for p in DG_datecounter:
        vprint(p,sum(DG_datecounter[p].values()))

    ## Onset times for each pattern
    ## Don't include OTHER or patterns that don't appear in sequence set
    onset = {m: min(DG_datecounter[m]) for m in mutants if DG_datecounter[m]}
    for m in onset:
        vprint("onset",m,onset[m])

    ordmin, ordmax, ordplotmin, ordplotmax = get_daterange(DG_datecounter,args.dates)
    Ndays = ordmax+1-ordmin
    vprint("Days:",Ndays,ordmin,ordmax)

    DG_cum=defaultdict(list)
    for ordval in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ordval)
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

    title = covid.get_title(args)
    title = title + ": %d sequences" % (n_sequences,)

    def mmakeplot(fraction=False,lineplot=False):
        embersplot.embersplot(DG_cum,fullnames, mcolors, ordmin,
                              ordplotrange = (ordplotmin, ordplotmax),
                              title=title,
                              legend=args.legend,legendtitle=legendtitle,
                              fraction=fraction,lineplot=lineplot)


    if args.lineplot:
        ## Line plot
        mmakeplot(fraction=True,lineplot=True)
        if args.output:
            plt.savefig(filename_prepend("line-",args.output))

    else:
        ## Stacked bar plots
        outfilename = args.output
        if args.daily == 7:
            outfilename = filename_prepend("wk-",outfilename)
        if args.daily == 1:
            outfilename = filename_prepend("dy-",outfilename)

        mmakeplot(fraction=True)
        if args.output:
            plt.savefig(filename_prepend("f-",outfilename))

        mmakeplot(fraction=False)
        if args.output:
            plt.savefig(filename_prepend("c-",outfilename))

    if not args.output:
        plt.show()

if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very-verbose print'''
        if _args.verbose and _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(_args)
