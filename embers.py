import sys
import re
from collections import Counter,defaultdict
import datetime

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import warnings
import argparse

import readseq
import sequtil
import intlist
import spikevariants
import covid

DESCRIPTION='''
Stacked barplots (also, optionally, line-plots) of variant counts over time
'''

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION,
                                 conflict_handler='resolve')
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--keepx",action="store_true",
        help="Keep sequences with X in the pattern")
    #paa("--fraction",action="store_true",
    #    help="Plot fractional instead of numerical values")
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
    paa("--colormut",
        help="read SpikeVariants structure from color-mut file")
    paa("--ctable",
        help="write a count table to this file")
    paa("--output","-o",help="write plot to file")
    paa("--verbose","-v",action="count",
        help="verbosity")
    args = ap.parse_args()
    if (args.daily and args.weekly):
        raise RuntimeError("Pick only one of daily or weekly")
    return args

def half_labels(labels,n=2):
    '''replace list of labels ["one","two","three","four","five"]
    with ["one", "", "three", "", "five"]
    '''
    hlabels=[]
    for i,label in enumerate(labels):
        hlabels.append( label if i%n==0 else "" )
    return hlabels


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

def date_friendly(dt):
    #return dt.strftime("%b %-d")
    mmm = dt.strftime("%b")
    dd = dt.strftime("%-d")
    return "%3s %2d" % (mmm, int(dd))

def date_from_seqname(sname):
    ## a not-very-robust way to get the date
    ## alternative would be to search for /d/d/d/d-/d/d-/d/d
    tokens = sname.split('.')
    try:
        date = date_fromiso(tokens[4])
    except IndexError:
        warnings.warn(f"No date found for: {sname}, tokens={tokens}")
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

def filename_prepend(pre,file):
    ## prepend a string to a file name; eg
    ## "pre","file" -> "prefile", but also
    ## "pre","dir/file" -> "dir/prefile"
    if not file:
        return file
    return re.sub(r"(.*/)?([^/]+)",r"\1"+pre+r"\2",file)

def relativepattern(master,mutant,dittochar='_'):
    s=""
    for a,b in zip(master,mutant):
        s += dittochar if a==b else b
    return s

def reltoabspattern(master,mutant,dittochar='_'):
    s=""
    for a,b in zip(master,mutant):
        s += a if b==dittochar else b
    return s

def get_regex_mutant(master,mutant):
    raise RuntimeError("dont use this")
    r=""
    for a,b in zip(master,mutant):
        if b == "!":
            r += "[^"+a+"]"
        else:
            r += b
    return r

def check_dups(xlist):
    '''see if there are any duplicates in the list xlist'''
    duplist=[]
    xset=set()
    for x in xlist:
        if x in xset:
            duplist.append(x)
        xset.add(x)
    return duplist

def lineage_counts(sitelist,master,voclist,cpatt,Nsequences):
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent Lineage"
    for voc in voclist[::-1]:
        patt = voc.pattern
        rpatt = relativepattern(master,patt)
        yield "%s %6d %6.2f%% %s" % (rpatt,cpatt[patt],100*cpatt[patt]/Nsequences,voc.name)

def missing_patterns_with_nearby(sitelist,master,voclist,xpatt,Nsequences):
    ''' for each of the seq patterns in xpatt, find nearby patterns in voclist '''
    ## routine yields lines that are meant to be printed
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent NearbyLineage(s)"
    for seq in sorted(xpatt,key=xpatt.get,reverse=True)[:50]:
        rseq = relativepattern(master,seq)
        dv = {v: sum(bool(x != y and y != ".")
                     for x,y in zip(seq,v.pattern))
              for v in voclist}
        vnearby_names = [v.name.strip() for v in dv if dv[v]<2]
        yield "%s %6d %6.2f%% %s" % (rseq,xpatt[seq],100*xpatt[seq]/Nsequences,
                             ", ".join(vnearby_names))

def main(args):

    OTHER="other"
    OTHERCOLOR='#dddddd' 

    svar = spikevariants.SpikeVariants()
    if args.colormut:
        svar.init_from_colormut(args.colormut)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar.init_from_defaults()

    master = svar.master
    sitelist = svar.sites
    voclist = svar.vocs

    mutants = [v.pattern for v in voclist]
    patterns = mutants + [OTHER]
    dups = check_dups(patterns)
    if dups:
        raise RuntimeError(f"Duplicated patterns {dups}")
    
    colors = [v.color for v in voclist] + [OTHERCOLOR]
    mcolors = {m:c for m,c in zip(patterns,colors)}
    dups = check_dups(colors)
    if dups:
        vprint("Duplicated colors:",dups)

    namelist = [v.name for v in voclist] + [OTHER]
    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    for v in voclist:
        v.name = namefmt % v.name
    mnames =  {m: namefmt%n for m,n in zip(patterns,namelist)}
    dups = check_dups(namelist)
    if dups:
        vprint("Duplicated names:",dups)
    
    def relpattern(mut):
        if mut == OTHER:
            return "." * len(master)
        return relativepattern(master,mut,dittochar='_')
    
    mrelpatt = {p: relpattern(p) for p in patterns}

    for p in patterns:
        vprint(mnames[p],mrelpatt[p],mcolors[p])

    ## at this point, we have mutants, patterns, mcolors, mnames, mrelpatt, sitelist

    vprint("ok, reading sequences now...",end="")
    seqlist = covid.read_seqfile(args)
    seqlist = covid.filter_seqs(seqlist,args)
    seqlist = list(seqlist)
    svar.checkmaster(seqlist[0].seq) ## ensure master agrees with first seqlist
    for s in seqlist:
        s.seq = "".join(s.seq[n-1] for n in sitelist)

    Nsequences = len(seqlist)-1  ## -1 not to count the reference sequence

    if not args.keepx:
        seqlist = [s for s in seqlist if "X" not in s.seq]
        vprint("Removed",Nsequences+1-len(seqlist),"sequences with X")
        Nsequences = len(seqlist)-1        
        
    ## How many of each sequence
    c = Counter(s.seq for s in seqlist[1:])

    ## How many of each mutant
    cpatt = Counter()
    xpatt = Counter()
    for seq in c:

        vocmatch = [voc for voc in voclist if voc.re_pattern.match(seq)]
        for voc in vocmatch:
            cpatt[voc.pattern] += c[seq]

        ## Ideally just one match, if zero or more than one, then...
        if len(vocmatch)==0:
            xpatt[seq] = c[seq]
        elif len(vocmatch)>1:
            warn_msg = f"\n{seq} (count={c[seq]}) matches\n"
            warn_msg += " and\n".join(f"{relpattern(v.pattern)}" for v in vocmatch)
            warnings.warn(warn_msg)

            
    vprint("Unmatched sequences:",sum(xpatt.values()))

    ## Write counts table to file
    if args.ctable:
        with open(args.ctable,"w") as fout:
            for line in lineage_counts(sitelist,master,voclist,cpatt,Nsequences):
                print(line,file=fout)
                
    ## Write x-counts table to file (sequences that don't match patterns)
    ## Include the nearby c-pattern that is closest
    if args.ctable:
        xctable = filename_prepend("x-",args.ctable)
        vprint("Unmatched sequence patterns in file:",xctable)
        with open(xctable,"w") as fout:
            for line in missing_patterns_with_nearby(sitelist,master,voclist,xpatt,Nsequences):
                print(line,file=fout)

    ## now go through the sequences and tally dates
    
    DG_datecounter = {m: Counter() for m in patterns} 
    for s in seqlist[1:]:
        if not args.keepx and "X" in s.seq:
            raise RuntimeError("X's should have already been filtered out")

        seqdate = date_from_seqname(s.name)
        if not seqdate:
            continue

        vocmatch = [voc for voc in voclist if voc.re_pattern.match(s.seq)]
        for voc in vocmatch:
            DG_datecounter[voc.pattern][seqdate] += 1
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
    onset = {m: min(DG_datecounter[m]) for m in mutants if DG_datecounter[m]}
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

    def makeplot(legend=0,fraction=False,linePlot=False):

        barPlot = not linePlot

        if legend == 0:
            plt.figure(figsize=(6,3))
        elif legend == 1:
            F = 4.7 if linePlot else 5.0
            plt.figure(figsize=(12,1+len(patterns)/F))
        elif legend == 2:
            F = 4.7 if linePlot else 4.5
            plt.figure(figsize=(12,1+len(patterns)/F))
        else:
            raise RuntimeError(f"legend should be 0, 1, or 2: not {legend}")

        title = covid.get_title(args)
        title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')
            
        DG_bottom = dict()
        DG_bottom_current = [0] * len(DG_cum[patterns[0]])
        for m in DG_cum:
            DG_bottom[m] = DG_bottom_current
            DG_bottom_current = [x+y for x,y in zip(DG_bottom_current,DG_cum[m])]

        ## plot a dummy level to get it into the legend as a title
        if legend > 1 and barPlot:
            dummylabel = " "*(maxnamelen+2)
            dummylabel += master
            plt.bar(range(Ndays),[0]*Ndays,width=1,
                    label=dummylabel,color="white")

        name_color_sofar = set()
        for m in  patterns[::-1]:
            
            name = mnames[m] ## mnames
            if legend>1:
                name += " " + mrelpatt[m] 
            name = " " + name ## hack! leading underscore doesn't make it to legend??

            if linePlot:
                kwargs = dict(color=mcolors[m],label=name)
            else:
                ## For bar plot, repeated name,color tuples only appear once in legend
                name_color = (name,mcolors[m])
                kwargs = dict(color=mcolors[m])
                if name_color not in name_color_sofar:
                    name_color_sofar.add(name_color)
                    kwargs['label']=name            
            
            if fraction:
                fm = [a/(b+0.001) for a,b in zip(DG_cum[m],DG_bottom_current)]
                bm = [a/(b+0.001) for a,b in zip(DG_bottom[m],DG_bottom_current)]
            else:
                fm = DG_cum[m]
                bm = DG_bottom[m]

            if linePlot:
                dy = np.array(range(Ndays))
                fm = np.array(fm)
                dy = dy[fm>0]
                fm = fm[fm>0]
                plt.semilogy(dy,fm,lw=2, **kwargs)
            else:
                plt.bar(range(Ndays),fm,bottom=bm,width=1,**kwargs)

        if fraction and not linePlot:
            plt.ylim([0,1.05])

        if legend:
            plt.legend(bbox_to_anchor=(1.02, 1),
                       #handlelength=3,
                       #markerfirst=False,
                       frameon=False,
                       handletextpad=0,
                       labelspacing=0.45,
                       loc='upper left', borderaxespad=0.,
                       prop={'family' : 'monospace'})

        plt.xlim(ordplotmin-ordmin-1,ordplotmax-ordmin+1)
        xticks = list(range(ordplotmin-ordmin,ordplotmax-ordmin+1,7)) ## was n+6
        xlabels = [datetime.date.fromordinal(int(ord+ordmin)) for ord in xticks]
        xlabels = [date_friendly(dt) for dt in xlabels]
        if len(xlabels) > 16:
            xlabels = half_labels(xlabels)
        plt.xticks(xticks,xlabels,fontsize='small',
                   rotation=45,ha='right',position=(0,0.01))
        #plt.xlabel("Date (2020)")

        if args.onsets:
            ylo,yhi = plt.gca().get_ylim()
            for m in patterns:
                if m not in onset:
                    continue
                if m == OTHER:
                    continue
                x = onset[m].toordinal() - ordmin
                kwargs=dict(lw=1,color=mcolors[m])
                if not fraction:
                    kwargs['zorder']=0
                plt.plot([x,x],[ylo,yhi],**kwargs)

        if linePlot:
            if legend==0:
                plt.subplots_adjust(bottom=0.15,right=0.95) ## hardcoded hack!
            else:
                #plt.subplots_adjust(bottom=0.15,right=0.8) ## still hardcoded!
                #plt.subplots_adjust(bottom=0.15,right=0.7) ## still hardcoded!
                pass

            def ytickformat(x,pos):
                if x>1:
                    return "" if fraction else "%g" % x
                else:
                    return "%.1g" % x

            plt.gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(ytickformat))

        
        plt.tight_layout()

    if args.lineplot:
        ## Line plot
        makeplot(legend=args.legend,fraction=True,linePlot=True)
        if args.output:
            plt.savefig(filename_prepend("line-",args.output))

    else:
        ## Stacked bar plots 
        outfilename = args.output
        if args.weekly:
            outfilename = filename_prepend("wk-",outfilename)
        if args.daily:
            outfilename = filename_prepend("dy-",outfilename)
            
        makeplot(legend=args.legend,fraction=True)
        if args.output:
            plt.savefig(filename_prepend("f-",outfilename))

        makeplot(legend=args.legend,fraction=False)
        if args.output:
            plt.savefig(filename_prepend("c-",outfilename))
        
    if not args.output:
        plt.show()

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose and args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

