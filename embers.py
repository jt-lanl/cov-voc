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
    paa("--nolegend",action="store_true",help="avoid putting legend on plot")
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
        date = date_fromiso(tokens[-2])
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

def main(args):

    OTHER="other"
    OTHERCOLOR='#dddddd' 

    svar = spikevariants.SpikeVariants()
    if args.colormut:
        svar.init_from_colormut(args.colormut)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar.init_from_defaults()

    #svar.append_other() ## there's probably a better way!
    #and "other" should be a variable OTHER that's maybe all spaces or something
    
    master = svar.master
    sitelist = svar.sites
    vocs = svar.vocs

    mutants = [v.pattern for v in vocs]
    patterns = mutants + [OTHER]
    
    colors = [v.color for v in vocs] + [OTHERCOLOR]
    mcolors = {m:c for m,c in zip(patterns,colors)}

    namelist = [v.name for v in vocs] + [OTHER]
    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    mnames =  {m: namefmt%n for m,n in zip(patterns,namelist)}
    
    #colors[-1] = '#dddddd' ## light gray

    def relpattern(mut):
        if mut == OTHER:
            return "." * len(master)
        return relativepattern(master,mut,dittochar='_')
    
    rmutants = {p: relpattern(p) for p in patterns}

    vprint("ok, reading sequences now...",end="")
    seqlist = covid.read_seqfile(args)
    seqlist = covid.filter_seqs(seqlist,args)
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

        vocmatch = [voc for voc in vocs if voc.re_pattern.match(seq)]
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
            for line in intlist.write_numbers_vertically(sitelist):
                print(line,file=fout)
            print(master,namefmt % "Name","Counts Percent",file=fout)
            for voc in vocs:
                patt = voc.pattern
                if patt == "other":
                    continue
                print("%s %s %6d  %5.2f%%" % (relpattern(patt),mnames[patt],
                                             cpatt[patt],100*cpatt[patt]/Nsequences),file=fout)

    ## Write x-counts table to file (sequences that don't match patterns)
    if args.ctable:
        xctable = filename_prepend("x-",args.ctable)
        vprint("Unmatched sequences in file:",xctable)
        with open(xctable,"w") as fout:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line,file=fout)
            print(master,"Counts Percent",file=fout)
            for seq in sorted(xpatt,key=xpatt.get,reverse=True)[:50]:
                print("%s %6d %6.2f%%" %
                      (relpattern(seq),xpatt[seq],100*xpatt[seq]/Nsequences),file=fout)

                
    DG_datecounter = {m: Counter() for m in patterns}
    DG_other = Counter()
    for s in seqlist[1:]:
        if not args.keepx and "X" in s.seq:
            raise RuntimeError("X's should have already been filtered out")

        s.date = date_from_seqname(s.name)
        if not s.date:
            continue

        vocmatch = [voc for voc in vocs if voc.re_pattern.match(s.seq)]
        for voc in vocmatch:
            DG_datecounter[voc.pattern][s.date] += 1
        if not vocmatch:
            DG_other[s.date] += 1
            #DG_datecounter['other'][s.date] += 1
            
    nmatches = sum(sum(DG_datecounter[p].values()) for p in patterns)
    vprint("matched sequences:",nmatches)
    if nmatches==0:
        raise RuntimeError("No sequences for pattern: " + " ".join(args.filterbyname))

    for p in DG_datecounter:
        vprint(p,sum(DG_datecounter[p].values()))

    onset = {m: min(DG_datecounter[m]) for m in mutants if DG_datecounter[m]}
    for m in onset:
        vprint("onset",m,onset[m])

    DG_datecounter['other'] = DG_other

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

    def stackedbarplot(DG_cum,bigLegend=False,fraction=False):

        patterns = list(DG_cum)
        
        if args.nolegend:
            plt.figure(figsize=(6,3))
        else:
            plt.figure(figsize=(12,1+len(patterns)/4.5))

        title = covid.get_title(args)
        title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')
            
        DG_bottom = dict()
        DG_bottom_current = [0] * len(DG_cum[patterns[0]])
        for m in DG_cum:
            DG_bottom[m] = DG_bottom_current
            DG_bottom_current = [x+y for x,y in zip(DG_bottom_current,DG_cum[m])]

        ## plot a dummy level to get it into the legend as a title
        dummylabel = " "*(maxnamelen+2)
        if bigLegend:
            dummylabel += master
        plt.bar(range(Ndays),[0]*Ndays,width=1,  ## Ndays
                label=dummylabel,color="white")

        name_color_sofar = set()
        for m in  patterns[::-1]:
            if m == OTHER:
                name = OTHER
                mr = "." * len(master)
            else:
                name = mnames[m] ## mnames
                mr = rmutants[m] ## rmutants

            ## various hacks
            if mr == master: ## master
                mr = relpattern(master)

            if bigLegend:
                name += " " + mr
            name = " " + name ## hack! leading underscore doesn't make it to legend??

            name_color = (name,mcolors[m])
            if name_color in name_color_sofar:
                kwargs = dict(color=mcolors[m])
            else:
                name_color_sofar.add(name_color)
                kwargs = dict(label=name,color=mcolors[m])
            
            
            if fraction:
                fm = [a/(b+0.001) for a,b in zip(DG_cum[m],DG_bottom_current)]
                bm = [a/(b+0.001) for a,b in zip(DG_bottom[m],DG_bottom_current)]
                plt.bar(range(Ndays),fm,width=1,bottom=bm,**kwargs)
            else:
                vprint("m,mr,mc:",m,mr,mcolors[m])
                plt.bar(range(Ndays),DG_cum[m],width=1,bottom=DG_bottom[m],**kwargs)

        if fraction:
            plt.ylim([0,1.05])

        if not args.nolegend:
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

        plt.tight_layout()

    def lineplot(DG_cum,fraction=True):

        mutants = list(DG_cum)
        
        if args.nolegend:
            plt.figure(figsize=(6,3))
        else:
            plt.figure(figsize=(8,max([3,len(mutants)/4])))
            #plt.figure(figsize=(12,1+len(mutants)/4.5))
        title = covid.get_title(args)
        title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')
        DG_bottom = dict()
        DG_bottom_current = [0] * len(DG_cum[mutants[0]])
        for m in mutants:
            DG_bottom[m] = DG_bottom_current
            DG_bottom_current = [x+y for x,y in zip(DG_bottom_current,DG_cum[m])]

        ## plot a dummy level to get it into the legend as a title
        #plt.bar(range(Ndays),[0]*Ndays,width=1,
        #        label=" "*(maxnamelen+2) + master,color="white")
        for m in  mutants[::-1]:
            #print(m,"Ndays=",Ndays,
            #      "len(DG_cum[m])=",len(DG_cum[m]),len(DG_bottom[m]))
            mr = rmutants[m]

            ## various hacks
            if mr == master:
                mr = relpattern(master)
            if mr == "other":
                mr = "." * len(master)
                #mr = 'SLVYALKNLYESNATQP', #G beige
                #mr = "......other......"
                #mr = "other            "
                #mr = ""

            name = mnames[m]
            #name = name + " " + mr
            name = " " + name ## hack! leading underscore doesn't make it to legend??
            plotkw = dict(lw=2, label=name, color=mcolors[m])
            if fraction:
                fm = [a/(b+0.001) for a,b in zip(DG_cum[m],DG_bottom_current)]
                bm = [a/(b+0.001) for a,b in zip(DG_bottom[m],DG_bottom_current)]
                dy = np.array(range(Ndays))
                fm = np.array(fm)
                dy = dy[fm>0]
                fm = fm[fm>0]
                plt.semilogy(dy,fm, **plotkw)
            else:
                vprint("m,mr,mc:",m,mr,mcolors[m])
                plt.semilogy(range(Ndays),DG_cum[m], **plotkw)

        #if args.fraction:
        #    plt.ylim([0,1.05])

        if not args.nolegend:
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

        if args.onsets:
            ylo,yhi = plt.gca().get_ylim()
            for m in mutants:
                if m not in onset:
                    continue
                if m == "other":
                    continue
                x = onset[m].toordinal() - ordmin
                if x == 0:
                    continue
                kwargs=dict(lw=1,color=mcolors[m])
                if not fraction:
                    kwargs['zorder']=0
                plt.plot([x,x],[ylo,1],**kwargs)
                plt.text(x,1 - 0.03*math.log(ylo),mnames[m],rotation='vertical',ha='center',size=6)
                plt.ylim(bottom = ylo,
                         top = (1/ylo)**(0.3))

        if args.nolegend:
            plt.subplots_adjust(bottom=0.15,right=0.95) ## hardcoded hack!
        else:
            #plt.subplots_adjust(bottom=0.15,right=0.8) ## still hardcoded!
            plt.subplots_adjust(bottom=0.15,right=0.7) ## still hardcoded!

        def ytickformat(x,pos):
            vvprint("x,pos=",x,pos)
            if x>1:
                return ""
            else:
                return "%.1g" % x

        #plt.gca().yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%.1g'))
        plt.gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(ytickformat))

    if args.lineplot:
        lineplot(DG_cum)
        if args.output:
            plt.savefig(filename_prepend("line-",args.output))

    else:
        outfilename = args.output
        if args.weekly:
            outfilename = filename_prepend("wk-",outfilename)
        if args.daily:
            outfilename = filename_prepend("dy-",outfilename)
            
        stackedbarplot(DG_cum,fraction=True)
        if args.output:
            plt.savefig(filename_prepend("f-",outfilename))

        stackedbarplot(DG_cum,fraction=False,bigLegend=False)
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
    

