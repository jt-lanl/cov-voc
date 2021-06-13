import sys
import re
from collections import Counter
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

def relativename(master,mutant,dittochar='_'):
    s=""
    for a,b in zip(master,mutant):
        s += dittochar if a==b else b
    return s

def reltoabsname(master,mutant,dittochar='_'):
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

    svar = spikevariants.SpikeVariants()
    if args.colormut:
        svar.init_from_colormut(args.colormut)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar.init_from_defaults()

    ## Would we prefer array of svar's where each svar has mutant, color, name (also a regex_mutant)
    ## Clearly master and sites would be separate arrays

    svar.append_other() ## there's probably a better way!
    
    mutants = svar.mutants 
    master = svar.master
    vocs = svar.vocs
    colors = svar.colors
    sitelist = svar.sites
    namelist = svar.names
    colors[0]  = '#eeeeee' ## very-light gray
    colors[-1] = '#dddddd' ## light gray

    print("master:",master)

    mcolors = {m:c for m,c in zip(mutants,colors)}

    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    mnames =  {m: namefmt%n for m,n in zip(mutants,namelist)}

    def relname(mut):
        return relativename(master,mut,dittochar='_')
    
    rmutants = {mut: relname(mut) for mut in mutants}

    ## mnames and rmutants could also be part of struct

    if args.verbose:
        for m in rmutants:
            vprint(m,rmutants[m])

    print("ok, reading sequences now...",flush=True)
    seqlist = covid.read_seqfile(args)
    seqlist = covid.filter_seqs(seqlist,args)
    print("ok, done reading now...",flush=True)
    svar.checkmaster(seqlist[0].seq) ## ensure master agrees with first seqlist
    for s in seqlist:
        s.seq = "".join(s.seq[n-1] for n in sitelist)
        
    Nsequences = len(seqlist)-1  ## -1 not to count the reference sequence
        
    ## How many of each sequence
    c = Counter(s.seq for s in seqlist)

    ## How many of each mutant
    cpatt = Counter()
    for seq in c:
        for voc in vocs:
            if re.search(voc.re_pattern,seq):
                vvprint(seq,voc.pattern,relname(voc.pattern),c[seq])
                cpatt[voc.pattern] += c[seq]
                #break

    ## make table of patterns and counts
    for line in intlist.write_numbers_vertically(sitelist):
        vprint(line,line)
    for voc in vocs:
        p = voc.pattern
        vprint(p, relname(p), cpatt[p])
    if args.ctable:
        with open(args.ctable,"w") as fout:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line,file=fout)
            print(master,namefmt % "Name","Counts Percent",file=fout)
            for voc in vocs:
                patt = voc.pattern
                if patt == "other":
                    continue
                print("%s %s %6d  %5.2f%%" % (relname(patt),mnames[patt],
                                             cpatt[patt],100*cpatt[patt]/Nsequences),file=fout)

    ## Add s.date and s.mutt attributes to each sequence
    x_count=0
    no_matches=[]
    for s in seqlist:
        if s.name == "master" or s.name == "NC_045512_spike_surface_glycoprotein":
            vprint("Special case:",s.seq,s.name)
            s.date = s.mutt = None
            continue
        if not args.keepx and "X" in s.seq:
            x_count += 1
            s.date = s.mutt = None
            continue
        ## a not-very-robust way to get the date
        tokens = s.name.split('.')
        try:
            s.date = date_fromiso(tokens[-2])
        except IndexError:
            vprint("seq:",s.seq,s.name,"tokens:",tokens)
            s.date = s.mutt = None
            continue
        s.mutt = None
        for voc in vocs:
            patt = voc.pattern
            rpatt = voc.re_pattern
            if re.search(rpatt,s.seq):
                if s.mutt:
                    warnings.warn(f"\n{s.seq} matches\n{relname(s.mutt)} and\n{relname(patt)} with regex\n{rpatt}")
                    continue
                s.mutt = patt
                #break
        if s.mutt is None:
            vvprint("No match",s.seq,s.name)
            no_matches.append(s.seq)
            s.mutt = 'other'
                
    vprint("Sequences with X:",x_count)
    vprint("Sequences without matches:",len(no_matches))
    for line in intlist.write_numbers_vertically(sitelist):
        vprint(line,line)
    seq_nomat = []
    for seq in set(no_matches):
        seq_nomat.append((seq,sum(1 for s in no_matches if s == seq)))
    for seq,nomat in sorted(seq_nomat,key=lambda x: -x[1])[:50]:
        vprint(seq,relname(seq),nomat)

    if args.ctable:
        xctable = re.sub(r"(.*/)?([^/]+)",r"\1x-\2",args.ctable)
        vprint("Missed sequences in file:",xctable)
        with open(xctable,"w") as fout:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line,file=fout)
            print(seqlist[0].seq,"Counts Percent",file=fout)
            for seq,nomat in sorted(seq_nomat,key=lambda x: -x[1])[:50]:
                print("%s %6d  %5.2f%%" % (relname(seq),nomat,100*nomat/Nsequences),file=fout)
                
            
    nmutt = sum(1 for s in seqlist if s.mutt and s.date)
    vprint("   mutt sequences:",nmutt)
    if nmutt==0:
        raise RuntimeError("No sequences for pattern: " + " ".join(args.filterbyname))

    DG_datecounter = {m: Counter() for m in mutants}
    for s in seqlist:
        if s.mutt and s.date:
            DG_datecounter[s.mutt][s.date] += 1
    for p in DG_datecounter:
        vprint(p,sum(DG_datecounter[p].values()))

    onset = {m: min(DG_datecounter[m]) for m in mutants if DG_datecounter[m]}
    for m in onset:
        vprint("onset",m,onset[m])
        
    all_dateset=set()
    for p in DG_datecounter:
        all_dateset.update(DG_datecounter[p])

    vprint("Range of dates:",min(all_dateset),max(all_dateset))
    ordmin = min(all_dateset).toordinal()
    ordmax = max(all_dateset).toordinal()

    ordplotmin = ordmin
    ordplotmax = ordmax
    if args.dates and args.dates[0] and args.dates[0] != ".":
        ordplotmin = date_fromiso(args.dates[0]).toordinal()
    if args.dates and args.dates[1] and args.dates[1] != ".":
        ordplotmax = date_fromiso(args.dates[1]).toordinal()

    ordmin = min([ordmin,ordplotmin])
    ordmax = max([ordmax,ordplotmax])
    
    Ndays = ordmax+1-ordmin
    vprint("Days:",Ndays,ordmin,ordmax)
   
    DG_cum=dict()
    for m in DG_datecounter:
        DG_cum[m] = []
        
    for ord in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord)
        for m in DG_datecounter:
            DG_cum[m].append( sum(DG_datecounter[m][dt] for dt in DG_datecounter[m] if dt <= day ) )

    if args.weekly or args.daily: ## Weekly averages:
        DAYSPERWEEK=7 if args.weekly else 1
        DG_weekly=dict()
        for m in DG_cum:
            DG_weekly[m] = [x for x in DG_cum[m]]
            for n in range(DAYSPERWEEK,Ndays):
                DG_weekly[m][n] = DG_cum[m][n] - DG_cum[m][n-DAYSPERWEEK]
            DG_cum[m] = DG_weekly[m]

    ## Only keep data that is within the specified date range
    ## That way, automatic scaling on the y-axis will be based on available data
    if args.dates:
        for m in mutants:
            ztmp = []
            for i,cnt in enumerate(DG_cum[m]):
                if ordplotmin <= ordmin+i <= ordplotmax:
                    ztmp.append(cnt)
            DG_cum[m] = ztmp
        ordmin = ordplotmin
        Ndays = ordplotmax+1-ordmin

    def stackedbarplot(DG_cum,bigLegend=False,fraction=False):

        mutants = list(DG_cum)
        
        if args.nolegend:
            plt.figure(figsize=(6,3))
        else:
            plt.figure(figsize=(12,1+len(mutants)/4.5))

        title = covid.get_title(args)
        title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')
            
        DG_bottom = dict()
        DG_bottom_current = [0] * len(DG_cum[mutants[0]])
        for m in DG_cum:
            DG_bottom[m] = DG_bottom_current
            DG_bottom_current = [x+y for x,y in zip(DG_bottom_current,DG_cum[m])]

        ## plot a dummy level to get it into the legend as a title
        dummylabel = " "*(maxnamelen+2)
        if bigLegend:
            dummylabel += master
        plt.bar(range(Ndays),[0]*Ndays,width=1,  ## Ndays
                label=dummylabel,color="white")

        #for m in  mutants[::-1]:
        for v in vocs[::-1]:
            m = v.pattern
            mr = rmutants[m] ## rmutants

            ## various hacks
            if mr == master: ## master
                mr = relname(master)
            if mr == "other":
                mr = "." * len(master)
                #mr = "other            "
                #mr = ""

            name = mnames[m] ## mnames
            if bigLegend:
                name += " " + mr
            name = " " + name ## hack! leading underscore doesn't make it to legend??
            if fraction:
                fm = [a/(b+0.001) for a,b in zip(DG_cum[m],DG_bottom_current)]
                bm = [a/(b+0.001) for a,b in zip(DG_bottom[m],DG_bottom_current)]
                plt.bar(range(Ndays),fm,width=1,bottom=bm,
                        label=name, color=mcolors[m])
            else:
                vprint("m,mr,mc:",m,mr,mcolors[m])
                plt.bar(range(Ndays),DG_cum[m],width=1,bottom=DG_bottom[m],
                        label=name, color=mcolors[m])

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
            for m in mutants:
                if m not in onset:
                    continue
                if m == "other":
                    continue
                x = onset[m].toordinal() - ordmin
                kwargs=dict(lw=1,color=mcolors[m])
                if not fraction:
                    kwargs['zorder']=0
                plt.plot([x,x],[ylo,yhi],**kwargs)

        plt.tight_layout()

    def lineplot(DG_cum,fraction=True):
        
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
                mr = relname(master)
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
            plt.savefig("lin-" + args.output)

    else:
        stackedbarplot(DG_cum,fraction=True)
        if args.output:
            plt.savefig("f-wk-"+ args.output)

        stackedbarplot(DG_cum,fraction=False)
        if args.output:
            plt.savefig("c-wk-"+ args.output)
        
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
    

