import sys
import re
import datetime
from collections import Counter
import matplotlib.pyplot as plt
import pickle
import warnings
import argparse

import readseq
import sequtil
import intlist
import spikevariants
import covid

DESCRIPTION='''
Stacked barplots of variant counts over time
'''

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION,
                                 conflict_handler='resolve')
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--keepx",action="store_true",
        help="Keep sequences with X in the pattern")
    paa("--fraction",action="store_true",
        help="Plot fractional instead of numerical values")
    paa("--weekly",action="store_true",
        help="Make weekly average plots instead of cumulative")
    paa("--daily",action="store_true",
        help="Make daily plots instead of cumulative")
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

def main(args):

    seqlist = covid.read_seqfile(args)
    seqlist = covid.filter_seqs(seqlist,args)

    svar = spikevariants.SpikeVariants()
    if args.colormut:
        svar.init_from_colormut(args.colormut,seqlist[0].seq)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar.init_from_defaults()
        
    mutants = svar.mutants
    master = svar.master
    colors = svar.colors
    sitelist = svar.sites
    namelist = svar.names
    colors[0]  = '#eeeeee' ## very-light gray
    colors[-1] = '#dddddd' ## light gray

    if len(mutants) != len(colors):
        print(len(mutants),"mutants")
        print(len(colors),"colors")
        return
    
    mcolors = {m:c for m,c in zip(mutants,colors)}

    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    mnames =  {m: namefmt%n for m,n in zip(mutants,namelist)}

    def relname(mut):
        return relativename(master,mut,dittochar='_')
    
    rmutants = {mut: relname(mut) for mut in mutants}

    if args.verbose:
        for m in rmutants:
            vprint(m,rmutants[m])

    for s in seqlist:
        s.seq = "".join(s.seq[n-1] for n in sitelist)
        
    Nsequences = len(seqlist)-1  ## -1 not to count the reference sequence
        
    ## How many of each sequence
    c = Counter(s.seq for s in seqlist)
    #print(c)

    ## How many of each mutant
    cpatt = Counter()
    for seq in c:
        for patt in mutants:
            if re.search(patt,seq):
                vprint(seq,patt,relname(patt),c[seq])
                cpatt[patt] += c[seq]
                #break

    ## make table to appearances
    for line in intlist.write_numbers_vertically(sitelist):
        vprint(line,line)
    for patt in mutants:
        vprint(patt,relname(patt),cpatt[patt])
    if args.ctable:
        with open(args.ctable,"w") as fout:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line,line,file=fout)
            for patt in mutants:
                if patt == "other":
                    continue
                print(patt,relname(patt),cpatt[patt],file=fout)

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
        tokens = s.name.split('.')
        try:
            s.date = date_fromiso(tokens[-2])
        except IndexError:
            vprint("seq:",s.seq,s.name,"tokens:",tokens)
            s.date = s.mutt = None
            continue
        s.mutt = None
        for patt in mutants:
            if re.search(patt,s.seq):
                if s.mutt:
                    warnings.warn(f"seq: {s.seq} matches {s.mutt} AND {patt}")
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
        seq_nomat.append((seq,len([s for s in no_matches if s == seq])))
    for seq,nomat in sorted(seq_nomat,key=lambda x: -x[1])[:50]:
        vprint(seq,relname(seq),nomat)

    nmutt = len([s for s in seqlist if s.mutt and s.date])
    vprint("   mutt sequences:",nmutt)
    if nmutt==0:
        raise RuntimeError("No sequences for pattern: " + args.filterbyname)

    DG_datelist={m: [] for m in mutants}
    for s in seqlist:
        if s.mutt and s.date:
            DG_datelist[s.mutt].append(s.date)
    for p in DG_datelist:
        vprint(p,len(DG_datelist[p]))

    all_datelist=[]
    for p in DG_datelist:
        all_datelist.extend(DG_datelist[p])

    onset = {m: min(DG_datelist[m]) for m in mutants if DG_datelist[m]}
    for m in onset:
        vprint("onset",m,onset[m])
        
    vprint("Range of dates:",min(all_datelist),max(all_datelist))
    ordmin = min(all_datelist).toordinal()
    ordmax = max(all_datelist).toordinal()

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
    for m in DG_datelist:
        DG_cum[m] = []
    for ord in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord)
        for m in DG_datelist:
            DG_cum[m].append( len([dt for dt in DG_datelist[m] if dt <= day]) )

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
   
    #for m in DG_cum:
    #    plt.plot(DG_cum[m],label=m)
    #plt.legend()

    if args.nolegend:
        plt.figure(figsize=(6,3))
    else:
        plt.figure(figsize=(12,1+len(mutants)/4.5))
    title = covid.get_title(args)
    title = title + ": %d sequences" % (Nsequences,)
    plt.title(title,fontsize='x-large')
    DG_bottom = dict()
    DG_bottom_current = [0] * len(DG_cum[mutants[0]])
    for m in mutants:
        DG_bottom[m] = DG_bottom_current
        DG_bottom_current = [x+y for x,y in zip(DG_bottom_current,DG_cum[m])]

    ## plot a dummy level to get it into the legend as a title
    plt.bar(range(Ndays),[0]*Ndays,width=1,
            label=" "*(maxnamelen+2) + master,color="white")
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

        name = mnames[m] + " " + mr
        name = " " + name ## hack! leading underscore doesn't make it to legend??
        if args.fraction:
            fm = [a/(b+0.001) for a,b in zip(DG_cum[m],DG_bottom_current)]
            bm = [a/(b+0.001) for a,b in zip(DG_bottom[m],DG_bottom_current)]
            plt.bar(range(Ndays),fm,width=1,bottom=bm,
                    label=name, color=mcolors[m])
        else:
            vprint("m,mr,mc:",m,mr,mcolors[m])
            plt.bar(range(Ndays),DG_cum[m],width=1,bottom=DG_bottom[m],
                    label=name, color=mcolors[m])

    if args.fraction:
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
            if not args.fraction:
                kwargs['zorder']=0
            plt.plot([x,x],[ylo,yhi],**kwargs)

    plt.tight_layout()
    if args.output:
        plt.savefig(args.output)
    else:
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
    

