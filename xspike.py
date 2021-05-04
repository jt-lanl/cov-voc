DESCRIPTION='''
XSPIKE: eXplore the SPIKE protein sequence in the SARS CoV-2 virus
'''

import sys
import re
from collections import Counter
import numpy as np
import itertools as it
import scipy.stats as sst

import matplotlib.pyplot as plt
import argparse


import readseq
from hamming import hamming
import sequtil
import intlist
import covid

from xspikeplots import circleplot,heatmap

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    covid.corona_args(ap)
    paa = ap.add_argument
    paa("--allsites",
        help="write fasta file with patterns involving /all/ sites")
    paa("--addsites",
        help="prepend this (comma-separated) list of sites to those being observed")
    paa("--sites",type=int,default=30,
        help="Number of highest-entropy sites to use")
    paa("--keepx",action="store_true",
        help="Keep sequences that have X's in them")
    paa("--cvthresh",type=int,default=3,
        help="Require this many sequences for each site in the pairwise correlation")
    paa("--plot",action="store_true",
        help="make plots")
    paa("--writeplot",
        help="Write plots to file (instead of showing them on screen)")
    paa("--stripdashcols",action="store_true",
        help="Strip dashes from master, and align rest, saving to specified file")
    paa("--restrictsites",
        help="Consider only these sites (RBD, NTD, NTD+RBD, or comma-separated list)")
    paa("--nomutlist",action="store_true",
        help="Dont make mutant list at end of lines")
    paa("--nopairs",action="store_true",
        help="do not look for pairwise correlations")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def xentropy(clist,keepx=False):
    cnt = Counter(clist)
    if not keepx:
        cnt.pop('X',None)
    return sst.entropy(list(cnt.values()))

def stripxs(alist,blist,badchar='X'):
    '''given a pair of lists (or strings), strip element n where 
    either alist[n]=='X' or blist[n]=='X'
    return pair of tuples
    '''
    ab = [ (a,b) for a,b in zip(alist,blist)
           if a!=badchar and b!=badchar ]
    alist,blist = zip(*ab)
    return alist,blist

def contingency_table(alist,blist,keepx=False,thresh=3):
    if not keepx:
        alist,blist = stripxs(alist,blist)
        
    acnt = Counter(alist); avals=list(acnt); A=len(acnt)
    bcnt = Counter(blist); bvals=list(bcnt); B=len(bcnt)
    table = np.empty((A,B),dtype=np.int)

    ab = Counter(zip(alist,blist))
    for i,j in it.product(range(A),range(B)):
        table[i,j] = ab[(avals[i],bvals[j])]

    hzero = [np.sum(table[:,j])==0 for j in range(table.shape[1])]
    vzero = [np.sum(table[i,:])==0 for i in range(table.shape[0])]
    if (np.any( hzero ) or np.any( vzero )):
        ## This should never happen
        print("---------------")
        print("len a,b:",len(alist),len(blist))
        print(acnt)
        print(bcnt)
        print(ab)
        print(hzero)
        print(vzero)
        print("Table;")
        print(table)
        raise RuntimeError("Contingency table has all zeros in a row or column")

    ## Stability term (remove rows and columns whose sums are < thresh)
    if thresh:
        hok = [np.sum(table[:,j])>=thresh for j in range(table.shape[1])]
        vok = [np.sum(table[i,:])>=thresh for i in range(table.shape[0])]
        table = table[:,hok][vok,:]
        
    return table    

def cramerv(table):
    k = min(table.shape)
    N = np.sum(table)
    try:
        chisq,_,_,_ = sst.chi2_contingency(table,lambda_='pearson')
    except ValueError:
        chisq=0
        print("ERROR: in chisq",k,N,file=sys.stderr)
        print(table)
    
    #if alist == blist:
    #    print("alist=blist:",k,N,chisq)
    #    print(table)
    return np.sqrt(chisq/(N*max([1,(k-1)])))

def mutinfo(table):
    return \
        sst.entropy(np.sum(table,axis=1)) + \
        sst.entropy(np.sum(table,axis=0)) - \
        sst.entropy(table.flatten())
    

#############################################################

def get_title(args,title=None):
    if not title:
        title = covid.get_title(args)
        if args.restrictsites:
            title = title + " (" + args.restrictsites + " only)"
    return title

def pairwise(args,esites,charsatsite,mutname,title=None):
    ''' Do pair-wise analysis to get cross-correaltions and mutual entropies '''

    title = get_title(args,title=title)

    ## Make Cramer V table, and Mutual Info table
    ne = len(esites)
    cvtable = np.zeros((ne,ne))
    mitable = np.zeros((ne,ne))
    for i,j in it.product(range(ne),repeat=2):
        ei,ej = esites[i],esites[j]
        if ei <= ej:
            table = contingency_table(charsatsite[ei],charsatsite[ej],
                                      keepx=args.keepx,thresh=args.cvthresh)
            cvtable[i,j] = cvtable[j,i] = cramerv(table)
            mitable[i,j] = mitable[j,i] = mutinfo(table)

            vvprint()
            vvprint("ei,ej:",ei+1,ej+1)
            vvprint("cv:",cvtable[i,j])
            vvprint("mi:",mitable[i,j])
            vvprint("ei:",Counter(charsatsite[ei]))
            vvprint("ej:",Counter(charsatsite[ej]))
            vvprint("ij:",Counter(zip(charsatsite[ei],
                                      charsatsite[ej])))
            vvprint(table)

    ## Make plots of pairwise correlations
    mnames = [mutname[e] for e in esites]
    nodevalues = np.diag(mitable)
    for plotter,pname in ((heatmap,   "hotmap"),
                          (circleplot,"circle")):
        for stat,sname in ((mitable,"mutinfo"),
                            (cvtable,"cramerv")):
            nv = nodevalues if sname == "cramerv" else None
            plotter(stat,mnames,nodevalues=nv)
            plt.title(sname + ": " + title)
            plt.tight_layout()
            if args.writeplot:
                plt.savefig("-".join([pname,sname,args.writeplot]))

    #### Make table of correlated pairs
    print("\nMost highly correlated site-pairs for:",title)
    print("               cramerV  mutInfo")

    lines=[]
    for i,j in it.combinations(range(ne),2):
        lines.append( (mnames[i],mnames[j],cvtable[i,j],mitable[i,j]) )
    lines = sorted(lines,     key=lambda w: (w[3],w[2]),reverse=True)
    for line in lines[:15]:
        print("%6s %6s  %7.4f  %7.4f" % line)
    print()
    lines = sorted(lines[15:],key=lambda w: (w[2],w[3]),reverse=True)
    for line in lines[:15]:
        print("%6s %6s  %7.4f  %7.4f" % line)

def main(args):

    ## Get title for plots and tables
    title = get_title(args)

    print("Running",title,file=sys.stderr)

    allseqs = covid.read_seqfile(args)
    allseqs = covid.filter_seqs_by_date(allseqs,args)
    seqs = covid.filter_seqs_by_pattern(allseqs,args)

    firstseq = seqs[0].seq
    seqs = seqs[1:]
    
    if len(seqs)==0:
        raise RuntimeError("No sequences match pattern")

    N = len(seqs)
    M = len(seqs[0].seq)

    #### SINGLE-SITE ENTROPY
    vprint("Single-site entropy...",end="")
    E = [xentropy(clist,keepx=args.keepx)
         for clist in sequtil.gen_columns_seqlist(seqs)]
    vprint("ok")


    ## Determine which sites will be employed
    esites = intlist.string_to_intlist(args.addsites)
    esites = list(dict.fromkeys(esites)) ## removes dups while maintaining order
    if args.sites:
        etop = [e+1 for e in np.argsort(E)[::-1]]
        if args.restrictsites:
            rsites = covid.spike_sites(args.restrictsites)
            etop = [e for e in etop if e in rsites]
        ## sites of the largest entropy
        esites.extend(etop[:args.sites])
        esites = list(dict.fromkeys(esites)) ## removes dups while maintaining order
        
    vprint("Observing",len(esites),title)
    vprint("Observing",len(esites),"sites:",esites)
    vprint("Observing",len(esites),"sites:",sorted(esites))

    ## Make entropy plot
    if args.plot or args.writeplot:
        plt.figure()
        ## entropy vs site for all sites ... in default color (blue)
        plt.plot(range(1,M+1),E)
        plt.xlabel("Site index")
        plt.ylabel("Entropy")
        plt.title(title)
        ## entropy at specified sites ... in red
        for e in esites:
            plt.plot([e,e],[0,E[e-1]],'r-')
        
        if args.writeplot:
            plt.savefig("entrpy-" + args.writeplot)

    ## get lists of characters for each site, and name of mutation
    charsatsite=dict()
    mutname = dict()
    for e in esites:
        ## don't strip x's quite yet
        charsatsite[e] = sequtil.getcolumn(seqs,e,keepx=True)
        mutname[e] = str(e)

    print("Highest entropy sites for:",title)
    print("  Site Entropy")
    for e in esites:
        print("%6s %7.4f" % (mutname[e],E[e-1]))

    #### PAIRWISE ANALYSIS
    if not args.nopairs:
        pairwise(args,esites,charsatsite,mutname,title=title)


    #### COMMON PATTERNS
    print("\nMost common patterns for:",title)
    esites = np.sort(esites)
    for lines in intlist.write_numbers_vertically(esites,plusone=0):
        print(lines)

    pattseqs = sequtil.copy_seqlist(seqs)
    for s in pattseqs:
        s.seq = "".join(s.seq[n-1] for n in esites)        
    cnt = Counter([s.seq for s in pattseqs])
    patternlist = sorted(cnt, key=cnt.get, reverse=True)

    #### Get counts for the various continents
    #### Use all the sequences, not just those that fit pattern
    allseqs = allseqs[1:] ## don't keep the reference sequence

    cont_cnt = dict() ## cnt's for the various continents
    Cxcx = covid.parse_continents()
    for cx,c,x in [("Global","Global","")] + Cxcx:
        cseqs = sequtil.filter_by_pattern(allseqs,c)
        if x:
            cseqs = sequtil.filter_by_pattern_exclude(cseqs,x)
        cont_cnt[c] = Counter(sequtil.multicolumn(cseqs,esites))
    
    master =  "".join(firstseq[n-1] for n in esites)
    print(master," Local Global",
          " ".join("%6s" % covid.ABBREV_CONTINENTS[cx] for cx,_,_ in Cxcx),
          "[Context]")
    for p in patternlist:
        if args.keepx == False and "X" in p:
            continue
        if cnt[p] <  2: #max([args.cvthresh,min([200,N//5000])]):
            break
        print("%s %6d %6d " % (sequtil.relativename(master,p),cnt[p],
                               cont_cnt["Global"][p]),end="")
        print(" ".join("%6d"% cont_cnt[c][p] for _,c,_ in Cxcx),end="")

        ## What is the context for this pattern
        if args.nomutlist:
            print()
        else:
            pseqnames = set(s.name for s in pattseqs if s.seq == p)
            pfullseqs = [s for s in seqs if s.name in pseqnames]
            pcnt = Counter([s.seq for s in pfullseqs])
            pcommonseq,_ = pcnt.most_common(1)[0]
            mutantstr = sequtil.mutantlist(firstseq,pcommonseq,returnstring=True)
            print("",mutantstr)
        

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose and args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    if args.plot:
        plt.show()
  

