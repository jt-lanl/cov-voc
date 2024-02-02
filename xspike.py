'''
XSPIKE: eXplore the SPIKE protein sequence in the SARS CoV-2 virus
'''

import re
from collections import Counter,defaultdict
import numpy as np
import random
import itertools as it
from scipy import stats

import matplotlib.pyplot as plt
import argparse

import warnings

from verbose import verbose as v
import readseq
from spikevariants import SpikeVariants
from hamming import hamming
import sequtil
import intlist
import wrapgen
import covid
import mutant

from xspikeplots import circleplot,heatmap

def getargs():
    ap = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(ap)
    paa = ap.add_argument
    paa("--pairs",action="store_true",
        help="analyze pairwise correlations")
    paa("--allsites",
        help="write fasta file with patterns involving /all/ sites")
    paa("--addsites",
        help="prepend this (comma-separated) list of sites to those being observed")
    paa("--sites",type=int,default=50,
        help="Number of highest-entropy sites to use")
    paa("--baseline",choices=('Wuhan','BA.2','BA.5'), default=None,
        help="Baseline sequences (default is Wuhan)")
    paa("--thresh",type=int,default=2,
        help="Only include patterns that appear at least this many times")
    paa("--entropysamples","-E",type=int,default=0,
        help="For single-site entropy, subsample sequences for faster estimates")
    paa("--cvthresh",type=int,default=3, ## should maybe just be hardcoded
        help="Require this many sequences for each site in the pairwise correlation")
    paa("--plot",action="store_true",
        help="make plots")
    paa("--writeplot",
        help="Write plots to file (instead of showing them on screen)")
    paa("--restrictsites",
        help="Consider only these sites (RBD, NTD, NTD+RBD, or comma-separated list)")
    paa("--nomutlist",action="store_true",
        help="Dont make mutant list at end of lines")
    paa("--colormut",
        help="name of color mutation file (mutation_string,lineage_name) are 2nd,3rd columns")
    paa("--title",
        help="Use this TITLE in plots")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    return args

def stripxs_orig(alist,blist,badchar='X'):
    '''given a pair of lists (or strings), strip element n where
    either alist[n]=='X' or blist[n]=='X'
    return pair of tuples
    '''
    ab = [ (a,b) for a,b in zip(alist,blist)
           if a!=badchar and b!=badchar ]
    alist,blist = zip(*ab)
    return alist,blist

def stripxs(alist,blist,badchar='X'):
    '''given a pair of lists (or strings), strip element n
    from both lists, if either alist[n]=='X' or blist[n]=='X'
    return pair of tuples
    '''
    if badchar in alist or badchar in blist:
        good = [i for i,(a,b) in enumerate(zip(alist,blist))
                if not(a==badchar or b==badchar)]
        alist = [alist[g] for g in good]
        blist = [blist[g] for g in good]
    return alist,blist

def contingency_table(alist,blist,thresh=3):
    '''convert pair of lists (alist,blist) into a contingency table'''
    alist,blist = stripxs(alist,blist)

    acnt = Counter(alist); avals=list(acnt); A=len(acnt)
    bcnt = Counter(blist); bvals=list(bcnt); B=len(bcnt)
    table = np.empty((A,B),dtype=np.int32)

    ab = Counter(zip(alist,blist))
    for i,j in it.product(range(A),range(B)):
        table[i,j] = ab[(avals[i],bvals[j])]

    def minsum():
        '''return the miniumum sum along all the rows and columns'''
        hsum = [np.sum(table[:,j]) for j in range(table.shape[1])]
        vsum = [np.sum(table[i,:]) for i in range(table.shape[0])]
        return min(hsum + vsum)

    if thresh:
        ## Stability (remove rows and columns whose sums are < thresh)
        otable = table.copy()
        while( minsum() < thresh ):
            hok = [np.sum(table[:,j])>=thresh for j in range(table.shape[1])]
            vok = [np.sum(table[i,:])>=thresh for i in range(table.shape[0])]
            table = table[:,hok][vok,:]

    ## Check that table is okay
    hzero = [np.sum(table[:,j])==0 for j in range(table.shape[1])]
    vzero = [np.sum(table[i,:])==0 for i in range(table.shape[0])]
    if (np.any( hzero ) or np.any( vzero )): # or min(table.shape)<2 ):
        ## This should never happen
        v.print("warning: possible issues with contingency table")
        v.print("len a,b:",len(alist),len(blist))
        v.print(acnt)
        v.print(bcnt)
        v.print(ab)
        v.print(hzero)
        v.print(vzero)
        v.print("Table;")
        v.print(table)
        if thresh:
            v.print("Original Table:")
            v.print(otable)

    return table

def cramerv(table):
    k = min(table.shape)
    N = np.sum(table)
    try:
        chisq,_,_,_ = stats.chi2_contingency(table,lambda_='pearson')
    except ValueError:
        chisq=0
        v.print(f"Warning: in chisq k={k}, N={N}")
        v.print("Table:\n",table)

    return np.sqrt(chisq/(N*max([1,(k-1)])))

def mutinfo(table):
    return \
        stats.entropy(np.sum(table,axis=1)) + \
        stats.entropy(np.sum(table,axis=0)) - \
        stats.entropy(table.flatten())


#############################################################

def filename_prepend(pre,file):
    ## prepend a string to a file name; eg
    ## "pre","file" -> "prefile", but also
    ## "pre","dir/file" -> "dir/prefile"
    if not file:
        return file
    return re.sub(r"(.*/)?([^/]+)",r"\1"+pre+r"\2",file)

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
                                      thresh=args.cvthresh)
            cvtable[i,j] = cvtable[j,i] = cramerv(table)
            mitable[i,j] = mitable[j,i] = mutinfo(table)

            v.vvprint()
            v.vvprint("ei,ej:",ei+1,ej+1)
            v.vvprint("cv:",cvtable[i,j])
            v.vvprint("mi:",mitable[i,j])
            v.vvprint("ei:",Counter(charsatsite[ei]))
            v.vvprint("ej:",Counter(charsatsite[ej]))
            v.vvprint("ij:",Counter(zip(charsatsite[ei],
                                      charsatsite[ej])))
            v.vvprint(table)

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
                plt.savefig(filename_prepend(pname + "-" + sname + "_", args.writeplot))

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
    v.print("Running xspike",title)

    allseqs = covid.read_seqfile(args)
    allseqs = vcount(allseqs,"All sequences:")
    allseqs = covid.filter_seqs_by_date(allseqs,args)
    allseqs = vcount(allseqs,"All sequences in date range:")
    allseqs = list(allseqs)
    seqs = covid.filter_seqs_by_pattern(allseqs,args)
    seqs = vcount(seqs,"Sequences after filtering by pattern:")

    seqs = list(seqs)
    firstseq = seqs[0].seq
    seqs = seqs[1:]

    if len(seqs)==0:
        raise RuntimeError("No sequences match pattern")

    N = len(seqs)
    M = len(seqs[0].seq)

    print(f"Evaluating {N} sequences of length {M}")
    print("Sampled from %s to %s."
          % covid.range_of_dates(seqs))
    if args.dates:
        print("Specified Date Range:",args.dates)

    #### SINGLE-SITE ENTROPY
    v.vprint("Single-site entropy...",end="")
    ## If args.entropysamples, then use subsampling of the sequences to estimate entropy
    sampleseqs = random.choices(seqs,k=args.entropysamples) if N>args.entropysamples>0 else seqs
    E = sequtil.chunked_entropy(sampleseqs)
    v.vprint("ok")

    T = mutant.MutationManager(firstseq)

    ## Determine which sites will be employed
    esites = intlist.string_to_intlist(args.addsites)
    esites = list(dict.fromkeys(esites)) ## removes dups while maintaining order
    if args.sites:
        ndxtop = [n for n in np.argsort(E)[::-1]]
        if args.restrictsites:
            rsites = covid.spike_sites(args.restrictsites)
            rindices = T.indices_from_sitelist(rsites)
            ndxtop = [n for n in ndxtop if n in rindices]
        ## convert top indices back into top sites
        etopsites = [T.site_from_index(n) for n in ndxtop]
        ## sites of the largest entropy
        esites.extend(etopsites[:args.sites])
        esites = list(dict.fromkeys(esites)) ## removes dups while maintaining order

    v.vprint("Observing",len(esites),title)
    v.vprint("Observing",len(esites),"sites:",esites)
    v.vprint("Observing",len(esites),"sites:",sorted(esites))

    ## Make entropy plot
    if args.plot or args.writeplot:
        plt.figure()
        ## entropy vs site for all sites ... in default color (blue)
        siterange = [T.site_from_index(n) for n in range(M)]
        plt.plot(siterange,E)
        #plt.plot(range(1,M+1),E)
        plt.xlabel("Site index")
        plt.ylabel("Entropy")
        plt.title(title)
        ## entropy at specified sites ... in red
        for e in esites:
            ndx = T.index_from_site(e)
            plt.plot([e,e],[0,E[ndx]],'r-')

        if args.writeplot:
            plt.savefig(filename_prepend("entrpy-",args.writeplot))

    print("Highest entropy sites for:",title)
    print("  Site Entropy")
    for e in esites:
        n = T.index_from_site(e)
        print("%6d %7.4f" % (e,E[n]))

    #### PAIRWISE ANALYSIS
    if args.pairs:
        charsatsite=dict()
        mutname = dict()
        for e in esites:
            ## don't strip x's quite yet
            n = T.index_from_site(e)
            charsatsite[e] = sequtil.getcolumn(seqs,n,keepx=True)
            mutname[e] = str(e)
        pairwise(args,esites,charsatsite,mutname,title=title)

    #### COMMON PATTERNS

    ## baseline mutation for mstrings
    base_mut = mutant.Mutation(covid.get_baseline_mstring(args.baseline))

    print()
    print("Most common patterns for local area, where Local =",title)
    if not args.nomutlist and args.baseline:
        print(f"Context is relative to baseline sequence = {args.baseline}: {base_mut}")
    esites = sorted(esites)
    for lines in intlist.write_numbers_vertically(esites,plusone=0):
        print(lines)

    pattseqs = sequtil.copy_seqlist(seqs)
    ndxsites = [T.index_from_site(e) for e in esites]
    for s in pattseqs:
        s.seq = "".join(s.seq[n] for n in ndxsites)
    cnt = Counter([s.seq for s in pattseqs])

    ## Do not include patterns with X's in them
    ## set count[patt]=0 if "X" in patt
    ## one-liner: cnt = { patt: cnt[patt] * bool("X" not in patt) for patt in cnt }
    ## alt: cnt = Counter({patt: cnt[patt] for patt in cnt if 'X' not in patt})
    for patt in cnt:
        if 'X' in patt:
            cnt[patt]=0

    patternlist = sorted(cnt, key=cnt.get, reverse=True)
    v.vvprint("Sums:",args.filterbyname,len(seqs),sum(cnt.values()))

    #### Get counts for the various continents
    #### Use all the sequences, not just those that fit pattern
    allseqs = list(allseqs)[1:] ## don't keep the reference sequence

    cont_cnt = dict() ## cnt's for the various continents
    cont_sum = dict()
    Cxcx = covid.parse_continents()
    for cx,c,x in [("Global","Global","")] + Cxcx:
        cseqs = sequtil.filter_by_pattern(allseqs,c)
        if x:
            cseqs = sequtil.filter_by_pattern_exclude(cseqs,x)
        cseqs = list(cseqs)
        #if c == "Global":
        #    print("Global cseqs=",len(cseqs))

        cont_cnt[c] = Counter(sequtil.multicolumn(cseqs,ndxsites,
                                                  keepx=True))
        cont_sum[c] = len(cseqs)
        v.vvprint("Sums:",c,cont_sum[c],sum(cont_cnt[c].values()))

    master =  "".join(firstseq[n] for n in ndxsites)
    
    print(master," Global",
          " ".join("%6s" % covid.ABBREV_CONTINENTS[cx] for cx,_,_ in Cxcx),
          "  Local",
          " Exact  Pct [Context]" if not args.nomutlist else "")

    ## Totals do not include sequences with X at any of the high-entropy sites
    print(" "*len(master),"%7d" % sum(cont_cnt['Global'].values()),
          " ".join(["%6d" % sum(cont_cnt[c].values()) for _,c,_ in Cxcx]),
          "%7d"% sum(cnt.values()),"<----------- Totals" ) ## Totals

    if args.colormut:
        svar = SpikeVariants.from_colormut(args.colormut,refseq=firstseq)
    else:
        svar = SpikeVariants.default(refseq=firstseq)
    def get_lineage(seq):
        vocs = svar.vocmatch(seq)
        return ", ".join(v.name for v in vocs)

    ## dict indexes list of full sequences based on pattseq
    pattseqdict=defaultdict(list)
    for ps,s in zip(pattseqs,seqs):
        pattseqdict[ps.seq].append(s)

    for p in patternlist:

        if "X" in p:
            continue
        if cnt[p] <  args.thresh:
            break
        print("%s %7d " % (sequtil.relativename(master,p),
                               cont_cnt["Global"][p]),end="")
        print(" ".join("%6d"% cont_cnt[c][p] for _,c,_ in Cxcx),end="")
        print(" %7d" % cnt[p],end="")

        ## What is the context for this pattern
        if args.nomutlist:
            print()
        else:
            ## count full seq's consistent with pattern p
            pcnt = Counter(s.seq for s in pattseqdict[p])
            [(pcommonseq,npcs)] = pcnt.most_common(1)
            lineage_name = get_lineage(pcommonseq)
            mut = T.get_mutation(pcommonseq)
            mstring = mut.relative_to(base_mut) if args.baseline else str(mut)
            print(f" {npcs:6d} {100*npcs//cnt[p]:3d}% {mstring} {lineage_name}")



if __name__ == "__main__":

    args = getargs()
    v.verbosity(args.verbose)

    def vcount(seqs,*p,**kw):
        if args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        else:
            return seqs

    main(args)
    if args.plot:
        plt.show()
