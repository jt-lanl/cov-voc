'''
match a sequence pattern to sequences in a seqfile
ad hoc variant of vfasta.py
requires full spike sequence to EXACTLY match the full sequence pattern
created when Wuhan is mutated as specified at the specified sites
'''
import sys
import re
from collections import Counter
import random
import argparse

import warnings

#import readseq
#import sequtil
import intlist
import mutant
import covid

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("--mutants","-m",
        help="mutant string, such as '[W152R,N439K,D614G,P681R]'")
    paa("--sites","-s",
        help="list of sites; eg 145-148,156,178-188")
    paa("--seqpattern",
        help="pattern for filtering by sequence")
    paa("-N",type=int,default=10,
        help="show at most this many variants")
    paa("--nlist",
        help="list of sequences; eg. 1-100, or 3-6")
    paa("--keepx",action="store_true",
        help="include sequences with X in them")
    paa("--mincount",type=int,default=0,
        help="Require varants to appear in this many sequences")
    #paa("--output","-o",type=Path,
    #    help="output fasta file")
    paa("--fullmatch",action="store_true",
        help="fullpattern = Wuhan + mutations; default: dots + mutations")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def main(args):

    seqlist = covid.read_seqfile(args)
    vprint(len(seqlist),"sequences read")
    seqlist = covid.filter_seqs(seqlist,args)
    vprint(len(seqlist),"sequences after filtering")

    if args.random:
        seqlist = seqlist[:1] + random.sample(seqlist[1:],k=len(seqlist[1:]))

    ## lop off that damn ending '$'
    bottomdollar = re.compile("\$$")
    for s in seqlist:
        s.seq = bottomdollar.sub("",s.seq)

    firstseq = seqlist[0].seq
    
    if args.mutants:
        assert(not args.sites)
        assert(not args.seqpattern)
        muts = mutant.Mutation().init_from_line(args.mutants)
        for mut in muts:
            if mut.ref != firstseq[mut.site-1]:
                vprint("Warning: at site",mut.site,":",
                       mut.ref,"should be",firstseq[mut.site-1])
        seqpattern = "".join(mut.mut for mut in muts)
        sites = [mut.site for mut in muts]
        
    if args.sites:
        assert(args.seqpattern)
        sites = intlist.string_to_intlist(args.sites)
            
    if args.seqpattern:
        assert(args.sites)
        seqpattern = args.seqpattern

    fullpatt = ["."] * len(seqlist[0].seq)
    for c,s in zip(seqpattern,sites):
        fullpatt[s-1] = c
    fullpatt = "".join(fullpatt)

    seqlist = seqlist[:1] + [s for s in seqlist[1:]
                             if re.match(fullpatt,s.seq)]
    vprint("Sequences:",len(seqlist)-1,"match seqpattern:",seqpattern)

    ## Count all the sequences consistent with this sequence
    cnt = Counter(s.seq for s in seqlist)

    ## Don't count sequences with X's in them
    if not args.keepx:
        xlist = [s for s in cnt if "X" in s]
        for x in xlist:
            del cnt[x]

    ## Don't count patterns with fewer than MINCOUNT appearances
    if args.mincount:
        xlist = [s for s in cnt if cnt[s] < args.mincount]
        for x in xlist:
            del cnt[x]
    
    print(len(cnt),"variants of the variant; here are the top",args.N)
    topcnt = sorted(cnt, key=cnt.get, reverse=True)
    for seq in topcnt[:args.N]:
        print(f"{cnt[seq]:6d}",mutant.Mutation().init_from_sequences(seqlist[0].seq,seq))

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    
