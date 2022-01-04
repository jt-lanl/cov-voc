'''Find sequences that are far (in Hamming distance) from reference sequence'''
import sys
import itertools as it
from operator import itemgetter
from collections import defaultdict
import argparse

import covid
import sequtil
import mutant
from hamming import hamming

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("-T",type=int,default=10,
        help="show the top T sequences, sorted by Hamming distance")
    paa("-N",type=int,default=0,
        help="only read N sequencesfrom file")
    paa("--ref","-r",
        help="mutant string to be used as reference sequence")
    paa("--nearby",action="store_true",
        help="look for small instead of large Hamming distances")
    paa("--mincount","-m",type=int,default=1,
        help="only show seqs that appear at least m times")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    vprint(args)
    
    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)
    
    if args.N:
        ## just grab the first N (for debugging w/ shorter runs)
        seqs = it.islice(seqs,args.N+1)
    
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)

    if args.ref:
        refmut = mutant.Mutation(args.ref)
        refseq = m_mgr.seq_from_mutation(refmut)
    else:
        refseq = first.seq

    ## deal with duplicates
    seq_cases = defaultdict(list)
    for s in seqs:
        seq_cases[s.seq].append(s)

    hamlist = [ (hamming(refseq,seq), seq_cases[seq]) for seq in seq_cases ]
    hamlist = sorted(hamlist,
                     key=itemgetter(0),
                     reverse=False if args.nearby else True)

    count=0
    for h in hamlist:
        if count >= args.T:
            break
        
        hamdist,seqlist = h
        n = len(seqlist)
        s = seqlist[0]
        if n < args.mincount:
            continue
        
        print("%3d %4d %s" % (hamdist,n,s.name))
        print("        ",m_mgr.get_mutation(s.seq))
        count += 1
        

if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    _main(_args)
