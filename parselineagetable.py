'''
create output text with three columns: sequence name, lineage designation, color
'''

import sys
import itertools as it
import re
import argparse

import warnings

import sequtil
import covid
import colornames

import lineagetable

def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--usehex",action="store_true",
        help="use six-character hex-codes instead of X11 colornames")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def main(args):
    '''main parselineagetable'''

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    first,seqs = sequtil.get_first_item(seqs)
 
    T = lineagetable.get_lineage_table(args.lineagetable)
    
    for s in seqs:
        voc = T.first_match(s.name)
        print(s.name,T.names[voc],T.colors[voc])

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

    main(_args)
