DESCRIPTION='''
"fix" input fasta (or tbl) sequence file by 
1/ removing columns with dashes in reference sequence,
2/ stripping final stop codons from each sequence,
3/ applying a date filter
'''
import sys
import re
from pathlib import Path
import itertools as it
import random
import argparse

import warnings

import readseq
import sequtil
import intlist
import mutant
import wrapgen
import covid

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("-N",type=int,default=0,
        help="keep at most this many sequences")
    paa("--badisls",
        help="File with list of EPI_ISL numbers to exclude")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def getisls(file):
    '''read list of ISL numbers from text file'''
    isls=[]
    with open(file,"r") as f:
        for line in f:
            m = re.match(r".*(EPI_ISL_\d+).*",line.strip())
            if m:
                isls.append(m[1])
    return isls

def main(args):

    def vcount(seqs,*p,**kw):
        if args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        else:
            return seqs

    seqs = covid.read_seqfile(args)
    seqs = vcount(seqs,"Total sequences read:")
    seqs = covid.filter_seqs(seqs,args)
    seqs = vcount(seqs,"Sequences after filtering:")

    if args.badisls:
        bads = getisls(args.badisls)
        seqs = (s for s in seqs
                   if not any( b in s.name for b in bads ))
        seqs = vcount(seqs,"Sequences after removing bad ISLs:")

    if args.random:
        seqs = list(seqs)
        seqs = seqs[:1] + random.sample(seqs[1:],k=len(seqs[1:]))

    if args.N:
        seqs = it.islice(seqs,args.N)
        seqs = vcount(seqs,"Sequences after truncation:")

    if args.output:
        readseq.write_seqfile(args.output,seqs)
    


if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

