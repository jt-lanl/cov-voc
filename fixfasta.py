DESCRIPTION='''
"fix" input fasta (or tbl) sequence file by 
1/ removing columns with dashes in reference sequence,
2/ stripping final stop codons from each sequence,
3/ applying a date filter
'''
import sys
import re
from pathlib import Path
import random
import argparse

import warnings

import readseq
import sequtil
import intlist
import mutant
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

    seqlist = covid.read_seqfile(args)
    seqlist = covid.filter_seqs(seqlist,args)
    vprint(len(seqlist),"sequences after filtering")

    if args.badisls:
        bads = getisls(args.badisls)
        seqlist = [s for s in seqlist
                   if not any( b in s.name for b in bads )]
        vprint(len(seqlist),"sequences after removing bad ISLs")


    if args.random:
        seqlist = seqlist[:1] + random.sample(seqlist[1:],k=len(seqlist[1:]))

    firstseq = seqlist[0].seq
    
    if args.N:
        seqlist = seqlist[:args.N]
        vprint(len(seqlist),"sequences after truncation")


    if args.output:
        readseq.write_seqfile(args.output,seqlist)
    


if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

