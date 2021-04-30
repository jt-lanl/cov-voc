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

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        help="input fasta file with all spike sequences")
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("-N",type=int,default=0,
        help="keep at most this many sequences")
    paa("--dates",nargs=2,
        help="range of dates (two dates, yyyy-mm-dd format)")
    paa("--keepdashcols",action="store_true",
        help="Do not strip columns with dash in reference sequence")
    paa("--keeplastchar",action="store_true",
        help="Do not strip final stop codon from end of sequences")
    paa("--fixsiteseventy",action="store_true",
        help="Replace --I and -I- with I-- at sites 68,69,70")
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

def fixsiteseventy(seqs):
    count=0
    for s in seqs:
        if s.seq[67:70] == "--I" or s.seq[67:70] == "-I-":
            count += 1
            s.seq = s.seq[:67] + "I--" + s.seq[70:]
            #print("fix70: ",s.name,file=sys.stderr)
    return count

def filterseqs(args,seqlist):
    ''' pull out a sublist of the sequence list, based on 
    various options in the args structure '''

    firstseq = seqlist[0].seq

    if args.dates:
        seqlist = sequtil.filter_by_date(seqlist,
                                         args.dates[0],args.dates[1],
                                         keepfirst=True)
        vprint(len(seqlist),"sequences in date range:",args.dates)

    if args.badisls:
        bads = getisls(args.badisls)
        seqlist = [s for s in seqlist
                   if not any( b in s.name for b in bads )]
        vprint(len(seqlist),"sequences after removing bad ISLs")

    if "-" in firstseq and not args.keepdashcols:
        vprint("Stripping sites with dashes in first sequence...",end="")
        sequtil.stripdashcols(firstseq,seqlist)
        vprint("ok")
        if "-" in seqlist[0].seq:
            raise RuntimeError("strip dash failed!")

    firstseq = seqlist[0].seq
    
    if "-" in firstseq:
        warnings.warn("dashes in reference sequence")

    return seqlist

def main(args):

    seqlist = readseq.read_seqfile(args.input)
    vprint(len(seqlist),"sequences read")

    seqlist = filterseqs(args,seqlist)
    vprint(len(seqlist),"sequences after filtering")

    if args.random:
        seqlist = seqlist[:1] + random.sample(seqlist[1:],k=len(seqlist[1:]))

    firstseq = seqlist[0].seq
    
    if args.N:
        seqlist = seqlist[:args.N]
        vprint(len(seqlist),"sequences after truncation")

    if args.fixsiteseventy:
        count = fixsiteseventy(seqlist)
        if count:
            vprint(count,"sequences fixed at site 70")

    if not args.keeplastchar and firstseq[-1]=="$":
        for s in seqlist:
            s.seq = s.seq[:-1]

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
    

