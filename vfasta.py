'''
view a segment of a sequence file
'''
import sys
import re
from pathlib import Path
import argparse

import warnings

import readseq
import sequtil
import intlist

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        help="input fasta file with all spike sequences")
    paa("--filtername","-f",
        help="pattern for filtering by name")
    paa("--filterfile",
        help="file with list of patterns for filtering by name")
    paa("--seqpattern",
        help="pattern for filtering by sequence")
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--nlist",
        help="list of sequences; eg. 1-100, or 3-6")
    paa("--sites","-s",
        help="list of sites; eg 145-148,156,178-188")
    paa("--dates",nargs=2,
        help="range of dates (two dates, yyyy-mm-dd format)")
    paa("--stripdashcols",action="store_true",
        help="Strip dash columns")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def main(args):

    seqlist = readseq.read_seqfile(args.input)
    vprint("Sequences:",len(seqlist),"read")

    firstseq = seqlist[0].seq

    if args.dates:
        seqlist = sequtil.filter_by_date(seqlist,
                                         args.dates[0],args.dates[1],
                                         keepfirst=True)
        vprint(len(seqlist),"sequences in date range:",args.dates)                 

    if args.stripdashcols and "-" in firstseq:
        vprint("Stripping sites with dashes in first sequence...",end="")
        sequtil.stripdashcols(firstseq,seqlist)
        vprint("ok")
        if "-" in seqlist[0].seq:
            raise RuntimeError("strip dash failed!")

    firstseq = seqlist[0].seq
    
    if "-" in firstseq:
        warnings.warn("dashes in reference sequence")

    if args.filtername:
        seqlist = seqlist[:1] + [s for s in seqlist[1:] if args.filtername in s.name]
        vprint("Sequences:",len(seqlist),"with name matching:",args.filtername)

    if args.filterfile:
        with open(args.filterfile,'r') as f:
            filternames = [line.strip() for line in f if line.strip()]
            vprint("filternames:")
            for fn in filternames:
                vprint("  ",fn)
            seqlist = seqlist[:1] + [s for s in seqlist[1:]
                                     if any([fname in s.name
                                             for fname in filternames])]
            vprint("Sequences:",len(seqlist),"with name in:",args.filterfile)
            
        

    if args.nlist:
        seqlist = [seqlist[n]
                   for n in intlist.string_to_intlist(args.nlist)]
        
    if args.N:
        seqlist = seqlist[:args.N]
        vprint("Sequences:",len(seqlist),"after truncation")

    if args.sites:
        sites = intlist.string_to_intlist(args.sites)
        for s in seqlist:
            s.seq = "".join(s.seq[n-1] for n in sites)
            
    if args.seqpattern:
        patt = re.compile(args.seqpattern)
        seqlist = seqlist[:1] + [s for s in seqlist[1:]
                                 if patt.match(s.seq)]
        vprint("Sequences:",len(seqlist),"match seqpattern:",args.seqpattern)


    if args.output:
        readseq.write_seqfile(args.output,seqlist)
    else:
        if args.sites:
            for line in intlist.write_numbers_vertically(sites):
                print(line)
        for s in seqlist[:10]:
            print(s.seq)    
    


if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

