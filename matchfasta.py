DESCRIPTION='''
find sequences in a fasta file that
match a sequence pattern (or a mutant string) 
'''
import os
import sys
import re
from pathlib import Path
import random
import itertools
import argparse

import warnings

import readseq
import sequtil
import intlist
import covid
import mutant
import wrapgen

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("--mutant","-m",
        help="mutant string, such as '[W152R,N439K,D614G,P681R]'")
    paa("--sites","-s",
        help="list of sites; eg 145-148,156,178-188")
    paa("--seqpattern",
        help="pattern for filtering by sequence")
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--nlist",
        help="list of sequences (before filtering); eg. 1-100, or 3-6")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--fullmatch",action="store_true",
        help="require all non-listed sites to match reference sequence")
    paa("--showmutants",action="store_true",
        help="show mutant string after sequence name")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def main(args):

    maxseqs = max(intlist.string_to_intlist(args.nlist))+1 \
        if args.nlist else None
        
    seqs = covid.read_seqfile(args,maxseqs=maxseqs)
    if args.nlist:
        nlist = [0] + intlist.string_to_intlist(args.nlist)
        seqs = (s for n,s in enumerate(seqs) if n in nlist)
        seqs = vcount(seqs,"Sequences in nlist:")
        
    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)
        seqs = vcount(seqs,"Sequences w/o X:")
        
    seqs = covid.filter_seqs(seqs,args)

    seqs = covid.checkseqlengths(seqs)
    if args.random:
        seqlist = list(seqs)
        seqls = seqlist[:1] + random.sample(seqlist[1:],k=len(seqlist[1:]))

    seqs = iter(seqs)
    first = next(seqs)
    firstseq = first.seq
    
    matchpatt=None
    
    if args.mutant:
        assert(not args.sites)
        assert(not args.seqpattern)
        muts = mutant.Mutation().init_from_line(args.mutant)
        assert( muts.checkref(firstseq,verbose=True) )
        matchpatt = muts.regex_pattern(firstseq,exact=bool(args.fullmatch))
        sites = [mut.site for mut in muts]
        
    if args.sites:
        sites = intlist.string_to_intlist(args.sites)
        if not args.seqpattern:
            args.seqpattern = "." * len(sites)
            
    if args.seqpattern:
        assert(args.sites)
        seqpattern = args.seqpattern
        
        matchpatt = list(firstseq) if args.fullmatch else ["."] * len(firstseq)
        for c,s in zip(seqpattern,sites):
            matchpatt[s-1] = c
        matchpatt = "".join(matchpatt)

    if matchpatt:
        rematchpatt = re.compile(matchpatt)
        seqs = (s for s in seqs if rematchpatt.match(s.seq))
        wrapgen.keepcount(seqs,"Matches:")
        
    #seqlist = list(seqs)
    #vprint("Sequences:",len(seqlist)-1,"match seqpattern:",matchpatt)
    #print("Matches: ",len(seqlist)-1)

    if args.N:
        seqs = itertools.islice(seqs,args.N+1)
        #vprint(len(seqlist)-1,"sequences after truncation")

    if args.output:
        readseq.write_seqfile(args.output,seqs)
    else:
        if matchpatt:
            for line in intlist.write_numbers_vertically(sites):
                print(line)
            for s in seqs:
                print("".join(s.seq[n-1] for n in sites),s.name,end="")
                if args.showmutants:
                    print("",mutant.Mutation((firstseq,s.seq)),end="")
                print()
        else:
            for s in seqs:
                print(s.name)

def mainwrapper(args):
    ''' avoids the bulky BrokenPipeError that arises if, eg,
    output is piped through 'head' '''
    #Traceback (most recent call last):
    #  File "matchfasta.py", line 164, in <module>
    #    main(args)
    #  File "matchfasta.py", line 149, in main
    #    print("".join(s.seq[n-1] for n in sites),s.name)
    #BrokenPipeError: [Errno 32] Broken pipe
    #Exception ignored in: <_io.TextIOWrapper name='<stdout>' mode='w' encoding='UTF-8'>
    #BrokenPipeError: [Errno 32] Broken pipe
    try:
        main(args)
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
            

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)
            
    def vcount(seqs,*p,**kw):
        if args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        else:
            return seqs

    mainwrapper(args)
    

