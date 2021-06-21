'''
match a sequence pattern to sequences in a seqfile
ad hoc variant of vfasta.py
requires full spike sequence to EXACTLY match the full sequence pattern
created when Wuhan is mutated as specified at the specified sites
'''
import os
import sys
import re
from pathlib import Path
import random
import argparse

import warnings

import readseq
import sequtil
import intlist
import covid
import mutant

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
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--nlist",
        help="list of sequences; eg. 1-100, or 3-6")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--fullmatch",action="store_true",
        help="fullpattern = Wuhan + mutations; default: dots + mutations")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def filterseqs(args,seqlist):
    ''' pull out a sublist of the sequence list, based on 
    various options in the args structure '''

    firstseq = seqlist[0].seq

    if args.nlist:
        seqlist = [seqlist[n]
                   for n in intlist.string_to_intlist(args.nlist)]
        vprint(len(seqlist),"sequences in:",args.nlist)

    if args.dates:
        seqlist = sequtil.filter_by_date(seqlist,
                                         args.dates[0],args.dates[1],
                                         keepfirst=True)
        vprint(len(seqlist),"sequences in date range:",args.dates)

    if args.filtername:
        seqlist = seqlist[:1] + [s for s in seqlist[1:]
                                 if args.filtername in s.name]
        vprint(len(seqlist),"sequences with pattern:",args.filtername)

    if args.stripdashcols and "-" in firstseq:
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
        assert( muts.checkref(firstseq,verbose=True) )
        matchpatt = muts.regex_pattern(firstseq,exact=bool(args.fullmatch))
        sites = [mut.site for mut in muts]
        
    if args.sites:
        assert(args.seqpattern)
        sites = intlist.string_to_intlist(args.sites)
            
    if args.seqpattern:
        assert(args.sites)
        seqpattern = args.seqpattern
        
        matchpatt = list(firstseq) if args.fullmatch else ["."] * len(firstseq)
        for c,s in zip(seqpattern,sites):
            matchpatt[s-1] = c
        matchpatt = "".join(matchpatt)

    rematchpatt = re.compile(matchpatt)
    seqlist = seqlist[:1] + [s for s in seqlist[1:]
                             if rematchpatt.match(s.seq)]
    
    vprint("Sequences:",len(seqlist)-1,"match seqpattern:",matchpatt)

    if args.N:
        seqlist = seqlist[:args.N+1]
        vprint(len(seqlist)-1,"sequences after truncation")

    if args.output:
        readseq.write_seqfile(args.output,seqlist)
    else:
        for line in intlist.write_numbers_vertically(sites):
            print(line)
        for s in seqlist:
            print("".join(s.seq[n-1] for n in sites),s.name)
            ## consider also printing: mutant.mkmutname(seqlist[0].seq,s.seq))

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


    mainwrapper(args)
    

