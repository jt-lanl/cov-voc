'''
find sequences in a fasta file that
match a mutant string
'''
import os
import sys
from pathlib import Path
import random
import itertools as it
import argparse

import readseq
import sequtil
import intlist
import covid
import mutant
import wrapgen

def getargs():
    '''get command-line arguments'''
    ap = argparse.ArgumentParser(description=__doc__)
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
    paa("--compact",action="store_true",
        help="write sitelist in a compact way with no room for insertions")
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--nlist",
        help="list of sequences (before filtering); eg. 1-100, or 3-6")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--exact",action="store_true",
        help="require all non-listed sites to match reference sequence")
    paa("--showmutants",action="store_true",
        help="show mutant string after sequence name")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def ndx_and_site_lists(MM,sites,compact=True):
    '''return list of indices and possibly-expanded sitelist'''
    ndxlist = []
    for site in sites:
        if compact:
            ndxlist.append(MM.index_from_site(site))
        else:
            ndxlist.extend(MM.indices_from_site(site))
    sitelist = [MM.site_from_index(ndx) for ndx in ndxlist]
    return ndxlist,sitelist


def read_seqfile(args):
    '''get seqs from file, and filter accordingly'''

    maxseqs = max(intlist.string_to_intlist(args.nlist))+1 \
        if args.nlist else None

    seqs = covid.read_seqfile(args,maxseqs=maxseqs)
    if args.nlist:
        nset = set([0] + intlist.string_to_intlist(args.nlist))
        seqs = (s for n,s in enumerate(seqs) if n in nset)
        seqs = vcount(seqs,"Sequences in nlist:")

    if not args.keepx:
        ## neeed [:-1] here because haven't 'fixed' data yet
        ## via the covid.filter_seqs() call
        seqs = (s for s in seqs if "X" not in s.seq[:-1])
        seqs = vcount(seqs,"Sequences w/o X:")

    seqs = covid.filter_seqs(seqs,args)

    seqs = sequtil.checkseqlengths(seqs)
    if args.random:
        seqlist = list(seqs)
        seqs = seqlist[:1] + random.sample(seqlist[1:],k=len(seqlist[1:]))

    return seqs

def main(args):
    '''main'''

    seqs = read_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs)

    MM = mutant.MutationManager(first.seq)

    mpatt = None
    sites = []

    if args.mutant:
        assert not args.seqpattern
        mpatt = mutant.Mutation.from_mstring(args.mutant)
        sites = sorted(set(ssm.site for ssm in mpatt))

    if args.sites:
        sites = intlist.string_to_intlist(args.sites)
        if args.seqpattern:
            ssms=[]
            for site,mut in zip(sites,args.seqpattern):
                r = MM.refval(site)
                ssms.append( mutant.SingleSiteMutation.from_ref_site_mut(r,site,mut) )
            mpatt = mutant.Mutation(ssms)

    if mpatt:
        seqs = MM.filter_seqs_by_pattern(mpatt,seqs,exact=args.exact)
        if args.verbose:
            seqs = wrapgen.keepcount(seqs,"Sequences matched pattern:")
        seqs = it.chain([first],seqs)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    if args.output:
        seqs = list(seqs)
        readseq.write_seqfile(args.output,seqs)
    #else:
    if 1:
        ndxlist,sitelist = ndx_and_site_lists(MM,sites,compact=args.compact)
        for line in intlist.write_numbers_vertically(sitelist):
            print(line)
        for s in seqs:
            print( "".join(s.seq[n] for n in ndxlist), s.name, end=" ")
            if args.showmutants:
                print(MM.get_mutation(s.seq),end="")
            print()

def _mainwrapper(args):
    '''
    avoids the bulky BrokenPipeError that arises, eg,
    if output is piped through 'head'
    '''
    try:
        main(args)
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        devnull = os.open(os.devnull, os.O_WRONLY)
        os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
    except KeyboardInterrupt:
        print("Keyboard Interrupt",file=sys.stderr)
        sys.exit(1) ## Just exit

if __name__ == "__main__":

    _args = getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''more verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    def vcount(seqs,*p,**kw):
        '''count items in generator as they go by'''
        if _args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        return seqs

    _mainwrapper(_args)
