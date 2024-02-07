'''
find sequences in a fasta file that
match a mutant string
'''
import os
import sys
import re
from pathlib import Path
import random
import itertools as it
import argparse

import verbose as v
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
    paa("--grep","-g",
        help=("match a string fragment sequence file; "
              "preceed with '/' to avoid leading with '-'"))
    paa("--sites","-s",
        help="list of sites; eg 145-148,156,178-188")
    paa("--seqpattern",
        help="pattern for filtering by sequence")
    paa("--compact",action="store_true",
        help="write sitelist in a compact way with no room for insertions")
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--count",action="store_true",
        help="Only output a count of the number of sequences that match")
    paa("--nlist",
        help="list of sequences (before filtering); eg. 1-100, or 3-6")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--exact",action="store_true",
        help="require all non-listed sites to match reference sequence")
    paa("--uniq",action="store_true",
        help="eliminate sequences with duplicate names")
    paa("--uniqisl",action="store_true",
        help="eliminate sequences with duplicate ISL numbers")
    paa("--showmutants",action="store_true",
        help="show mutant string after sequence name")
    paa("--jobno",type=int,default=1,
        help="job number if using parallel")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    return args


def ndx_and_site_lists(m_mgr,sites,compact=True):
    '''return list of indices and possibly-expanded sitelist'''
    ndxlist = []
    for site in sites:
        if compact:
            ndxlist.append(m_mgr.index_from_site(site))
        else:
            ndxlist.extend(m_mgr.indices_from_site(site))
    sitelist = [m_mgr.site_from_index(ndx) for ndx in ndxlist]
    return ndxlist,sitelist

def read_seqfile(args,firstisref=True):
    '''get seqs from file, and filter accordingly'''

    maxseqs = None
    if args.nlist:
        maxseqs = max(intlist.string_to_intlist(args.nlist))
        if firstisref:
            maxseqs += 1
    seqs = covid.read_seqfile(args,maxseqs=maxseqs)
    if args.nlist:
        nset = set(intlist.string_to_intlist(args.nlist))
        if firstisref:
            nset.add(0)
        seqs = (s for n,s in enumerate(seqs) if n in nset)
        seqs = vcount(seqs,"Sequences in nlist:")

    ## comment: i think covid.filter_seqs also filters if args.keepx is set
    ## but this is here, i guess, so that the vcount will deal with keepx=False
    if not args.keepx:
        ## neeed [:-1] here because haven't 'fixed' data yet
        ## via the covid.filter_seqs() call
        seqs = (s for s in seqs if "X" not in s.seq[:-1])
        seqs = vcount(seqs,"Sequences w/o X:")

    seqs = covid.filter_seqs(seqs,args,firstisref=firstisref)

    #seqs = sequtil.checkseqlengths(seqs)
    if args.random:
        seqlist = list(seqs)
        if firstisref:
            first,seqlist = seqlist[0],seqlist[1:]
        seqlist = random.sample(seqlist,k=len(seqlist))
        if firstisref:
            seqlist = [first]+seqlist
        seqs = seqlist

    return seqs

def keepuniq(seqs):
    '''yield input sequences but with duplicate named sequences filtered out'''
    sname_set = set()
    for s in seqs:
        if s.name in sname_set:
            v.vprint('Dup: ',s.name)
            continue
        sname_set.add(s.name)
        yield s

def keepuniqisl(seqs):
    '''yield input sequences but with duplicate ISL numbers filtered out'''
    isl_set = set()
    for s in seqs:
        isl = covid.get_isl(s.name)
        if isl in isl_set:
            v.vprint('Dup: ',s.name)
            continue
        isl_set.add(isl)
        yield s

def main(args):
    '''main'''

    seqs = read_seqfile(args,firstisref=True)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    m_mgr = mutant.MutationManager(first.seq)

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
                r = m_mgr.refval(site)
                ssm = mutant.SingleSiteMutation.from_ref_site_mut(r,site,mut)
                ssms.append(ssm)
            mpatt = mutant.Mutation(ssms)

    if mpatt:
        seqs = m_mgr.filter_seqs_by_pattern(mpatt,seqs,exact=args.exact)
        if args.verbose:
            seqs = wrapgen.keepcount(seqs,"Sequences matched pattern:")

    if args.grep:
        args.grep = re.sub(r'/','',args.grep) ## you can use /---I--TT/ as a pattern
        seqs = (s for s in seqs
                if args.grep in s.seq)

    if args.uniq:
        seqs = keepuniq(seqs)

    if args.uniqisl:
        seqs = keepuniqisl(seqs)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    ## how many seqs still in the generator after all that filtering/matching
    if args.verbose:
        ## how many sequences match
        seqs = wrapgen.keepcount(seqs,"Sequences match:")
    if args.count:
        ## if --count, then send the count to stdout (w/o msg prefix)
        seqs = wrapgen.keepcount(seqs,None,file=sys.stdout)

    if args.output:
        if args.jobno==1:
            seqs = it.chain([first],seqs)
        readseq.write_seqfile(args.output,seqs)

    if not (args.count or args.output):
        ## write summary to stdout
        ndxlist,sitelist = ndx_and_site_lists(m_mgr,sites,
                                              compact=args.compact)
        if args.jobno == 1:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line)
            seqs = it.chain([first],seqs)
        for s in seqs:
            print( "".join(s.seq[n] for n in ndxlist), s.name, end=" ")
            if args.showmutants:
                print(m_mgr.get_mutation(s.seq),end="")
            print()

    ## Remaining sequences? send them to oblivion...
    ## (we need to do this in order for wrapgen-based count to work)
    for s in seqs:
        pass

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
    v.verbosity(_args.verbose)

    def vcount(seqs,*p,**kw):
        '''count items in generator as they go by'''
        if _args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        return seqs

    _mainwrapper(_args)
