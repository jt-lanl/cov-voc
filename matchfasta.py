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
import breakpipe

import readseq
import sequtil
import intlist
import covid
import mutant
import wrapgen

def getargs():
    '''get command-line arguments'''
    ap = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(ap)
    paa = ap.add_argument_group("Matching Options").add_argument
    paa("--mutant","-m",
        help="mutant string, such as '[W152R,N439K,D614G,P681R]'")
    paa("--extended","-e",action="store_true",
        help="allow mutant string to include wildcards (slower)")
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
        help="show at most this many matched sequences")
    paa("--exact",action="store_true",
        help="require all non-listed sites to match reference sequence")
    paa("--uniq",action="store_true",
        help="eliminate sequences with duplicate names")
    paa("--uniqisl",action="store_true",
        help="eliminate sequences with duplicate ISL numbers")
    paa("--showmuts",action="store_true",
        help="show mutant string after sequence name")
    paa("--jobno",type=int,default=1,
        help="job number if using parallel")
    paa("--count",action="store_true",
        help="Print count of the number of sequences that match (no other output!)")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    if args.count and args.output:
        v.print("Warning: No output if '--count' is invoked")
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

## should this maybe be in mutant.py ??
def mut_pattern_match(m_mgr,mpatt,seqs,exact=False):
    '''alternative to (s for s in seqs if mutant.seq_fits_pattern(mpatt,s.seq))'''
    ## re-express the mpatt as an explicit full-length sequence
    mpatt_seq = m_mgr.seq_from_mutation(mpatt)
    if exact:
        for s in seqs:
            if s.seq == mpatt_seq:
                yield s
    elif 0:
        ## obtain indexlist from mutation sites
        ndxlist = []
        for ssm in mpatt:
            ndxlist.extend( m_mgr.indices_from_site(ssm.site) )
        for s in seqs:
            if all(s.seq[ndx] == mpatt_seq[ndx]
                   for ndx in ndxlist):
                yield s
    else:
        ## convert list of indices to list of ranges
        ndxlist = []
        for ssm in mpatt:
            ndxs = m_mgr.indices_from_site(ssm.site)
            if ssm.ref == "+" and len(ndxs)>1:
                ndxs = ndxs[1:] ## if leading with "+" don't need site iteslf
            ndxlist.extend( ndxs )
        ndxrangelist = intlist.intlist_to_rangelist(ndxlist)
        ## pre-compute pattern for each range
        mpatt_seqlist = [mpatt_seq[lo:hi] for lo,hi in ndxrangelist]
        for s in seqs:
            ## see if s.seq matches mpatt_seq over ALL the ranges
            ## equiv: fails to match for ANY range
            ismatch = True
            for (lo,hi),patt in zip(ndxrangelist,mpatt_seqlist):
                if s.seq[lo:hi] != patt:
                    ismatch = False
                    break
            if ismatch:
                yield s

@breakpipe.no_broken_pipe
def main(args):
    '''main'''

    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,keepfirst=False)

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
                ref = m_mgr.refval(site)
                ssm = mutant.SingleSiteMutation(ref,site,mut)
                ssms.append(ssm)
            mpatt = mutant.Mutation(ssms)

    if mpatt and args.extended:
        seqs = m_mgr.filter_seqs_by_pattern(mpatt,seqs,exact=args.exact)
        if args.verbose:
            seqs = wrapgen.keepcount(seqs,"Sequences matched pattern:")

    if mpatt and not args.extended:
        seqs = mut_pattern_match(m_mgr,mpatt,seqs,exact=args.exact)

    if args.grep:
        args.grep = re.sub(r'/','',args.grep) ## you can use /---I--TT/ as a pattern
        greprange = slice(None,None)
        if sites:
            ndxlo = m_mgr.index_from_site(min(sites))
            ndxhi = m_mgr.index_from_site(max(sites)+1)
            greprange = slice(ndxlo,ndxhi)
        seqs = (s for s in seqs if args.grep in s.seq[greprange])

    if args.uniq:
        seqs = keepuniq(seqs)

    if args.uniqisl:
        seqs = keepuniqisl(seqs)

    if args.N:
        seqs = it.islice(seqs,args.N)

    ## how many seqs still in the generator after all that filtering/matching
    seqs = vcount(seqs,"Sequences match:")

    ## Only one of the three options below (either count, or output, or neither)
    ## All three of thse options consume the generator "seqs"

    if args.count:
        print( sum(1 for s in seqs) )

    elif args.output:
        if args.jobno==1:
            seqs = it.chain([first],seqs)
        sequtil.write_seqfile(args.output,seqs)

    else:
        ## write summary to stdout
        ndxlist,sitelist = ndx_and_site_lists(m_mgr,sites,
                                              compact=args.compact)
        if args.jobno == 1:
            for line in intlist.write_numbers_vertically(sitelist):
                print(line)
            seqs = it.chain([first],seqs)
        for s in seqs:
            print( "".join(s.seq[n] for n in ndxlist), s.name, end=" ")
            if args.showmuts:
                print(m_mgr.get_mutation(s.seq),end="")
            print()

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    def vcount(seqs,*p,**kw):
        '''count items in generator as they go by'''
        if _args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        return seqs

    main(_args)
