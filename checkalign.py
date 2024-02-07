'''check alignemnt consistency and provide suggested alignment tweaks'''
## But: unlike fixalign, do not bother trying to actually align!

import sys
import re
from functools import lru_cache
import argparse

import verbose as v
from readseq import xopen
import sequtil
import mutant
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    ## Generic options first
    gaa = argparser.add_argument
    gaa("--verbose","-v",action="count",default=0,
        help="verbose")
    ## Corona (Input) options
    covid.corona_args(argparser)
    ## Alignment options
    paa = argparser.add_argument_group('Alignment Options').add_argument
    paa("--mstringpairs","-M",
        help="write file of mstring pairs")
    paa("--xavoid",action="store_true",
        help="Avoid X's in the mutation inconsistencies")
    paa("--windowsize","-w",type=int,default=0,
        help="subsequence window size")
    paa("--range","-R",type=int,nargs=2,#default=(float('-inf'),float('inf')),
        help="Restrict attention to this range of sites in the sequence")
    paa("--na",action="store_true",
        help="set for nucleotide alignment (default is amino acid alignment)")

    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    ## fix windowsize
    if not args.windowsize:
        args.windowsize = 60 if args.na else 20
    if args.na:
        args.windowsize = 3*(args.windowsize//3)
    return args

DEDASH = re.compile("-")
@lru_cache(maxsize=None)
def de_gap(seq):
    '''remove '-'s from sequence'''
    return DEDASH.sub("",seq)

def check_subsequences(subseqset):
    '''returns a set of inconsistencies; or empty set if alignment is okay'''

    ## strategy is to ensure that if two "gapped" subseqs (gseqs)
    ## have the same de-gapped string (dseq), then they are identical

    v.vvvprint('distinct subsequences:',len(subseqset))

    lenrange = [len(seq) for seq in subseqset]
    if min(lenrange) != max(lenrange):
        v.print('Un equal lengths in input sequences!!:',
                min(lenrange),max(lenrange))

    inconsistent=set()
    gseq_with_dseq = dict()
    for gseq in subseqset:
        dseq = de_gap(gseq)
        v.vvvprint_only(5,'gseq/dseq:',f'gseq={gseq} dseq={dseq}')
        if dseq in gseq_with_dseq:
            gseq_already = gseq_with_dseq[dseq]
            if gseq != gseq_already:
                gseq_a,gseq_b = sorted([gseq,gseq_already])
                inconsistent.add((gseq_a,gseq_b))
                v.vprint(f'{gseq_a} != {gseq_b} : dseq={dseq}')
                if len(gseq_a) != len(gseq_b):
                    v.print("Unequal sequence lengths!!!")
        else:
            gseq_with_dseq[dseq]=gseq
    return inconsistent


def trimfront(ga,gb):
    '''given two strings, trim from the front if there are identical
    characters; eg: ABCLMNOP,ABCDEFOP -> LMNOP,DEFOP
    '''
    for n,(a,b) in enumerate(zip(ga,gb)):
        if a!=b:
            return (ga[n:],gb[n:])
    ## if strings equal, return empty string
    return ("","")

def trimseqs(ga,gb):
    '''Given two seqs (eg, 'FL-G-V--TS' and 'FL-G-VT-S-')
    return trimmed sequences that illustrate the "core differeces"
    (eg, '--TS' and 'T-S-') by trimming characters off the
    front and back
    '''
    ## trim off the front then the back (trim-reverse, do it twice)
    assert len(ga) == len(gb)
    for _ in range(2):
        ga,gb = trimfront(ga,gb)
        ga = ga[::-1]
        gb = gb[::-1]
    return ga,gb

def substr_to_mut(firstseq,lo,ndxlo,gseq):
    '''given a substring (such as 'R---T-T-') and the site value we are starting with,
    return an mstring associated with the gseq'''
    ndxhi = ndxlo + len(gseq)
    v.vvprint('first:',firstseq[ndxlo:ndxhi])
    v.vvprint(' gseq:',gseq)
    mut_mgr = mutant.MutationManager(firstseq[ndxlo:ndxhi])
    mut = mut_mgr.get_mutation(gseq)
    for ssm in mut:
        ssm.site += lo-1
    return mut


def mutpair_to_mstringpair(mut_mgr,badmut,goodmut):
    '''
    given a pair of Mutation objects, return a pair of m-strings...
    making adjustments to snsure the mstrings can be used to
    tweak future sequences
    '''

    def mstringify(mutlike):
        return str(mutant.Mutation(mutlike))

    ## Adjustment 1: remove common ssm's
    bmut = set(badmut)
    gmut = set(goodmut)
    commonmut = bmut & gmut
    bmut = bmut - commonmut
    gmut = gmut - commonmut

    if len(bmut|gmut) == 0:
        v.print('Problematic mutations: ',
                mstringify(badmut),mstringify(goodmut))
        return None,None

    ## Adjustment 2: add back common ssm's that are in range
    ## (note: at this piont bmut,gmut are _sets_ of ssm's, not Mutation objects)
    minsite = min(mut.site for mut in (bmut|gmut) if mut.ref != "+")
    maxsite = max(mut.site for mut in (bmut|gmut))
    for mut in commonmut:
        if minsite < mut.site < maxsite:
            bmut.add(mut)
            gmut.add(mut)

    ## Adjustment 3: add back explicit identity ssm's (eg, 'S247S')
    for site in range(minsite,maxsite+1):
        ref = mut_mgr.refval(site)
        if ref=='-':
            v.print(f'Warning: ref value at site {site} is {ref} -- this should never happen')
            continue
        if site not in [ssm.site for ssm in bmut]:
            bmut.add(mutant.SingleSiteMutation(f'{ref}{site}{ref}'))
        if site not in [ssm.site for ssm in gmut]:
            gmut.add(mutant.SingleSiteMutation(f'{ref}{site}{ref}'))

    return mstringify(bmut),mstringify(gmut)

def get_inconsistent_mstringpair(mut_mgr,lo,ndxlo,gseqa,gseqb):
    '''convert sequence fragments gseqa,gseqb into a pair of mstrings'''
    ma = substr_to_mut(mut_mgr.refseq,lo,ndxlo,gseqa)
    mb = substr_to_mut(mut_mgr.refseq,lo,ndxlo,gseqb)
    mstr_a,mstr_b = mutpair_to_mstringpair(mut_mgr,ma,mb)
    return mstr_a,mstr_b

def show_inconsistency(lo,hi,
                       gseqa,gseqb,
                       mstr_a,mstr_b,
                       file=sys.stderr):
    '''write a summary of the inconsistency'''
    def vprint(*args,**kwargs):
        print(*args,**kwargs,file=file)

    ga,gb = trimseqs(gseqa,gseqb)
    dseq = de_gap(gseqa)
    dg = de_gap(ga)
    vprint(f'{lo:>4d}:{hi:<4d}: {gseqa} : {ga}  :{dg}')
    vprint(f'           {gseqb} : {gb}  :{dseq}')
    vprint(f'{mstr_a} {mstr_b}')

def _main(args):
    '''checkalign main'''
    v.vprint(args)

    seqs = sequtil.read_seqfile(args.input)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    xlator = mutant.SiteIndexTranslator(first.seq)
    mut_mgr = mutant.MutationManager(first.seq)

    ## range (default is whole sequence)
    rlo,rhi = args.range if args.range else (1,xlator.topsite+1)
    ## pad range by window size
    rlo = max(1,rlo-args.windowsize)
    rhi = min(xlator.topsite+1,rhi+args.windowsize)
    ## convert from site number of index number
    ndxmin = min(xlator.indices_from_site(rlo))
    ndxmax = max(xlator.indices_from_site(rhi-1))+1

    v.vvprint('Reading sequence file...',end='')
    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    first.seq = first.seq[ndxmin:ndxmax]
    seqset = set(s.seq[ndxmin:ndxmax] for s in seqs)
    ## these new subsequences have lost their name
    seqs = [sequtil.SequenceSample('',seq) for seq in seqset]
    seqs = [first]+seqs ## the first can keep its name!
    v.vvprint('ok')

    ## Having finished setup, go look for inconsistencies
    bad_intervals=[]
    for lo in range(rlo,rhi):
        hi = min(rhi,lo+args.windowsize)
        ndxlo = min(xlator.indices_from_site(lo))      -ndxmin
        ndxhi = max(xlator.indices_from_site(hi-1))+1  -ndxmin
        subseqset = set(s.seq[ndxlo:ndxhi] for s in seqs)
        if args.xavoid:
            subseqset = set(seq for seq in subseqset if 'X' not in seq)
            if len(subseqset)==0:
                continue
        for xhi in range(hi,lo,-1):
            ## decrease subseq length, by truncating last site (may be several characters)
            xndxhi = max(xlator.indices_from_site(xhi-1))+1  -ndxmin
            if ndxmin == 0:
                v.vvvprint(f'Checking sites {lo}:{xhi} / indices {ndxlo}:{xndxhi}')
            else:
                v.vvvprint(f'Checking sites {lo}:{xhi} / '
                           f'local indices {ndxlo}:{xndxhi} / '
                           f'global indices {ndxlo+ndxmin}:{xndxhi+ndxmin}')
            if xhi<hi:
                subseqset = set(seq[:xndxhi-ndxlo] for seq in subseqset)
            inconsistent = check_subsequences(subseqset)
            ## put a verbose test here
            if args.verbose > 1:
                for gsa,gsb in inconsistent:
                    msa,msb = get_inconsistent_mstringpair(mut_mgr,lo,ndxlo,gsa,gsb)
                    show_inconsistency(lo,hi,gsa,gsb,msa,msb)

            bad_intervals.extend((lo,xhi,ndxlo+ndxmin,xndxhi+ndxmin,gsa,gsb)
                                 for gsa,gsb in inconsistent)

    v.vprint('Bad intervals:',len(bad_intervals),f'in range {rlo}-{rhi}')
    mstringpairs = set()
    for lo,hi,ndxlo,ndxhi,gseqa,gseqb in bad_intervals:
        mstr_a,mstr_b = get_inconsistent_mstringpair(mut_mgr,lo,ndxlo,gseqa,gseqb)
        if args.verbose:
            show_inconsistency(lo,hi,gseqa,gseqb,mstr_a,mstr_b)
        mstringpairs.add((mstr_a,mstr_b))

    v.print('Distinct inconsistencies:',len(mstringpairs),f'in range {rlo}-{rhi}')
    for ma,mb in sorted(mstringpairs):
        v.print(f'{ma} {mb}')

    if args.mstringpairs and mstringpairs:
        ## write file only if mstringpairs is not empty
        with xopen(args.mstringpairs,'w') as fout:
            for ma,mb in sorted(mstringpairs):
                print(f'{ma} {mb}',file=fout)

if __name__ == '__main__':

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
