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
    paa("--range","-R",type=int,nargs=2,
        help="Restrict attention to this range of sites in the sequence")
    paa("--na",action="store_true",
        help="set for nucleotide alignment (default is amino acid alignment)")
    paa("--speedhack",type=int,default=0,
        help="for full-range spike, 70 is a good nubmer")
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
    ## if not, then they are inconsistent

    inconsistent=set()
    gseq_with_dseq = dict()
    for gseq in subseqset:
        dseq = de_gap(gseq)
        if dseq in gseq_with_dseq:
            gseq_already = gseq_with_dseq[dseq]
            if gseq != gseq_already:
                assert len(gseq) == len(gseq_already)
                gseq_a,gseq_b = sorted([gseq,gseq_already])
                inconsistent.add((gseq_a,gseq_b))
                v.vvprint(f'{gseq_a} != {gseq_b} : dseq={dseq}')
        else:
            gseq_with_dseq[dseq]=gseq
    return inconsistent


## sequence trimming routines basically cosmetic
## just used so human reader acn more readily see
## what the relevant differences are in two sequence segments
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
    '''Given two seqs (eg, 'FL-G-V--TSR-F' and 'FL-G-VT-S-R-F')
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
    '''
    given a gappy substring gaseq (such as 'R---T-T-')
    and the site value we are starting with,
    return the Mutation object associated with the gseq
    '''
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
            raise RuntimeError(f'ref value at site {site} is {ref} '
                               '-- this should never happen')
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
    def fprint(*args,**kwargs):
        print(*args,**kwargs,file=file)

    dseq = de_gap(gseqa)
    assert dseq == de_gap(gseqb)

    ga,gb = trimseqs(gseqa,gseqb)
    dg = de_gap(ga)
    fprint(f'{lo:>4d}:{hi:<4d}: {gseqa} : {ga}  :{dg}')
    fprint(f'           {gseqb} : {gb}  :{dseq}')
    fprint(f'{mstr_a} {mstr_b}')

## fcns used for profiling (and assessing good values of speedhack)
def setify_a(gen):
    return set(gen)

def setify_b(gen):
    return set(gen)

def setify_c(gen):
    return set(gen)

def _main(args):
    '''checkalign main'''
    v.vprint(args)

    speedhack=args.speedhack

    seqs = sequtil.read_seqfile(args.input)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    xlator = mutant.SiteIndexTranslator(first.seq)
    mut_mgr = mutant.MutationManager(first.seq)

    ## range (default is whole sequence)
    rlo,rhi = args.range if args.range else (1,xlator.topsite+1)
    rlo = max(rlo,1)
    rhi = min(rhi,xlator.topsite+1)
    ## convert from site number of index number
    ndxmin = min(xlator.indices_from_site(rlo))
    ## pad on the right (high) end by windowsize
    rhipad = rhi + args.windowsize
    rhipad = min(rhipad,xlator.topsite+1)
    ndxmax = max(xlator.indices_from_site(rhipad-1))+1

    v.vvprint('Reading sequence file...',end='')
    ## Read sequences, trim to range,
    ## and discard duplicates by keeping them in a set()
    seqset = setify_a(s.seq[ndxmin:ndxmax] for s in seqs)
    v.vvprint('ok')

    ## Having finished setup, go look for inconsistencies
    ## ie, window-length intervals with inconsistent substrings
    bad_intervals=[]
    for lo in range(rlo,rhi):
        hi = min(rhi,lo+args.windowsize)
        ndxlo = min(xlator.indices_from_site(lo))
        if speedhack and lo%speedhack == 0:
            ## trim from the left (not necessary, seems to speed up)
            seqset = setify_b(seq[ndxlo-ndxmin:] for seq in seqset)
            ndxmin = ndxlo
        ## now trim from the right
        ndxhi = max(xlator.indices_from_site(hi-1))+1
        subseqset = setify_b(seq[ndxlo-ndxmin:ndxhi-ndxmin] for seq in seqset)
        if args.xavoid:
            subseqset = set(seq for seq in subseqset if 'X' not in seq)
            if len(subseqset)==0:
                continue
        for xhi in range(hi,lo,-1):
            ## decrease subseq length, by truncating last site
            ## (may be several characters)
            xndxhi = max(xlator.indices_from_site(xhi-1))+1
            v.vvvprint(f'Checking sites {lo}:{xhi} / '
                       f'local indices {ndxlo-ndxmin}:{xndxhi-ndxmin} / '
                       f'global indices {ndxlo}:{xndxhi}')
            subseqset = setify_c(seq[:xndxhi-ndxlo] for seq in subseqset)
            inconsistent = check_subsequences(subseqset)
            bad_intervals.extend((lo,xhi,ndxlo,xndxhi,gsa,gsb)
                                 for gsa,gsb in inconsistent)
            ## show inconsistencies _as_ they are found
            if args.verbose > 1:
                for gsa,gsb in inconsistent:
                    msa,msb = get_inconsistent_mstringpair(mut_mgr,
                                                           lo,ndxlo,
                                                           gsa,gsb)
                    show_inconsistency(lo,hi,gsa,gsb,msa,msb)


    v.vprint('Bad intervals:',len(bad_intervals),f'in range {rlo}-{rhi}')
    mstringpairs = set()
    for lo,hi,ndxlo,ndxhi,gseqa,gseqb in bad_intervals:
        mstr_a,mstr_b = get_inconsistent_mstringpair(mut_mgr,
                                                     lo,ndxlo,
                                                     gseqa,gseqb)
        if args.verbose == 1:
            # if args.verbose > 1, we've already shown them!
            show_inconsistency(lo,hi,gseqa,gseqb,mstr_a,mstr_b)
        mstringpairs.add((mstr_a,mstr_b))

    v.print('Distinct inconsistencies:',len(mstringpairs),
            f'in range {rlo}:{rhi}')
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
