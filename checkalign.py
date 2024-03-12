'''check alignemnt consistency and provide suggested alignment tweaks'''
## But: unlike fixalign, do not bother trying to actually align!

import sys
import re
from collections import Counter,defaultdict
from functools import lru_cache
import argparse

import verbose as v
from breakpipe import no_broken_pipe
from xopen import xopen
import mutant
import covid
from tweak import IndexTweak,de_gap

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
    paa("--bysite",action="store_true",
        help="Only look for alignment anomalies that are on site indices")
    paa("--mstringpairs","-M",
        help="write file of mstring pairs [NOT IMPLEMENTED]")
    paa("--xavoid",action="store_true",
        help="Avoid X's in the mutation inconsistencies")
    paa("--windowsize","-w",type=int,default=0,
        help="subsequence window size")
    paa("--fracrange","-F",type=int,nargs=2,default=(1,1),
        help="Restrict attention to this range of sites in the sequence")
    paa("--na",action="store_true",
        help="set for nucleotide alignment (default is amino acid alignment)")
    paa("--speedhack",type=int,default=0,
        help="for full-range spike, 50 is a good nubmer")
    paa("--speedhackhelp",action="store_true",
        help="Use this option to print out a discussion of speedhack")
    paa("--viz",
        help="Write vizzy tables into this file")
    paa("--output","-o",
        help="Write output lines that can be used as input for applytweaks")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    ## fix fracrange
    num,den = args.fracrange
    num = max(num,1)
    den = max(den,1)
    assert num <= den
    args.fracrange = [num,den]
    ## fix windowsize
    if not args.windowsize:
        args.windowsize = 90 if args.na else 30
    if args.na:
        args.windowsize = 3*(args.windowsize//3)
    return args

SPEEDHACK_HELP='''
Runtime for this routine is dominated by three processes:
1/ reading the sequence file
2/.string slicing (eg, seq[ndxlo:ndxhi]) to create subsequencs
3/ converting generators to sets

Note that the identification of which subsequences match, while
central to checking the alignment, does not contribute substantially
to the runtime of the code.

Also note, the reason we use sets (instead of lists, say) is that
by eliminating duplicates, we have fewer subsequences to work
with. Uses less memory, and there's less data to process, so uses
less runtime.  (In fact, since we read the data into a generator
instead of a list to begin with, we /never/ require the memory
needed to hold the entire fasta file's worth of data.)

The speedhack deals with the interaction between #2 and #3. The
more aggressively we use sets to reduce the data size, the fewer
subsequecnes we need to be slicing through.  Using speedhack=0
(default), we do not separately left-truncate and right-truncate,
so there are fewer calls to slice and fewer set-ify.  But the number
of subsequences to which those slices and set'ifings are applied is
larger.  With speedhack=1, we left-truncate and right-truncate
separately.  More calls to slice and to set-ify, but the slicing and
set-ifying is done to fewer subsequences. Sometimes that's a win,
sometimes not so much.

But if we take speedhack=10, say, then we will left-truncate only on
every 10th loop.  So the amount of extra slicing and set-ifying
is marginal, but (especially for a long sequence such as spike protein),
the effect is to continue reducing the size of the sequence set, so
reducing how much slicing and set-ifying is done at each call.

For spike protein, I find that speedhack=50-ish works well, and
delivers a 2-3x speedup, depending on the data.
'''

def make_pairs(gseqlist,ctr):
    '''given N gseqs in gseqlist, return a list of N-1 pairs,
    where the second element of every pair is tht gseq with
    the highest count'''
    gseqs = sorted(gseqlist,key=ctr.get,reverse=True)
    return [(gseq,gseqs[0]) for gseq in gseqs[1:]]

def check_subsequences(subseqctr):
    '''Check for inconsistent subsequences
    Two gappy subsequences are inconsistent if:
        1/ they are not identical themselves
        2/ they share the same de-gapped sequence
    Args:
       subseqctr: Counter() object number of seqs in which subseq has appeared
    Returns:
        tweaklist: list of "tweaks"; ie, incnsistent pairs of subsequences

    '''
    tweaklist=[]
    gseq_with_dseq = defaultdict(list)
    for gseq in subseqctr:
        dseq = de_gap(gseq)
        gseq_with_dseq[dseq].append(gseq)
    for dseq,gseqlist in gseq_with_dseq.items():
        if len(gseqlist) > 1:
            tweaks = make_pairs(gseqlist,subseqctr)
            tweaklist.extend(tweaks)
    return tweaklist

def _trimfront(ga,gb):
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
        ga,gb = _trimfront(ga,gb)
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
    minsite = min(ssm.site for ssm in (bmut|gmut) if ssm.ref != "+")
    maxsite = max(ssm.site for ssm in (bmut|gmut))
    for ssm in commonmut:
        if minsite <= ssm.site <= maxsite:
            bmut.add(ssm)
            gmut.add(ssm)

    ## Adjustment 3: add back explicit identity ssm's (eg, 'S247S')
    for site in range(minsite,maxsite+1):
        ref = mut_mgr.refval(site)
        if ref=='-':
            raise RuntimeError(f'ref value at site {site} is "{ref}" '
                               '-- this should never happen')
        if site not in [ssm.site for ssm in bmut if ssm.ref != "+"]:
            bmut.add(mutant.SingleSiteMutation(ref,site,ref))
        if site not in [ssm.site for ssm in gmut if ssm.ref != "+"]:
            gmut.add(mutant.SingleSiteMutation(ref,site,ref))

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

def update_seq_counter(inctr,lo,hi=None):
    '''truncate seqs to hi:lo, and
    return a new counter that keeps count of the truncated seqs'''
    ctr = Counter()
    lohi = slice(lo,hi)
    for seq,cnt in inctr.items():
        ctr[seq[lohi]] += cnt
    return ctr

@no_broken_pipe
def _main(args):
    '''checkalign main'''
    v.vprint(args)

    if args.speedhackhelp:
        print(SPEEDHACK_HELP)
        return

    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,keepfirst=True)
    mut_mgr = mutant.MutationManager(first.seq)

    num,den = args.fracrange
    ndxminlo = ndxmin = (num-1)*len(first.seq)//den
    ndxmaxlo = num*len(first.seq)//den
    ndxmaxhi = ndxmaxlo + args.windowsize
    ndxmaxhi = min(ndxmaxhi,len(first.seq))

    sites = sorted(set(mut_mgr.site_from_index(ndx)
                       for ndx in range(ndxminlo,ndxmaxhi)))
    site_indexes = [mut_mgr.index_from_site(site)
                    for site in sites]
    ## also, use site indexes that are +1 from site, enables [+251V,...] style tweaks
    site_indexes.extend(mut_mgr.index_from_site(site)+1
                        for site in sites)
    site_indexes = set(site_indexes)

    ## Fill seqctr with seqs
    v.vprint("Reading input, filling Counter...",end="")
    seqctr = Counter(s.seq[ndxminlo:ndxmaxhi] for s in seqs)
    v.vprint("ok. Distinct (sub)sequences",len(seqctr))

    inconsistencies=[]
    for ndxlo in range(ndxminlo,ndxmaxlo):
        if args.bysite and ndxlo not in site_indexes:
            continue
        ndxhi = min(ndxlo + args.windowsize,ndxmaxhi)
        if args.speedhack and (ndxlo-ndxmin+1)%args.speedhack == 0:
            ## trim the main seq counter from the left
            seqctr = update_seq_counter(seqctr,ndxlo-ndxmin,hi=None)
            ndxmin = ndxlo
        ## make secondary subseq counter from main seq counter
        v.vvvprint('=slice: ',ndxlo-ndxmin,ndxhi-ndxmin)
        subseqctr = update_seq_counter(seqctr,ndxlo-ndxmin,ndxhi-ndxmin)
        for aux_ndxhi in range(ndxhi,ndxlo+1,-1):
            ## Loop over ranges lo:lo+2, lo:lo+3, ... ln:hi
            ## But do it backward, truncating subseq's at each step
            ## next two lines commented out; enables use to find
            ## more bysite tweaks
            #if args.bysite and aux_ndxhi-1 not in site_indexes:
            #    continue
            if args.bysite:
                v.vvvprint('sites:'
                           f'{mut_mgr.site_from_index(ndxlo)}-'
                           f'{mut_mgr.site_from_index(aux_ndxhi-1)}',
                           end=" ")
                v.vvvprint(f' ndx: {ndxlo}-{aux_ndxhi-1}',end="")
            v.vvvprint(f' slice: 0:{aux_ndxhi-ndxlo}')
            subseqctr = update_seq_counter(subseqctr,0,aux_ndxhi-ndxlo)
            tweaklist = check_subsequences(subseqctr)
            for (gsa,gsb) in tweaklist:
                if args.xavoid and 'X' in gsa:
                    continue
                inconsistencies.append(IndexTweak(ndxlo,aux_ndxhi,
                                                  gsa,gsb,
                                                  ca=subseqctr[gsa],
                                                  cb=subseqctr[gsb]))

    v.vprint(f'Inconsistencies: {len(inconsistencies)}',
             f'in range {ndxminlo}:{ndxmaxlo}')

    for inc in inconsistencies:
        v.vprint(inc)

    if args.bysite:
        minimal_incs = IndexTweak.get_minimal(inconsistencies)
    else:
        minimal_incs = [inc for inc in inconsistencies
                        if inc.is_trim()]

    if args.bysite:
        for inc in minimal_incs:
            lo = mut_mgr.site_from_index(inc.ndxlo)
            if mut_mgr.index_from_site(lo) < inc.ndxlo:
                ## grit-teeth, kind of a hack...
                ## needed for [+251V,...] style mutations
                lo += 1
            muta = substr_to_mut(first.seq,lo,inc.ndxlo,inc.sa)
            mutb = substr_to_mut(first.seq,lo,inc.ndxlo,inc.sb)
            mstra,mstrb = mutpair_to_mstringpair(mut_mgr,muta,mutb)
            inc.ma = mstra
            inc.mb = mstrb
            v.vprint(f'{inc} {mstra} {mstrb}')

    if args.viz and minimal_incs:
        with xopen(args.viz,"w") as vizout:
            print("\n\n".join(inc.viz(mut_mgr,showcontext=True)
                              for inc in minimal_incs),
                  file=vizout)

    if args.output:
        with xopen(args.output,"w") as fout:
            for inc in minimal_incs:
                print(inc,file=fout)

if __name__ == '__main__':

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
