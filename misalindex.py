'''check (mis)alignemnt (in)consistency based on index instead of site number'''

import sys
from collections import Counter
import argparse

import verbose as v
import intlist
from readseq import xopen
import mutant
import covid

import checkalign as ca

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
    paa("--msegpairs","-M",
        help="write file of msegment pairs")
    paa("--xavoid",action="store_true",
        help="Avoid X's in the mutation inconsistencies")
    paa("--windowsize","-w",type=int,default=0,
        help="subsequence window size")
    paa("--fracrange","-F",type=int,nargs=2,default=(1,1),
        help="Fractional range: eg, -R 3 4 restricts attention"
        "to third quarter of sequence")
    paa("--speedhack",type=int,default=0,
        help="for full-range spike, 50 is a good nubmer")
    paa("--speedhackhelp",action="store_true",
        help="Use this option to print out a discussion of speedhack")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    ## fix windowsize
    if not args.windowsize:
        args.windowsize = 90
    return args


def show_inconsistency(lo,gseqa,gseqb,na,nb,file=sys.stderr):
    '''write a summary of the inconsistency'''
    def fprint(*args,**kwargs):
        print(*args,**kwargs,file=file)

    hi = lo+len(gseqa)
    dseq = ca.de_gap(gseqa)
    assert dseq == ca.de_gap(gseqb)

    ga,gb = ca.trimseqs(gseqa,gseqb)
    dg = ca.de_gap(ga)
    fprint(f'ndx {lo:>6d}:{hi:<6d}: {gseqa} : {ga}  :{dg}')
    fprint(f'cnt {na:>6d},{nb:<6d}: {gseqb} : {gb}  :{dseq}')

def slice_updatecounts(subseqs,lo,hi):
    '''slice the subseqs Counter, and update
    with counts for the sliced subseqs
    '''
    counter=Counter()
    for seq,cnt in subseqs.items():
        counter[ seq[lo:hi] ] += cnt
    return counter


def _main(args):
    '''checkalign main'''
    v.vprint(args)

    if args.speedhackhelp:
        print(ca.SPEEDHACK_HELP)
        return

    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    xlator = mutant.SiteIndexTranslator(first.seq)
    ndxtop = len(first.seq)

    num,den = args.fracrange
    num = max(num,1)
    assert num <= den
    ndxmin = (num-1)*len(first.seq)//den
    ndxhi = num*len(first.seq)//den
    ndxmax = ndxhi + args.windowsize
    ndxmax = min(ndxmax,ndxtop)
    v.vvprint(f'{ndxmax=}, {ndxmin=} {ndxhi=},{ndxmax=}')
    v.vvprint('Reading input file...',end="")
    seqset = Counter(s.seq[ndxmin:ndxmax] for s in seqs)
    v.vvprint('ok')

    ## Having finished setup, go look for inconsistencies
    ## ie, window-length intervals with inconsistent substrings
    bad_intervals=[]
    for lo in range(ndxmin,ndxhi):
        hi = lo+args.windowsize
        hi = min(ndxtop,hi)
        if args.speedhack and (lo-ndxmin+1)%args.speedhack == 0:
            ## slice from the left (not necessary, seems to speed things up)
            seqset = slice_updatecounts(seqset,lo-ndxmin,0)
            ndxmin = lo
        ## now slice from the right
        subseqset = slice_updatecounts(seqset,lo-ndxmin,hi-ndxmin)
        if args.xavoid:
            xkeys = set(seq for seq in subseqset if 'X' in seq)
            for xkey in xkeys:
                subseqset.pop(xkey)
            if len(subseqset)==0:
                continue
        for xhi in range(hi,lo+1,-1):
            ## decrease subseq length, by truncating last site
            ## (may be several characters)
            subseqset = slice_updatecounts(subseqset,0,xhi-lo)
            v.vvvprint(f'local indices {lo-ndxmin}:{xhi-ndxmin} / '
                       f'global indices {lo}:{xhi} / '
                       f'sample: {len(subseqset)} {next(iter(subseqset))}')

            inconsistent = ca.check_subsequences(subseqset)
            bad_intervals.extend((lo,gsa,gsb,subseqset[gsa],subseqset[gsb])
                                 for gsa,gsb in inconsistent)
            ## show inconsistencies _as_ they are found
            if args.verbose > 1:
                for gsa,gsb in inconsistent:
                    show_inconsistency(lo,gsa,gsb,subseqset[gsa],subseqset[gsb])

    v.vprint('Bad intervals:',len(bad_intervals),f'in range {ndxmin}:{ndxhi}')
    uniq_segments = set()
    for lo,gseqa,gseqb,cnta,cntb in bad_intervals:
        if cnta > cntb:
            ## show the rarer form first
            cnta,cntb = cntb,cnta
            gseqa,gseqb = gseqb,gseqa
        if args.verbose == 1:
            # if args.verbose > 1, we've already shown them!
            show_inconsistency(lo,gseqa,gseqb,cnta,cntb)
        ga,gb = ca.trimseqs(gseqa,gseqb)
        if ga == gseqa:
            uniq_segments.add((lo,gseqa,gseqb,cnta,cntb))
            hi = lo+len(ga)
            sitelist = [xlator.site_from_index(ndx)
                        for ndx in range(lo,hi)]
            for line in intlist.write_numbers_vertically(sitelist):
                print(line)
            print(first.seq[lo:hi],"ref")
            print(gseqa,cnta)
            print(gseqb,cntb)
            print()

    if args.msegpairs and uniq_segments:
        ## write file only if uniq_segments is not empty
        with xopen(args.msegpairs,'w') as fout:
            if num==1:
                print("Index:Range  Site~Range (counts,counts) "
                      "segment / segment",file=fout)
            for lo,ga,gb,cnta,cntb in sorted(uniq_segments):
                hi = lo+len(ga)
                losite = xlator.site_from_index(lo)
                hisite = xlator.site_from_index(hi)
                print(f'{lo:>5d}:{hi:<5d} '
                      f'{losite:>5d}~{hisite:<5d} '
                      f'({cnta:>6d},{cntb:<6d}) {ga} / {gb}',file=fout)

if __name__ == '__main__':

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
