'''Tweak alignment according to specifed mutant strings'''

import itertools as it
from collections import Counter
import argparse

import verbose as v
import intlist
import covid
import sequtil
import mutant
import tweak as tku
from tweak import IndexTweak

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--jobno",type=int,default=1,
        help="job number if using tweakalign in parallel")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    covid.corona_args(argparser)
    paa = argparser.add_argument_group("Apply Tweak Options").add_argument
    paa("--mstrings","-m",nargs=2,action="append",
        help="pair of from/to mstrings")
    paa("--mfile","-M",action="append",
        help="read from/to mstring pairs from a file")
    paa("--tfile","-T",action="append",
        help="read indexed tweaks from a file")
    paa("--checktweaks",action="store_true",
        help="Check for dependencies among the tweaks")
    paa("--reverse",action="store_true",
        help="reverse the effect of tweak: b->a instead of a->b")
    paa("--output","-o",
        help="output tweaked fasta file")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def showtweakpair(xlator,ts,to):
    '''given a pair of tweaks (ts,to), show the result
    of applying the tweaks; show two results, depending
    on which tweaks is applied first'''
    ndxmin = min(ts.ndxlo,to.ndxlo)
    ndxmax = max(ts.ndxhi,to.ndxhi)

    oslice = slice(to.ndxlo-ndxmin,to.ndxhi-ndxmin)
    sslice = slice(ts.ndxlo-ndxmin,ts.ndxhi-ndxmin)
    ttslices = [(ts,sslice),(to,oslice)]

    untweaked = ["."] * (ndxmax-ndxmin)
    for tt,tslice in ttslices:
        untweaked[tslice] = tt.sa
    tweaked_so = untweaked.copy()
    for tt,tslice in ttslices:
        if tweaked_so[tslice] == list(tt.sa):
            tweaked_so[tslice] = tt.sb
    tweaked_os = untweaked.copy()
    for tt,tslice in reversed(ttslices):
        if tweaked_os[tslice] == list(tt.sa):
            tweaked_os[tslice] = tt.sb

    sites = [xlator.site_from_index(ndx)
             for ndx in range(ndxmin,ndxmax)]

    ## following is roughly cut-and-pasted from IndexTweak.viz method
    lines = []
    lines.extend([f'n   {line}'
                  for line in intlist.write_numbers_vertically(range(ndxmin,ndxmax))])
    lines.extend([f's   {line}'
                  for line in intlist.write_numbers_vertically(sites)])
    lines.append( f' R: {xlator.refseq[ndxmin:ndxmax]} Reference')
    lines.append( f' U: {"".join(untweaked)} Untweaked' )
    lines.append( f' T: {"".join(tweaked_so)} if first {ts.ma} {ts.mb} then {to.ma} {to.mb}')
    lines.append( f' T: {"".join(tweaked_os)} if first {to.ma} {to.mb} then {ts.ma} {ts.mb}')
    return "\n".join(lines)

def get_interacting_tweaks(tweaklist):
    '''Look for pairs of tweaks in the tweaklist
    such that the effect of applying the two tweaks
    will depend on the order in which they are applied
    '''
    for ts,to in it.combinations(tweaklist,2):
        act = ts.interacts_with(to)
        if act:
            yield ts,to,act

def checktweaks(tweaklist):
    '''Look for pairs of tweaks in the tweaklist
    such that the effect of applying the two tweaks
    will depend on the order in which they are applied
    '''
    return [f'{ts} + {to} : {act[0]}->[{act[1]} != {act[2]}]'
            for ts,to,act in get_interacting_tweaks(tweaklist)]

def show_checkedtweaks(xlator,tweaklist):
    '''provides a more elaborate/informative description of what the
    differences are when tweak pairs are applied in different order
    '''
    showlines = []
    for ts,to,_ in get_interacting_tweaks(tweaklist):
        showlines.append(showtweakpair(xlator,ts,to))
    return "\n\n".join(showlines)

def apply_mstringpairs(seqs,mstringpairs,change_counter=None):
    '''apply tweaks, as defined by mseqpairs, to seqs'''
    for s in seqs:
        if covid.test_isref(s):
            ## if it's a ref seq; then
            ## we re-initialize  everything
            mut_mgr = mutant.MutationManager(s.seq)
            xpand = tku.expansion_needed(mut_mgr,mstringpairs)
            if xpand:
                s.seq = tku.expand_seq(xpand,s.seq)
                mut_mgr = mutant.MutationManager(s.seq)
            tweaklist = [tku.tweak_from_mstringpair(mut_mgr,ma,mb)
                         for ma,mb in mstringpairs]
            yield s
            continue

        if xpand:
            s.seq = tku.expand_seq(xpand,s.seq)
        for tweak in tweaklist:
            s.seq,wastweaked = tweak.apply_to_seq(s.seq)
            if wastweaked and change_counter is not None:
                change_counter[tweak] += 1
        yield s

def _main(args):
    '''tweakalign main'''
    v.vprint(args)

    changes_by_tweak = Counter()

    ## Begin reading sequences
    seqs = covid.read_filter_seqfile(args)

    ## Get all the mstring pairs
    mstringpairs = tku.get_mstring_pairs(args.mfile,args.mstrings)
    if args.reverse:
        mstringpairs = [(mb,ma) for ma,mb in mstringpairs]
    ## Get all the index-based tweaks
    if args.tfile:
        tweaklist = IndexTweak.tweaks_from_filelist(args.tfile)
        if args.reverse:
            for tweak in tweaklist:
                tweak.swap_ab()

    if args.checktweaks:
        if mstringpairs:
            first,seqs = covid.get_first_item(seqs,keepfirst=True)
            mut_mgr = mutant.MutationManager(first.seq)
            tweaklist = [tku.tweak_from_mstringpair(mut_mgr,ma,mb)
                         for ma,mb in mstringpairs]
        warnlist = checktweaks(tweaklist)
        v.print("\n".join(warnlist))

        first,seqs = covid.get_first_item(seqs,keepfirst=True)
        xlator = mutant.SiteIndexTranslator(first.seq)
        v.print( show_checkedtweaks(xlator,tweaklist) )

        return

    if mstringpairs:
        seqs = apply_mstringpairs(seqs,mstringpairs,changes_by_tweak)
    elif args.tfile:
        tweaklist = IndexTweak.tweaks_from_filelist(args.tfile)
        for s in seqs:
            for tweak in tweaklist:
                s.seq,wastweaked = tweak.apply_to_seq(s.seq)
                if wastweaked:
                    changes_by_tweak[tweak] += 1

    if args.output:
        sequtil.write_seqfile(args.output,seqs)

    for _ in seqs:
        ## in case seq not written, empty the generator
        pass

    if changes_by_tweak:
        v.vprint("Changes by tweak:")
        for tweak,cnt in changes_by_tweak.items():
            v.vprint(f'{tweak} (count={cnt})')

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
