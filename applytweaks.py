'''Tweak alignment according to specifed mutant strings'''
import itertools as it
from collections import Counter
import argparse

import verbose as v
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
        help="Check for dependencies among the tweaks [EXPERIMENTAL!]")
    paa("--reverse",action="store_true",
        help="reverse the effect of tweak: b->a instead of a->b")
    paa("--output","-o",
        help="output tweaked fasta file")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args


def checktweaks(tweaklist):
    '''Look for pairs of tweaks in the tweaklist
    such that the effect of applying the two tweaks
    will depend on the order in which they are applied
    '''
    warninglist=[]
    for ts,to in it.combinations(tweaklist,2):
        act = ts.interacts_with(to)
        if act:
            msg = f'{ts} + {to} : {act[0]}->[{act[1]} != {act[2]}]'
            warninglist.append(msg)
    return warninglist

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

    if changes_by_tweak:
        v.vprint("Changes by tweak:")
        for tweak,cnt in changes_by_tweak.items():
            v.vprint(f'{tweak} (count={cnt})')

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
