'''Tweak alignment according to specifed mutant strings'''
import itertools as it
from collections import Counter
import random
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
    paa("--shuffletweaks",action="store_true",
        help="apply the mstring-pair tweaks in random order")
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

def _main(args):
    '''tweakalign main'''
    v.vprint(args)

    ## Initial set of tweaks from tweakfile(s)
    tweaklist = IndexTweak.tweaks_from_filelist(args.tfile)

    ## Begin reading sequences
    seqs = covid.read_filter_seqfile(args)

    ## Get all the mstring pairs
    mstringpairs = tku.get_mstring_pairs(args.mfile,args.mstrings)
    seqs,xtra_dashes = tku.add_needed_dashes(seqs,mstringpairs)
    if xtra_dashes:
        v.print("Sequences modified, added dashes:",xtra_dashes)
    if xtra_dashes and args.tfile:
        raise RuntimeError(f'Since seqs expanded, cannot trust {args.tfile}')

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    mut_mgr = mutant.MutationManager(first.seq)

    tweaklist = [tku.tweak_from_mstringpair(mut_mgr,ma,mb)
                 for ma,mb in mstringpairs]

    if args.shuffletweaks:
        tweaklist = random.sample(tweaklist,k=len(tweaklist))

    ## remove duplicates (preserve order in v3.7+)
    tweaklist = list(dict.fromkeys(tweaklist))

    if args.checktweaks:
        warninglist = checktweaks(tweaklist)
        if warninglist:
            v.print("\n".join(warninglist))
        v.print(f'Found {len(warninglist)} interacting pairs '
                f'among {len(tweaklist)} tweaks')
        return

    v.vprint("Reading all sequences...")
    seqs = list(seqs)
    v.vprint("Finished reading sequences")

    changes_by_tweak = Counter()
    for s in seqs:
        for tweak in tweaklist:
            s.seq,wastweaked = tweak.apply_to_seq(s.seq)
            if wastweaked:
                changes_by_tweak[tweak] += 1

    if changes_by_tweak:
        v.vprint("Changes by tweak:")
        for tweak,cnt in changes_by_tweak.items():
            v.vprint(f'{tweak} (count={cnt})')

    if args.output:
        seqs = it.chain([first],seqs)
        sequtil.write_seqfile(args.output,seqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
