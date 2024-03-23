'''Find sequences that are far (in Hamming distance) from reference sequence'''

from operator import itemgetter
from collections import defaultdict
import argparse

import covid
import sequtil
import mutant
import verbose as v
from hamming import hamming


def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    covid.corona_args(argparser)
    paa = argparser.add_argument_group('Hamming Options').add_argument
    paa("-T",type=int,default=10,
        help="show the top T sequences, sorted by Hamming distance")
    paa("--ref","-r",
        help="mutant string (or baseline) to be used as reference sequence")
    paa("--nearby",action="store_true",
        help="look for small instead of large Hamming distances")
    paa("--mincount","-m",type=int,default=1,
        help="only show seqs that appear at least m times")
    paa("--output","-o",
        help="put the large hamming distances sequences into a fasta file")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def earliest_sequence(seqs):
    '''from a list of sequences, return the one with the earliest date'''
    slist = [(covid.get_date(s.name),s) for s in seqs]
    slist = [(d,s) for d,s in slist if d]
    if not slist:
        ## if none of the sequences has a good date
        ## then just return the first sequence in the original list
        return seqs[0]
    slist = sorted(slist,key=itemgetter(0))
    return slist[0][1]

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)

    if args.ref:
        refmstring = covid.BASELINE_MSTRINGS.get(args.ref,args.ref)
        refmut = mutant.Mutation.from_mstring(refmstring)
        refseq = m_mgr.regex_from_mutation(refmut,exact=True)
    else:
        refseq = first.seq

    ## deal with duplicates
    seq_cases = defaultdict(list)
    for s in seqs:
        seq_cases[s.seq].append(s)

    hamlist = [ (hamming(refseq,seq), seq_cases[seq]) for seq in seq_cases ]
    hamlist = sorted(hamlist,
                     key=itemgetter(0),
                     reverse=not bool(args.nearby))

    outseqs=[first]
    count=0
    for (hamdist,seqlist) in hamlist:
        if count >= args.T:
            break

        seq_len = len(seqlist)
        if seq_len < args.mincount:
            continue

        s = earliest_sequence(seqlist)
        print(f'{hamdist:3d} {seq_len:4d} {s.name}')
        print("        ",m_mgr.seq_to_mutation(s.seq))

        if args.output:
            for s in seqlist:
                s.name = f'H={hamdist}:{s.name}'
                outseqs.append(s)

        count += 1

        if args.output:
            sequtil.write_seqfile(args.output,outseqs)


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
