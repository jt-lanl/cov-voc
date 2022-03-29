'''Find sequences that are far (in Hamming distance) from reference sequence'''

from operator import itemgetter
from collections import defaultdict
import argparse

import covid
import sequtil
import mutant
from verbose import verbose as v
from hamming import hamming


def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("-T",type=int,default=10,
        help="show the top T sequences, sorted by Hamming distance")
    paa("--ref","-r",
        help="mutant string to be used as reference sequence")
    paa("--nearby",action="store_true",
        help="look for small instead of large Hamming distances")
    paa("--mincount","-m",type=int,default=1,
        help="only show seqs that appear at least m times")
    paa("--output","-o",
        help="put the large hamming distances sequences into a fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def earliest_sequence(seqs):
    '''from a list of sequences, return the one with the earliest date'''
    slist = [(covid.date_from_seqname(s.name),s) for s in seqs]
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
        refmut = mutant.Mutation(args.ref)
        refseq = m_mgr.seq_from_mutation(refmut)
    else:
        refseq = first.seq

    ## deal with duplicates
    seq_cases = defaultdict(list)
    for s in seqs:
        seq_cases[s.seq].append(s)

    hamlist = [ (hamming(refseq,seq), seq_cases[seq]) for seq in seq_cases ]
    hamlist = sorted(hamlist,
                     key=itemgetter(0),
                     reverse=False if args.nearby else True)

    outseqs=[]
    count=0
    for (hamdist,seqlist) in hamlist:
        if count >= args.T:
            break

        seq_len = len(seqlist)
        if seq_len < args.mincount:
            continue

        s = earliest_sequence(seqlist)
        print("%3d %4d %s" % (hamdist,seq_len,s.name))
        print("        ",m_mgr.get_mutation(s.seq))

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
