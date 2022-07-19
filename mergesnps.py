'''Find all the mutations (relative to first seq) in a set of sequences
   and merge them into a single sequence that is appended to the sequence list
'''
## Because this uses apobec.get_all_mutsites(), it only works for DNA sequences

import argparse
from collections import defaultdict,Counter
import scipy.stats as stats

import verbose as v
import sequtil
import apobec as apo

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input sequence file")
    paa("--output","-o",
        help="output sequence file, "
        "with merged SNPs as second sequence in list")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = sequtil.read_seqfile(args.input)
    seqs = list(seqs)
    v.vprint(f'Read {len(seqs)} sequence files')
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    all_muts = defaultdict(list)
    for s in seqs:
        mut_sites = apo.get_all_mutsites(first.seq,s.seq)
        for n in mut_sites:
            all_muts[n].append(s.seq[n])

    v.vprint(f'Fount {len(all_muts)} mutation sites')
    muts = dict()
    for n,mutlist in all_muts.items():
        cnt = Counter(mutlist)
        [(cmut,_)] = cnt.most_common(1)
        muts[n] = cmut

    ## Now create the merged sequence that contains the most common
    ## mutation at each site

    merged_list = list(first.seq)
    for n,cmut in muts.items():
        merged_list[n] = cmut
    merged_str = "".join(merged_list)
    merged_seq = sequtil.SequenceSample('Merged-SNPs',merged_str)

    if args.output:
        v.vprint(f"Writing {2+len(seqs)} sequences")
        sequtil.write_seqfile(args.output,
                              [first] + [merged_seq] + seqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
