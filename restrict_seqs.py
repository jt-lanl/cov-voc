'''Restrict a sequence (eg fasta) file to those sequences for which ISL numbers are in a separate list (perhaps drawn from names of another file)'''

import argparse
import verbose as v
import sequtil
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--reference","-r",
        help="Name of reference list of ISL #s")
    paa("--output","-o",
        help="Write restricted seqeunces to this file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    v.vprint(f'Reading reference file: {args.reference}')
    refseqs = sequtil.read_seqfile(args.reference)
    _,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    ref_isls = set(covid.get_isl(s.name) for s in refseqs)
    v.vprint(f'Obtained {len(ref_isls)} ISL numbers from {args.reference}')

    outseqs = [first]
    orig_seqcount=1
    for s in seqs:
        orig_seqcount += 1
        if covid.get_isl(s.name) in ref_isls:
            outseqs.append(s)

    v.vprint(f'Read {orig_seqcount} sequences from {args.input}')
    v.vprint(f'Writing {len(outseqs)} sequences to {args.output}')

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
