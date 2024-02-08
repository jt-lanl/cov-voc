'''Insert dashes at a given site for these fasta files'''

import argparse
import verbose as v
import intlist
import readseq

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input fasta file")
    paa("--output","-o",
        help="output fasta file")
    paa("-n",type=int,default=1,
        help="number of dashes to insert")
    paa("--position","-p",type=int,default=1,
        help="position before which to insert dashes")
    paa("--seqs","-s",
        help="list of which sequences to add dashes (1=first, default=all)")
    paa("--dash",default="-",
        help="dash character")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = readseq.read_seqfile(args.input)
    seqs = list(seqs)

    whichseqs = range(1,1+len(seqs))
    if args.seqs:
        whichseqs = intlist.string_to_intlist(args.seqs)
    whichseqs = set(whichseqs)
    
    p = args.position-1
    for n,s in enumerate(seqs,start=1):
        if n in whichseqs:
            s.seq = s.seq[:p] + args.dash * args.n + s.seq[p:]

    if args.output:
        readseq.write_seqfile(args.output,seqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
