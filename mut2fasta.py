DESCRIPTION='''
mut2fasta: takes one or more mutant strings (either from a file or from the command line), 
each of which looks something like "[A222V,A262S,S494P,D614G]", and produces a fasta file, 
each sequence of which corresponds to a mutant specified by the string, relative to the 
reference sequence, which is the first sequence in the specified reference fasta file.
'''
import sys
import re
import warnings
import argparse

import readseq
import mutant

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    paa("--input","-i",
        help="input fasta file with reference sequence first")
    paa("--skipref",action="store_true",
        help="do not include reference sequence in fasta output file")
    paa("--mutstrings","-m",nargs='+',
        help="one or more strings of the form '[L5F,...,Q957R]'")
    paa("--mutfile",
        help="input file with list of mutant strings")
    paa("--nodash",action="store_true",
        help="remove dashes from the output sequences")
    paa("--output","-o",
        help="ouptut fasta file of mutated seqeunces")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def rd_mutfile(filename):
    with open(filename) as f:
        for line in f:
            line = re.sub("#.*","",line)
            line = line.strip()
            if not line:
                continue
            yield mutant.Mutation.from_mstring(line)

def main(args):

    def de_dash(s):
        if args.nodash:
            s.seq = re.sub("-","",s.seq)
        return s

    mutlist = []
    if args.mutfile:
        mutlist.extend( rd_mutfile(args.mutfile) )
    if args.mutstrings:
        mutlist.extend( mutant.Mutation(mstring)
                        for mstring in args.mutstrings)

    vprint(len(mutlist),"mutant strings have been read")
    for m in mutlist:
        vvprint("  ",m)

    if len(mutlist) == 0:
        raise RuntimeError("Must specify mutlist with either --mutfile for --mutstirngs or -m")

    ## Grab the first sequence from a reference file
    seqs = readseq.read_seqfile(args.input,maxseqs=1)
    seqs = list(seqs)
    wuhan = seqs[0].seq

    #if "-" in wuhan:
    #    warnings.warn("Dashes have been stripped from reference sequence")
    #    wuhan = re.sub("-","",wuhan)
    #    seqs[0].seq = wuhan

    MM = mutant.MutationManager(wuhan)
    
    mutseqs = []
    if not args.skipref:
        mutseqs.append(seqs[0])
    for mut in mutlist:
        mutseq = MM.seq_from_mutation(mut)
        mutseqs.append( readseq.SequenceSample( str(mut), mutseq ) )

    mutseqs = (de_dash(s) for s in mutseqs)        

    if args.output:
        readseq.write_seqfile(args.output,mutseqs)
    
    

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

