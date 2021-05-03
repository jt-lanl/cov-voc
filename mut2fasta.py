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
    paa("--reference","-r",required=True,
        help="input fasta file with reference sequence (eg Wuhan form)")
    paa("--includeref",action="store_true",
        help="include reference sequence in fasta output file")
    paa("--mutstrings","-m",nargs='+',
        help="one or more strings of the form '[L5F,...,Q957R]'")
    paa("--mutfile",
        help="input file with list of mutant strings")
    paa("--output","-o",
        help="ouptut fasta file of mutated seqeunces")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def rd_mutfile(filename):
    mutants=[]
    with open(filename) as f:
        for line in f:
            line = re.sub("#.*","",line)
            line = line.strip()
            if not line:
                continue
            mutants.append( mutant.Mutation(line) )
    return mutants

def main(args):

    mutants = []
    if args.mutfile:
        mutants.extend( rd_mutfile(args.mutfile) )
    if args.mutstrings:
        mutants.extend( mutant.Mutation(mstring)
                        for mstring in args.mutstrings)

    vprint(len(mutants),"mutant strings have been read")
    for m in mutants:
        vvprint("  ",m)

    if len(mutants) == 0:
        raise RuntimeError("Must specify mutants with either --mutfile for --mutstirngs or -m")

    ## Grab the first sequence from a reference file
    seqs = readseq.read_seqfile(args.reference,maxseqs=1)
    wuhan = seqs[0].seq

    if "-" in wuhan:
        warnings.warn("Dashes have been stripped from reference sequence")
        wuhan = re.sub("-","",wuhan)
        seqs[0].seq = wuhan

    mutseqs = []
    if args.includeref:
        mutseqs.append(seqs[0])
    for mut in mutants:
        refseq = list(wuhan)
        for ssm in mut:
            assert( refseq[ssm.site-1] == ssm.ref )
            refseq[ssm.site-1] = ssm.mut
        seq = "".join(refseq)
        mutseqs.append( readseq.SequenceSample( str(mut), seq ) )

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
    

