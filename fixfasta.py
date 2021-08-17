'''
"fix" input fasta (or tbl) sequence file by
1/ removing columns with dashes in reference sequence,
2/ stripping final stop codons from each sequence,
3/ applying a date filter
'''
import sys
import re
from pathlib import Path
import itertools as it
import random
import argparse

import readseq
import wrapgen
import covid

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("-N",type=int,default=0,
        help="keep at most this many sequences")
    paa("--badisls",
        help="File with list of EPI_ISL numbers to exclude")
    paa("--translate",action="store_true",
        help="Translate from nucleotides to amino acids")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

codon_to_aa_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'X', 'TAG':'X',
    'TGC':'C', 'TGT':'C', 'TGA':'X', 'TGG':'W',
    '---':'-',
}

def translate_to_aa(seqs):
    ''' for a list of nucleotide sequences, yield amino acid sequences '''
    for s in seqs:
        aalist = []
        for j in range(0,len(s.seq)-2,3):
            codon = s.seq[j:j+3]
            aa = codon_to_aa_table.get(codon,"X")
            aalist.append( aa )
        s.seq = "".join(aalist)
        yield s

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)

def translate_to_aa_alt(seqs):
    ''' for a list of nucleotide sequences, yield amino acid sequences '''
    codon_to_aa_table_alt = dict()
    for c,aa in codon_to_aa_table.items():
        codon_to_aa_table_alt[(c[0],c[1],c[2])] = aa
    for s in seqs:
        codonlist = grouper(s.seq,3,'-')
        #aalist = [codon_to_aa_table.get("".join(codon),"X") for codon in codonlist]
        aalist = [codon_to_aa_table_alt.get(codontuple,"X") for codontuple in codonlist]
        s.seq = "".join(aalist[:-1])  ## don't keep the stop codon == X
        yield s

def getisls(file):
    '''read list of ISL numbers from text file'''
    isls=[]
    with open(file,"r") as isl_file:
        for line in isl_file:
            m = re.match(r".*(EPI_ISL_\d+).*",line.strip())
            if m:
                isls.append(m[1])
    return isls

def main(args):
    '''fixfasta main'''
    seqs = covid.read_filter_seqfile(args)

    if args.badisls:
        bads = getisls(args.badisls)
        seqs = (s for s in seqs
                   if not any( b in s.name for b in bads ))
        seqs = vcount(seqs,"Sequences after removing bad ISLs:")

    if args.random:
        seqs = list(seqs)
        seqs = seqs[:1] + random.sample(seqs[1:],k=len(seqs[1:]))

    if args.N:
        seqs = it.islice(seqs,args.N)
        seqs = vcount(seqs,"Sequences after truncation:")

    if args.translate:
        seqs = translate_to_aa_alt(seqs)
        seqs = vcount(seqs,"Sequences translated")

    if args.output:
        readseq.write_seqfile(args.output,seqs)

if __name__ == "__main__":

    _args = getargs()
    def vprint(*p,**kw):
        "verbose print"
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        "very verbose print"
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    def vcount(seqs,*p,**kw):
        "count items in the generator as they go by"
        return wrapgen.keepcount(seqs,*p,**kw) if _args.verbose else seqs

    main(_args)
