'''
"fix" input fasta (or tbl, mase, etc) sequence file by (optionally)
 - removing columns with dashes in reference sequence,
 - removing columns with dashes in all sequences,
 - padding sequences with dashes so all sequences are the same length
 - stripping final stop codons from each sequence,
 - applying a date filter
 - applying any name-based filter (eg, country or lineage)
 - codon aligning DNA sequence
 - translating DNA sequence to amino-acids
 - snipping out a specified range of sites
and writing the sequenecs to an output fasta file
'''

import re
from pathlib import Path
import itertools as it
import random
import argparse

import verbose as v
import readseq
import sequtil
import wrapgen
import intlist
import covid

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("--badisls",
        help="File with list of EPI_ISL numbers to exclude")
    paa("--translate",action="store_true",
        help="Translate from nucleotides to amino acids")
    paa("--codonalign",action="store_true",
        help="Add dashes to keep seqeunces aligned in triples")
    paa("--padlength",action="store_true",
        help="Pad seqs with dashes so all are the same length")
    paa("--rmgapcols",action="store_true",
        help="Remove gap-only columns")
    paa("--snipsites",
        help="List of sites (eg 1-4,7,9-240) to be sent to output")
    paa("--sniprange",type=int,nargs=2,
        ## like snipsites, less flexible, but maybe faster?
        help="site range given as two integers, lo and hi")
    paa("--rmdash",action="store_true",
        help="remove all the dashes in each sequence")
    paa("--splitout",
        help="split the output into a separate file for each sequence")
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

def codon_align_indices(refseq):
    '''return indices that require appending a dash'''
    ## assume refseq contains triples of non-dash characters (eg, TAA, CGT, etc)
    ## interspersed with dashes
    ## if refseq constains strings like "-TA-", this violates assumption
    ## return list of indices for which adding a dash will lead to triples of dashes as well
    ndxlist = []
    ndx = 0
    while ndx < len(refseq):
        if refseq[ndx] != '-':
            ## we are assuming refseq[ndx+1] != '-'
            ##             and refseq[ndx+2] != '-'
            ndx += 3
        elif refseq[ndx:ndx+3] == '---':
            ndx += 3
        elif refseq[ndx:ndx+2] == '--':
            ndxlist.append(ndx)
            ndx += 2
        else:
            ## if single dash at this site, then
            ## append two copies of ndx, so two dashes will be appended
            ndxlist.append(ndx)
            ndxlist.append(ndx)
            ndx += 1
    print('ndxlist:',ndxlist)
    return ndxlist

def codon_align_seqs(seqs,ndxlist=None):
    '''replace sequences with codon-aligned sequences'''
    if ndxlist is None:
        first,seqs = sequtil.get_first_item(seqs)
        ndxlist = codon_align_indices(first.seq)
    for s in seqs:
        slist = list(s.seq)
        for ndx in ndxlist:
            slist[ndx] += "-"
        s.seq = "".join(slist)
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

def get_max_length(seqs):
    '''get length of longest sequence'''
    return max(len(s.seq) for s in seqs)

def pad_to_length(seqs,length=None):
    '''pads all sequences to a common length'''
    seqs=list(seqs) ## will consume seqs if iterator
    if not length:
        length = get_max_length(seqs)
    for s in seqs:
        s.seq = s.seq + "-"*(length-len(s.seq))
    return seqs

def snipsites(seqs,sitestring,offset=1):
    '''replace sequence with shorter sequence including only sites in list'''
    ## offset=1 for 1-based counting, eg 1'st site is 1
    sitelist = intlist.string_to_intlist(sitestring)
    for s in seqs:
        s.seq = "".join(s.seq[offset+n] for n in sitelist)
        yield s

def sniprange(seqs,lohi,offset=1):
    '''replace sequence with shorter sequence in range lo:hi+1'''
    v.vprint('lohi:',lohi)
    lo,hi = lohi
    lo,hi = lo+offset,hi+offset
    for s in seqs:
        s.seq = s.seq[lo:hi+1]
        yield s

def rmdashes(seqs):
    '''remove dashes in sequences'''
    for s in seqs:
        s.seq = re.sub("-","",s.seq)
        yield s

def main(args):
    '''fixfasta main'''
    v.vprint(args)

    def vcount(seqs,*p,**kw):
        "count items in the generator as they go by"
        return wrapgen.keepcount(seqs,*p,**kw) if args.verbose else seqs

    seqs = covid.read_filter_seqfile(args)

    if args.badisls:
        bads = getisls(args.badisls)
        seqs = (s for s in seqs
                   if not any( b in s.name for b in bads ))
        seqs = vcount(seqs,"Sequences after removing bad ISLs:")

    if args.random:
        seqs = list(seqs)
        seqs = seqs[:1] + random.sample(seqs[1:],k=len(seqs[1:]))

    if args.codonalign:
        seqs = codon_align_seqs(seqs)
        seqs = vcount(seqs,"Sequences codon aligned")

    if args.translate:
        seqs = translate_to_aa_alt(seqs)
        seqs = vcount(seqs,"Sequences translated")

    if args.padlength:
        seqs = pad_to_length(seqs)

    if args.snipsites:
        seqs = snipsites(seqs,args.snipsites)

    if args.sniprange:
        seqs = sniprange(seqs,args.sniprange)

    if args.rmgapcols:
        first,seqs = sequtil.get_first_item(seqs)
        initlen = len(first.seq)
        seqs = sequtil.remove_gap_columns(seqs)
        first,seqs = sequtil.get_first_item(seqs)
        v.vprint(f"Removed {initlen-len(first.seq)} dashes")

    if args.rmdash:
        seqs = rmdashes(seqs)

    if args.output:
        readseq.write_seqfile(args.output,seqs)

    if args.splitout:
        if "." in args.splitout:
            base,ext = args.splitout.split(".")
        else:
            base = args.splitout
            ext = ".fasta"

        for n,s in enumerate(seqs,start=1):
            readseq.write_seqfile("%s-%03d.%s" % (base,n,ext), [s])

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
