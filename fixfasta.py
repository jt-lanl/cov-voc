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
import sequtil
import wrapgen
import mutant
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
    paa("--keepisls",
        help="File with list of EPI_ISL numbers to keep")
    paa("--translate",action="store_true",
        help="Translate from nucleotides to amino acids")
    paa("--codonalign",action="store_true",
        help="Add dashes to keep seqeunces aligned in triples")
    paa("--padlength",action="store_true",
        help="Pad seqs with dashes so all are the same length (as longest seq)")
    paa("--toomanyx",type=int,default=0,
        help="Remove sequences that have too many total X's")
    paa("--toomanygaps",type=int,default=0,
        help="Remove sequences with too many gaps in a single stretch")
    paa("--rmdash",action="store_true",
        help="remove all the dashes in each sequence")
    paa("--rmgapcols",action="store_true",
        help="Remove gap-only columns")
    paa("--stripdashcols",action="store_true",
        help="Strip columns with dash in reference sequence")
    paa("--islreplace",
        help="file of sequences used to replace existing seqs based on ISL")
    paa("--keepsites",
        help="List of sites (eg 1-4,7,9-240) to be sent to output")
    paa("--output","-o",type=Path,
        help="output fasta file")
    paa("--jobno",type=int,default=1,
        help="when running multiple jobs in parallel, use --jobno {#}")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
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

def isls_from_file(file):
    '''read list of ISL numbers from text file'''
    isls=[]
    with open(file,"r") as isl_file:
        for line in isl_file:
            m = re.match(r".*(EPI_ISL_\d+).*",line.strip())
            if m:
                isls.append(m[1])
    v.vprint(f'ISL file {file} has {len(set(isls))} distinct ISL numbers')
    return isls

def seqs_indexed_by_isl(file):
    '''read file and make dict of sequences indexed by ISL number'''
    seqdict = dict()
    seqs = sequtil.read_seqfile(file,badchar='X')
    for s in seqs:
        isl = covid.get_isl(s.name)
        seqdict[isl] = s
    return seqdict

def islreplace(seqdict,seqs):
    '''filter seqs, replacing any that are in seqdict with the seqdict value'''
    for s in seqs:
        isl = covid.get_isl(s.name)
        if isl in seqdict:
            lin_a = covid.get_lineage(s.name)
            lin_b = covid.get_lineage(seqdict[isl].name)
            v.vprint(f"updated ISL: {isl} {lin_a} -> {lin_b}")
            v.vprint(f"             {s.name}")
            v.vprint(f"         ->  {seqdict[isl].name}")
        yield seqdict.get(isl,s)

def pad_to_length(seqs,length=None):
    '''pads all sequences to a common length'''
    if not length:
        seqs=list(seqs) ## will consume seqs if iterator
        length = max(len(s.seq) for s in seqs)
    for s in seqs:
        if len(s.seq) > length:
            ## if seq too long, truncate
            s.seq = s.seq[:length]
        elif len(s.seq) < length:
            s.seq = s.seq + "-"*(length-len(s.seq))
        yield s


def rmdashes(seqs):
    '''remove dashes in sequences'''
    for s in seqs:
        s.seq = re.sub("-","",s.seq)
        yield s

def rm_toomanygaps(toomany,seqs):
    '''filter out sequnces that have too many gaps'''
    count_removed=0
    if toomany == 0:
        yield from seqs
    else:
        biggap = "-" * toomany
        bigN = "N" * toomany
        for s in seqs:
            if biggap in s.seq or bigN in s.seq:
                count_removed += 1
                continue
            yield s
    if count_removed:
        v.vprint(f'Removed {count_removed} sequences with gap > {toomany}')

def main(args):
    '''fixfasta main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)

    ## First set of filters here, we want to apply to ALL sequences
    ## including the first (which, for now, is included in the seqs array)

    if args.codonalign:
        v.print('Warning: need to keep first for this too...')
        seqs = codon_align_seqs(seqs)
        seqs = vcount(seqs,"Sequences codon aligned")

    if args.translate:
        v.print('Warning: need to also translate first')
        seqs = translate_to_aa_alt(seqs)
        seqs = vcount(seqs,"Sequences translated")

    if args.padlength:
        seqs = pad_to_length(seqs)

    if args.keepsites:
        first,seqs = covid.get_first_item(seqs,keepfirst=True)
        mut_mgr = mutant.MutationManager(first.seq)
        seqs = covid.keepsites(mut_mgr,seqs,args.keepsites)

    if args.rmdash:
        seqs = rmdashes(seqs)

    ## next batch of filters are not meant to be applied
    ## to the reference sequence, so take it out
    first,seqs = covid.get_first_item(seqs,keepfirst=False)

    if args.stripdashcols:
        ## removes all the columns that have dashes in first.seq
        ## actually, we /do/ want first among the seqs here
        seqs = it.chain([first],seqs)
        seqs = sequtil.stripdashcols(first.seq,seqs)
        ## and now we take the first back out again...
        first,seqs = covid.get_first_item(seqs,keepfirst=False)

    if args.islreplace:
        isl_seqs = seqs_indexed_by_isl(args.islreplace)
        seqs = islreplace(isl_seqs,seqs)

    if args.badisls:
        bads = isls_from_file(args.badisls)
        seqs = (s for s in seqs
                if not any( b in s.name for b in bads ))

    if args.keepisls:
        keepers = set(isls_from_file(args.keepisls))
        seqs = (s for s in seqs
                if covid.get_isl(s) in keepers)
        seqs = vcount(seqs,"Sequences found with specified ISLs:")

    if args.toomanygaps:
        seqs = rm_toomanygaps(args.toomanygaps,seqs)

    if args.toomanyx:
        seqs = (s for s in seqs
                if s.seq.count('X') < args.toomanyx)

    ## these next filters will require that we FULL list of sequence
    ## be available. so we cannot be parallelizing the sequence list
    ## which means that jobno != 1 will be an error
    ## if jobno==1, then the first has already been taken out of the seqs

    if args.random:
        v.vprint("Randomizing sequence order...",end="")
        assert args.jobno == 1
        seqs = list(seqs)
        seqs = random.sample(seqs,k=len(seqs))
        v.vprint("ok")

    if args.rmgapcols:
        v.vprint("Removing columns that are all-gaps...",end="")
        assert args.jobno == 1
        initlen = len(first.seq)
        first,seqs = sequtil.remove_gap_columns(first,seqs)
        v.vprint(f"ok. Removed {initlen-len(first.seq)} dashes")

    if args.output:
        if args.jobno == 1:
            ## we need to put that first seq back in
            seqs = it.chain([first],seqs)
        sequtil.write_seqfile(args.output,seqs)

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    def vcount(seqs,*p,**kw):
        "count items in the generator as they go by"
        return wrapgen.keepcount(seqs,*p,**kw) if _args.verbose else seqs

    main(_args)
