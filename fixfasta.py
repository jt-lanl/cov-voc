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

import warnings

import verbose as v
import breakpipe
from xopen import xopen
from lineagenotes import LineageNotes
import sequtil
import wrapgen
import mutant
import covid

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    paa("--jobno",type=int,default=1,
        help="when running multiple jobs in parallel, use --jobno {#}")
    covid.corona_args(ap)
    paa = ap.add_argument_group("Fix Sequence(fasta) Options").add_argument
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
    paa("--stripdashcols",action="store_true",
        help="Strip columns with dash in reference sequence")
    paa("--islreplace",
        help="file of sequences used to replace existing seqs based on ISL")
    paa("--keepsites",
        help="List of sites (eg 1-4,7,9-240) to be sent to output")
    paa("--pangoreplace",
        help="Name of file with pango name for each sequence")
    paa("--fclades",nargs='+',
        help="Filter sequences to include those in the fclade(s)")
    paa("--xclades",nargs='+',
        help="Filter sequences to exclude those in the xclade(s)")
    paa("--rmgapcols",action="store_true",
        help="Remove gap-only columns")
    paa("--notesfile", help="lineage_notes.txt")
    paa("--keyfile",help="alias_key.json")
    paa("--random",action="store_true",
        help="randomize input data order")
    paa("--reverse",action="store_true",
        help="reverse input data order")
    paa("--sortbydate",action="store_true",
        help="sort input sequences by date")
    paa("--sortbyisl",action="store_true",
        help="sort input sequences by ISL number")
    paa("--output","-o",type=Path,
        help="output fasta file")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    if args.codonalign and (args.rmgapcols or args.stripdashcols):
        v.print(args)
        raise RuntimeError('If --codealign, then cannot have --rmgapcols or --stripdashcols')
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

def grouper(iterable, chunksize, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * chunksize
    return it.zip_longest(*args, fillvalue=fillvalue)

def translate_to_aa_alt(seqs):
    ''' for a list of nucleotide sequences, yield amino acid sequences '''
    codon_to_aa_table_alt = dict()
    for nuc,aa in codon_to_aa_table.items():
        codon_to_aa_table_alt[(nuc[0],nuc[1],nuc[2])] = aa
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
    v.vprint('codon align ndxlist:',ndxlist)
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
            mat = re.match(r".*(EPI_ISL_\d+).*",line.strip())
            if mat:
                isls.append(mat[1])
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

def get_pangoreplace(pangofile):
    '''Use USHER-provided pango names
    Args:
       pangofile (Path or str): name of file with pango names
    Returns:
       pangocat (dict): catalog maps ISL numbers to pango names'''
    ## lower-case (possibly preceded by an underscore) indicates
    ## non-standard pango name (eg, proposed, misc, dropout, _rev, _quad, etc)
    ## in that case, removed the nonstandard part of the name
    RE_NONPANGO=re.compile('_?[a-z_].*')
    pangocat = dict()
    if not pangofile:
        return pangocat
    with xopen(pangofile,"r") as fpin:
        non_std_count = 0
        total_count = 0
        for line in fpin:
            tokens = line.strip().split()
            if len(tokens) != 2 or tokens[0]=="name":
                continue
            isl = covid.EPI_ISL_REGEX.search(tokens[0])
            if not isl:
                v.print_only(3,'badisl:',f'in seqname={tokens[0]}')
            isl = isl.group(0)
            pango = tokens[1]
            pango = RE_NONPANGO.sub('',pango)
            total_count += 1
            if not pango:
                non_std_count += 1
                continue
            pangocat[isl]=pango
        if non_std_count:
            v.vprint(f'Nonstandard pango names '
                     f'in usher catalog: {non_std_count}/{total_count}')
    return pangocat

def filterclades(seqs,lin_notes,cladelist,exclude=False):
    '''keep (or exclude) seqs whose lineage is in one of the clades'''
    cladelist = [lin_notes.get_fullname(clade)
                 for clade in cladelist or []]
    v.vvvprint('exclude:',exclude,'Clades:',cladelist)
    if not cladelist:
        yield from seqs
    else:
        for s in seqs:
            lin = covid.get_lineage(s)
            flin = lin_notes.get_fullname(lin)
            lin_in_clades = any(flin.startswith(clade) for clade in cladelist)
            v.vvprint(f'{exclude=},{lin=},{flin=},{lin_in_clades=}')
            if ((lin_in_clades and not exclude) or
                (not lin_in_clades and exclude)):
                yield s

def apply_pangocat(seqs,pangocat):
    '''update sequences lineages using pangocat'''
    if not pangocat:
        yield from seqs
    for s in seqs:
        isl = covid.get_isl(s)
        lin = covid.get_lineage(s)
        new_lin = pangocat.get(isl,lin)
        if new_lin != lin:
            v.vprint_only(3,'New Pango',f'{lin} -> {new_lin}')
            s.name = re.sub(lin+r'$',new_lin,s.name)
        yield s

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

def sort_bydate(seqs):
    '''re-order seqs by date'''
    seqs = list(seqs).sort(key=covid.get_date)
    return seqs

    

        
def filter_all_seqs(args,seqs):
    '''filter all sequences, even the ref sequence'''
    if args.codonalign:
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

    return seqs

def filter_nonref_seqs(args,seqs):
    '''filter all sequences, but assume ref seq is taken out'''

    if args.pangoreplace:
        pangocat = get_pangoreplace(args.pangoreplace)
        seqs = apply_pangocat(seqs,pangocat)

    if args.islreplace:
        isl_seqs = seqs_indexed_by_isl(args.islreplace)
        seqs = islreplace(isl_seqs,seqs)

    if args.badisls:
        bads = set(isls_from_file(args.badisls))
        seqs = (s for s in seqs
                if covid.get_isl(s) not in bads)

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

    return seqs

@breakpipe.no_broken_pipe
def main(args):
    '''fixfasta main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)

    ## First set of filters here, we want to apply to ALL sequences
    ## including the first (which, for now, is included in the seqs array)

    seqs = filter_all_seqs(args,seqs)

    ## next batch of filters are not meant to be applied
    ## to the reference sequence, so take it out
    first,seqs = covid.get_first_item_ifref(seqs)

    seqs = filter_nonref_seqs(args,seqs)

    if args.stripdashcols:
        if not first:
            raise RuntimeError('Cannot --stripdashcols becuz first sequence of '
                               f'input={args.input} is not a reference sequence')
        ## removes all the columns that have dashes in first.seq
        ## actually, we /do/ want first among the seqs here
        seqs = it.chain([first],seqs)
        seqs = sequtil.stripdashcols(first.seq,seqs)
        ## and now we take the first back out again...
        first,seqs = covid.get_first_item(seqs,keepfirst=False)

    if args.fclades or args.xclades:
        lin_notes = LineageNotes.from_file(args.notesfile,args.keyfile,fix=True)
        seqs = filterclades(seqs,lin_notes,args.fclades)
        seqs = filterclades(seqs,lin_notes,args.xclades,exclude=True)

    ## these next filters will require that the FULL list of sequences
    ## be available. so we cannot be parallelizing the sequence list
    ## which means that jobno != 1 will be an error
    ## if jobno==1, then the first has already been taken out of the seqs

    if args.random:
        assert args.jobno == 1
        seqs = list(seqs)
        v.vprint("Randomizing sequence order...",end="")
        seqs = random.sample(seqs,k=len(seqs))
        v.vprint("ok")

    if args.reverse:
        assert args.jobno == 1
        v.vprint("Reversing sequences order ..",end="")
        seqs = reversed(list(seqs))
        v.vprint("ok")

    if args.sortbydate:
        assert args.jobno == 1
        v.vprint("Sorting sequences by date...",end="")
        seqs = list(seqs)
        seqs.sort(key=lambda s:
                  covid.get_date(s,as_string=True,check_date=False))
        v.vprint("ok")

    if args.sortbyisl:
        assert args.jobno == 1
        v.vprint("Sorting sequences by ISL...",end="")
        seqs = list(seqs)
        seqs.sort(key = covid.get_isl_number)
        v.vprint("ok")

    if args.rmgapcols:
        v.vprint("Removing columns that are all-gaps...",end="")
        assert args.jobno == 1
        if not first:
            ## it is not required for first seq to be reference seq
            ## but if it's not, something might be wrong
            warnings.warn('first seqeunce is not a reference sequence')
            ## we'll take the first sequence anyway!
            first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
        initlen = len(first.seq)
        first,seqs = sequtil.remove_gap_columns(first,seqs)
        v.vprint(f"ok. Removed {initlen-len(first.seq)} dashes")

    if args.output:
        if args.jobno == 1 and first is not None:
            ## put that first seq back in
            seqs = it.chain([first],seqs)
        sequtil.write_seqfile(args.output,seqs)

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    def vcount(seqs,*p,**kw):
        "count items in the generator as they go by"
        return wrapgen.keepcount(seqs,*p,**kw) if _args.verbose else seqs

    main(_args)
