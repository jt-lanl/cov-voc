'''
Read mutant strings from a file.
Identify sequences from input sequence file that match those stringa.
Get the ISL numbers corresponding to those matching sequences.
Open /another/ file (eg a DNA file instead of protein).
Find all the sequences in that file with the given ISL numbers.
Determine the most common sequences among the matches.
Output an example (with mstring pattern + ISL number)
of a seqeunce that exhibits that most common variant.
Final output is a fasta file with a sequence for each mutant string
'''
## Took about 45 minutes to do 274 rows
## Using --skipx, or about 10M AA sequences, and 16M DNA sequences
## About 22 minutes with parallel -k 50


import re
from collections import Counter,namedtuple
import argparse

import warnings

import verbose as v
import breakpipe
import xopen
import sequtil
import covid
import mutant
import mstringfix

def _getargs():
    '''get arguments from command line'''
    aparser = argparse.ArgumentParser(description=__doc__)
    paa = aparser.add_argument
    covid.corona_args(aparser)
    paa("--mutfile","-M",
        help="file with list of mutant strings")
    paa("--dnainput",
        help="sequence file with reference (eg, DNA) sequences")
    paa("--tweakfile",
        help="use file for tweaking mstrings")
    paa("--output","-o",
        help="output fasta file with reference sequences")
    paa("--isloutput",
        help="write tsv table of lineage names and ISL names")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = aparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

NamedMutation = namedtuple("NamedMutation",['nom','mut'])

def read_mutantfile(filename):
    '''read mutant file and return mstring list'''
    ## assume tab separated columns
    nom_mut_list = []
    with xopen.xopen(filename) as fptr:
        for line in xopen.nonempty_lines(fptr):
            try:
                nom,mstring = line.split('\t')
                nom = re.sub(r' ','',nom)
                mstring = mstringfix.mstring_brackets(mstring)
                nom_mut_list.append(NamedMutation(nom,mstring))

            except ValueError:
                warnings.warn(f"Invalid line: {line}")
                continue

    return nom_mut_list

def read_fix_mutantfile(mutfile,tweakfile=None):
    '''read mutant file and make minor fixes
    (plus any fixes suggested by the tweakfile)
    '''
    nom_mut_list = read_mutantfile(mutfile)

    fixer = mstringfix.MStringFixer(tweakfile)
    fixer.append(r'\s*ancestral\s*','')
    fixer.append(r'G142[GD_]','G142.')
    v.vprint(f"FIXER:\n{fixer}")
    return [NamedMutation(nom_mut.nom,fixer.fix(nom_mut.mut))
            for nom_mut in nom_mut_list]

def get_refseqdict(dnainputfile,isl_matches):
    '''read DNA file, make a dict of sequences indexed by ISL number,
    only including those in the isl_matches dictionary'''
    isl_setofall = set()
    for islmatches in isl_matches.values():
        isl_setofall.update(islmatches)

    refseqs = sequtil.read_seqfile(dnainputfile)
    firstref,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    refseqdict = dict()
    for s in refseqs:
        isl_name = covid.get_isl(s)
        if isl_name in isl_setofall:
            refseqdict[isl_name] = s
    v.vprint("Read",len(refseqdict),"reference sequences")
    return firstref,refseqdict

def get_isl_matches(seqs,nom_mut_list,fmt_mstring=None):
    '''return a dict with a list of ISL numbers 
    corresponding to seqs that match the mstring pattern'''

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)
    seqs = list(seqs)

    isl_matches = dict() ## list of isl names for seq's that match pattern
    for nom,mstring in nom_mut_list:
        mregex = m_mgr.regex_from_mstring(mstring,exact=True)
        re_regex = re.compile(mregex)
        islnames = [covid.get_isl(s) for s in seqs
                    if re_regex.match(s.seq)]
        islnames = sorted(islnames,
                          key=lambda isl: int(re.sub('EPI_ISL_','',isl)))
        if len(islnames)==0:
            warnings.warn(f"No matches found for {nom}: {mstring}")
        v.vprint(fmt_mstring[mstring] if fmt_mstring else mstring,
                 " ".join(islnames[:3]),
                 "..." if len(islnames)>3 else "")
        isl_matches[mstring] = islnames
    v.vprint('Matched',sum(len(islnames) for islnames in isl_matches.values()),
             'ISL names in',len(isl_matches),'mstrings')
    return isl_matches

RE_DASH = re.compile('-')

@breakpipe.no_broken_pipe
def _main(args):
    '''mutlineage main'''

    nom_mut_list = read_fix_mutantfile(args.mutfile,args.tweakfile)

    max_len_mstring = min(50,max(len(nm.mut) for nm in nom_mut_list))
    fmt_mstring = {mstring: f'{mstring:>{max_len_mstring}s}'
                   for _,mstring in nom_mut_list}
    for nom,mstring in nom_mut_list:
        v.vprint(fmt_mstring[mstring],nom)

    seqs = covid.read_filter_seqfile(args)
    isl_matches = get_isl_matches(seqs,nom_mut_list,fmt_mstring)

    if not args.dnainput:
        v.print('DNA input file not specified')
        return

    firstref,refseqdict = get_refseqdict(args.dnainput,isl_matches)
    if not refseqdict:
        v.print(f'No matching DNA sequences in {args.dnainput}')
        return

    dna_seqs = []
    if firstref:
        dna_seqs.append(firstref)

    isloutputlines = []

    for nom,mstring in nom_mut_list:
        islnames = isl_matches[mstring]
        matchseqs = [refseqdict[islname]
                     for islname in islnames
                     if islname in refseqdict]
        v.vprint('islnames:',len(islnames),islnames[:5])
        v.vprint('matchseq:',len(matchseqs))

        try:
            [(cseq,_)] = Counter(RE_DASH.sub('',s.seq)
                                 for s in matchseqs).most_common(1)
        except ValueError:
            v.print(f"ISL not found for {nom}; mstring={mstring}")
            ## append dummy seq ("xxx")
            dna_seqs.append( sequtil.SequenceSample(nom,'xxx') )
            continue

        for islname in islnames:
            if islname not in refseqdict:
                continue
            if cseq == RE_DASH.sub('',refseqdict[islname].seq):
                v.print(fmt_mstring[mstring],islname,nom)
                sname = f'{nom}__{islname}'
                dna_seqs.append( sequtil.SequenceSample(sname,refseqdict[islname].seq) )
                isloutputlines.append(f'{nom}\t{islname}\t{sname}')
                break ## just grab the first (lowest) one

    if args.isloutput:
        with open(args.isloutput,'w') as fp_isl:
            print("\n".join(isloutputlines),file=fp_isl)

    if args.output:
        sequtil.write_seqfile(args.output,dna_seqs)



if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
