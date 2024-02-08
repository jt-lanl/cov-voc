'''
Read mutant strings from a file.
Identify sequences from input sequence file that match those stringa.
Get the ISL numbers corresponding to those matching sequences.
Open /another/ file (eg a DNA file instead of protein).
Find all the sequences in that file with the given ISL numbers.
Determine the most common sequences among the matches.
Output an example (with mstring pattern + sequence name, including ISL number)
of a seqeunce that exhibits that most common variant.
Final output is a fasta file with a sequence for each mutant string
'''
import re
from collections import Counter
import argparse

import warnings

import sequtil
import covid
import verbose as v
import mutant
import pseq
import mstringfix

def _getargs():
    '''get arguments from command line'''
    aparser = argparse.ArgumentParser(description=__doc__)
    paa = aparser.add_argument
    covid.corona_args(aparser)
    paa("--mutfile","-M",
        help="file with list of mutant strings")
    paa("-j",
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
    return args

def read_mutantfile(filename):
    '''read mutant file and return mstring list'''
    ## assume tab separated columns
    mutlist=[]
    nomlist=[]
    with open(filename) as fptr:
        for line in fptr:
            line = line.strip()
            line = re.sub('#.*','',line)
            if not line:
                continue
            try:
                nom,mstring = line.split('\t')

                mstring = mstringfix.mstring_brackets(mstring)
                mutlist.append(mstring)

                nom = nom.strip()
                nom = re.sub(r' ','',nom)
                nomlist.append(nom)

            except ValueError:
                warnings.warn(f"Invalid line: {line}")
                continue

    return nomlist,mutlist

def _main(args):
    '''mutlineage main'''

    nomlist,mutlist = read_mutantfile(args.mutfile)

    fixer = mstringfix.MStringFixer(args.tweakfile)
    fixer.append(r'\s*ancestral\s*','')
    fixer.append(r'G142[GD_]','G142.')
    v.vprint(f"FIXER:\n{fixer}")
    mutlist = [fixer.fix(mstring) for mstring in mutlist]


    mutfmt = "%%%ds" % max(len(mstring) for mstring in mutlist)
    for nom,mstring in zip(nomlist,mutlist):
        v.vprint(mutfmt % mstring,nom)

    seqs = covid.read_filter_seqfile(args)

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)

    #seqs = list(seqs)
    seqs = [pseq.ProcessedSequence(m_mgr,s) for s in seqs]

    isl_matches = dict() ## list of isl names for seq's that match pattern
    isl_setofall = set() ## set of all islnames
    for nom,mstring in zip(nomlist,mutlist):
        mpatt = mutant.Mutation.from_mstring(mstring,exact=True)
        #matches = m_mgr.filter_seqs_by_pattern(mpatt,seqs)
        matches = pseq.filter_pseqs_by_pattern(m_mgr,mpatt,seqs,exact=True)
        islnames = [s.ISL for s in matches]
        islnames = sorted(islnames,
                          key=lambda isl: int(re.sub('EPI_ISL_','',isl)))
        if len(islnames)==0:
            warnings.warn(f"No matches found for {nom}: {mstring}")
        v.vprint(mutfmt % mstring," ".join(islnames[:3]),
               "..." if len(islnames)>3 else "")
        isl_matches[mstring] = islnames
        isl_setofall.update(islnames)

    if not args.j:
        return

    refseqs = sequtil.read_seqfile(args.j)
    firstref,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    refseqdict = dict()
    for s in refseqs:
        isl_name = covid.get_isl(s.name)
        if isl_name in isl_setofall:
            refseqdict[isl_name] = s
    v.vprint("Read",len(refseqdict),"reference sequences")
    if not refseqdict:
        return

    outseqs = []
    if firstref:
        outseqs.append(firstref)

    if args.isloutput:
        fp_isl = open(args.isloutput,'w')

    for nom,mstring in zip(nomlist,mutlist):
        islnames = isl_matches[mstring]
        matchseqs = (refseqdict.get(islname,None)
                     for islname in islnames )
        matchseqs = filter(lambda x: x is not None,matchseqs)
        matchseqs = list(matchseqs)
        v.vprint('islnames:',len(islnames),islnames[:5])
        v.vprint('matchseq:',len(matchseqs))
                 
        try:
            [(cseq,_)] = Counter(re.sub('-','',s.seq)
                                 for s in matchseqs).most_common(1)
        except ValueError:
            print(mutfmt % mstring,f"ISL not found for {nom}")
            outseq = sequtil.SequenceSample(nom,'xxx')
            outseqs.append(outseq)
            continue

        for islname in islnames:
            if islname not in refseqdict:
                continue
            if cseq == re.sub('-','',refseqdict[islname].seq):
                print(mutfmt % mstring,islname,nom)
                sseq = refseqdict[islname].seq
                sname = f'{nom}__{islname}'
                outseq = sequtil.SequenceSample(sname,sseq)
                outseqs.append(outseq)
                break ## just grab the first (lowest) one
        if args.isloutput:
            print(f'{nom}\t{islnames[0]}\t{islname}',file=fp_isl)

    if args.isloutput:
        fp_isl.close()

    if args.output:
        sequtil.write_seqfile(args.output,outseqs)



if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
