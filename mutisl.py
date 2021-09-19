'''
read mutant strings from a file
identify sequences from input sequence file that match those stringa
get the ISL numbers corresponding to those matching sequences
Open /another/ file (eg a DNA file instead of protein)
find all the sequences in that file with the given IDL numbers,
determint the most common sequences among the matches
output an example (with mstring pattern + sequence name, including ISL number)
of a seqeunce that exhibits that most common variant
Final output is a fasta file with a sequence for each mutant string
'''
import sys
import re
from collections import Counter
import itertools as it
import argparse

import warnings

import sequtil
import covid
import mutant

def _getargs():
    '''get arguments from command line'''
    aparser = argparse.ArgumentParser(description=__doc__)
    paa = aparser.add_argument
    covid.corona_args(aparser)
    paa("--mutfile","-M",
        help="file with list of mutant strings")
    paa("-j",
        help="sequence file with reference (eg, DNA) sequences")
    paa("-N",type=int,default=0,
        help="just look at the first N sequences in input file")
    paa("--output","-o",
        help="output fasta file with reference sequences")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = aparser.parse_args()
    return args

def mstring_brackets(mstring):
    '''make sure mstring has brackets around it'''
    mstring = re.sub(r'\s','',mstring) ## remove extra spaces between ssms
    mstring = re.sub(r'[\[\]]','',mstring) ## remove if already there
    mstring = f'[{mstring}]'
    return mstring

def read_mutantfile(filename):
    '''read mutant file and return mstring list'''
    mutantfilepatt = re.compile(r'(.*)\s+(\[?.+\]?)')
    mutlist=[]
    nomlist=[]
    with open(filename) as fptr:
        for line in fptr:
            line = line.strip()
            if not line:
                continue
            mpatt = mutantfilepatt.search(line)
            if mpatt:
                mstring = mstring_brackets(mpatt[2])
                mutlist.append(mstring)
                nomlist.append(mpatt[1])
            else:
                warnings.warn(f"Invalid line: {line}")
    return nomlist,mutlist

def _main(args):
    '''mutlineage main'''

    nomlist,mutlist = read_mutantfile(args.mutfile)
    mutfmt = "%%%ds" % max(len(mstring) for mstring in mutlist)
    for nom,mstring in zip(nomlist,mutlist):
        vprint(mutfmt % mstring,nom)

    seqs = covid.read_filter_seqfile(args)
    if args.N:
        ## just grab the first N (for debugging w/ shorter runs)
        seqs = it.islice(seqs,args.N+1)

    seqs = list(seqs)
    first,seqs = sequtil.get_first_item(seqs)
    m_mgr = mutant.MutationManager(first.seq)
    isl_matches = dict() ## list of isl names for seq's that match pattern
    isl_setofall = set() ## set of all islnames
    for mstring in mutlist:
        mpatt = mutant.Mutation.from_mstring(mstring,exact=True)
        matches = m_mgr.filter_seqs_by_pattern(mpatt,seqs)
        islnames = [covid.get_isl(s.name) for s in matches]
        islnames = sorted(islnames)
        if len(islnames)==0:
            warnings.warn(f"No matches found for {mstring}")
        vprint(mutfmt % mstring," ".join(islnames[:3]),
               "..." if len(islnames)>3 else "")
        isl_matches[mstring] = islnames
        isl_setofall.update(islnames)

    if not args.j:
        return
        
    refseqs = sequtil.read_seqfile(args.j)
    firstref,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    refseqdict = {covid.get_isl(s.name): s
                  for s in refseqs
                  if covid.get_isl(s.name) in isl_setofall}
    vprint("Read",len(refseqdict),"reference sequences")

    outseqs = []
    if firstref:
        outseqs.append(firstref)

    for nom,mstring in zip(nomlist,mutlist):
        islnames = isl_matches[mstring]
        if refseqdict:
            matchseqs = (refseqdict.get(islname,None)
                         for islname in islnames )
            matchseqs = filter(lambda x: x is not None,matchseqs)
            [(cseq,_)] = Counter(s.seq for s in matchseqs).most_common(1)
            for islname in islnames:
                if islname not in refseqdict:
                    continue
                if cseq == refseqdict[islname].seq:
                    print(mutfmt % mstring,islname)
                    s = refseqdict[islname]
                    sname = ":".join([nom,s.name,mstring])
                    outseq = sequtil.SequenceSample(sname,s.seq)
                    outseqs.append(outseq)
                    break ## just grab the first one


    if args.output:
        sequtil.write_seqfile(args.output,outseqs)            
            

if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    _main(_args)
