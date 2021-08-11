'''Check for (and fix) alignment consistency'''
import sys
import re
from collections import defaultdict,Counter
import argparse
from functools import lru_cache

from seqsample import SequenceSample
import readseq
import mutantx as mutant
import covid

def getargs():
    '''parse options from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--output","-o",
        help="output inconsistent sequences")
    paa("--fix",
        help="fix sequences and write output to this fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

dedash = re.compile("-")

@lru_cache(maxsize=None)
def de_gap(seq):
    '''remove '-'s from sequence'''
    return dedash.sub("",seq)

def mutant_diffs(mref,mseq):
    '''
    given input Mutation's mref and mseq,
    pull out comment elements, and just report the differences; eg
    input: mref=[H146-,D614G] and mseq=[Y144-,H146Y,D614G];
    output: <H146-> and <Y144-,H146Y>
    '''
    msetref = set(mref)
    msetseq = set(mseq)

    mutref = mutant.Mutation(list(msetref-msetseq)).sort()
    mutseq = mutant.Mutation(list(msetseq-msetref)).sort()
    summary = ("<"+",".join(str(ssm) for ssm in mutref)+">",
               "<"+",".join(str(ssm) for ssm in mutseq)+">")
    return summary

def choose_alignment(MM,seqs):
    '''
    given a set of sequences whose de-gapped versions are identical;
    choose one (the best one?)
    '''
    ## first, make sure that de-gapped sequences are identical
    sseqs = [de_gap(s.seq) for s in seqs]
    assert all(sseq == sseqs[0] for sseq in sseqs[1:])
    ## find forms nearest to refseq
    mutlen = {s: len(MM.get_mutation(s.seq))
              for s in seqs}
    minmutlen = min(mutlen.values())
    candidateseqs = [s for s in seqs if mutlen[s] == minmutlen]
    ## find most common form
    cnt = Counter(s.seq for s in candidateseqs)
    [(cseq,_)] = cnt.most_common(1)
    return cseq

def write_summary(mismatches,diffs):
    '''summarize inconsistent sequences'''
    ## mismatches, diffs both Counter()'s

    print("Long summary of differences:")
    for (g,b),count in mismatches.most_common():
        asterisk = " *" if g==b else "" ## indicates equivalent sequences
        print("%6d %s\n       %s%s" % (count,g,b,asterisk))
    print()
    print("Short summary of differences:")
    gmaxlen = max(len(g) for g,b in diffs)
    gfmt = "%%6d %%%ds vs %%s" % (gmaxlen,)
    for (g,b),count in diffs.most_common():
        print(gfmt % (count,g,b))

def write_output(outputfile,first,fix_table,shareshortseq,inconsistent_shortseqs):
    '''ouput fasta file showing all the inconsistent sequences'''
    if outputfile:
        iseqs = [first]
        for n,sseq in enumerate(inconsistent_shortseqs,start=1):
            fullseqs = shareshortseq[sseq]
            goodseq = fix_table.get(fullseqs[0].seq,fullseqs[0].seq)
            iseqs.append(SequenceSample(f"Case_{n}",
                                        "-"*len(goodseq)))
            iseqs.append(SequenceSample(f"BaseSequence_{n}",
                                        goodseq))
            for s in fullseqs:
                if s.seq != goodseq:
                    iseqs.append(s)

        readseq.write_seqfile(outputfile,iseqs)


def main(args):
    '''checkalignment main'''

    args.keepdashcols = True
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs)

    seqs=list(seqs)

    ## sharedshortseq[key] is a list of seqs whose de-gapped seqstring = key
    shareshortseq = defaultdict(list)
    for s in seqs[1:]:
        shareshortseq[de_gap(s.seq)].append(s)

    ## find the shortseqs with inconsistent long sequence strings
    inconsistent_shortseqs = set()
    for sseq in shareshortseq:
        fullseqs = shareshortseq[sseq]
        if len(fullseqs) == 1:
            continue
        if any(s.seq != fullseqs[0].seq for s in fullseqs[1:]):
            inconsistent_shortseqs.add(sseq)

    fix_table = dict()
    summarize_mismatches = Counter()
    summarize_diffs = Counter()
    MM = mutant.MutationMaker(first.seq)
    countbad = 0
    for sseq in inconsistent_shortseqs:
        fullseqs = shareshortseq[sseq]
        goodseq = choose_alignment(MM,fullseqs)
        goodmut = MM.get_mutation(goodseq)
        badseqs = [s for s in fullseqs if s.seq != goodseq]
        countbad += len(badseqs)
        badcntr = Counter(s.seq for s in badseqs)
        for badseq in badcntr:
            fix_table[badseq] = goodseq
            badmut = MM.get_mutation(badseq)
            summarize_mismatches[(str(goodmut),str(badmut))] += badcntr[badseq]
            gdif,bdif = mutant_diffs(goodmut,badmut)
            summarize_diffs[(gdif,bdif)] += badcntr[badseq]

    print(f"Found {countbad} inconsistent sequences")
    print(f"Found {len(inconsistent_shortseqs)} cases of inconsistency")
    print()
    write_summary(summarize_mismatches,summarize_diffs)
    write_output(args.output,first,fix_table,shareshortseq,inconsistent_shortseqs)

    ## Note: number of cases in output, and number of differences in summary
    ## might not agree, since multiple badseq's for a single goodseq
    ## distinct cases = distinct goodseq's
    ## distinct diffs = distinct goodseq,badseq tuples >= distinct cases
    
    if args.fix:
        for s in seqs:
            s.seq = fix_table.get(s.seq,s.seq)
        readseq.write_seqfile(args.fix,seqs)

if __name__ == "__main__":

    _args = getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(_args)
