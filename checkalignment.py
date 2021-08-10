'''Check for alignment consistency'''
import sys
import re
from collections import defaultdict,Counter
import argparse

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

def de_gap(seq):
    '''remove '-'s from sequence'''
    return re.sub("-","",seq)

def mutstring(refseq,s):
    '''convert sequence into mutation string'''
    return mutant.Mutation().init_from_sequences(refseq,s.seq)

def report_inconsistency(refseq,sref,s):
    '''if s is inconsistent with sref, then describe inconsistency'''
    print(f"inconsistency: {sref.name} vs {s.name}")
    mref = mutstring(refseq,sref)
    mseq = mutstring(refseq,s)
    if mref != mseq:
        print(f"      mstring: {mref} vs {mseq}")
        msetref = set(str(ssm) for ssm in mref)
        msetseq = set(str(ssm) for ssm in mseq)
        print("           ie: ",
              "<"+",".join(sorted(msetref-msetseq))+">","vs",
              "<"+",".join(sorted(msetseq-msetref))+">")

def choose_alignment(refseq,seqs):
    '''
    given a set of sequences whose de-gapped versions are identical;
    choose one (the best one?)
    '''
    ## first, just make sure that de-gapped sequences are identical
    sseqs = [de_gap(s.seq) for s in seqs]
    assert all(sseq == sseqs[0] for sseq in sseqs[1:])
    ## find forms nearest to refseq
    mutlen = {s: len(mutant.Mutation().init_from_sequences(refseq,s.seq))
              for s in seqs}
    minmutlen = min(mutlen.values())
    candidateseqs = [s for s in seqs if mutlen[s] == minmutlen]
    ## find most common form
    cnt = Counter(s.seq for s in candidateseqs)
    [(cseq,_)] = cnt.most_common(1)
    return cseq

def fix_alignment(seqs,cseq):
    '''alter sequence alignemnt for list of sequences to all be the common sequence cseq'''
    for s in seqs:
        s.seq = cseq

def main(args):
    '''checkalignment main'''

    args.keepdashcols = True
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs)
    if args.fix:
        seqs=list(seqs)
    longseq = dict() ## key=shortseq, value= first s for which degap(s.seq)=shortesq
    inconsistents = defaultdict(list)
    for s in seqs:
        shortseq = de_gap(s.seq)
        if shortseq not in longseq:
            longseq[shortseq] = s
        sref = longseq[shortseq]
        if s.seq != sref.seq:
            inconsistents[sref].append(s)
            report_inconsistency(first.seq,sref,s)

    if len(inconsistents):
        slist = [first]
        for n,sref in enumerate(inconsistents,start=1):
            slist.append(SequenceSample(f"Inconsistency {n:3d}","-"*len(sref.seq)))
            slist.append(sref)
            for s in inconsistents[sref]:
                slist.append(s)
        if args.output:
            readseq.write_seqfile(args.output,slist)
            print(f"Wrote {len(slist)} seqeunces to file {args.output}")
    print(sum(len(inconsistents[s]) for s in inconsistents),"inconsistencies found")

    if args.fix:
        fix_table=dict()
        for sref in inconsistents:
            goodseq = choose_alignment(first.seq,[sref] + inconsistents[sref])
            for s in [sref] + inconsistents[sref]:
                if s.seq != goodseq:
                    fix_table[s.seq] = goodseq
        for s in seqs:
            if s.seq in fix_table:
                s.seq = fix_table[s.seq]
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
