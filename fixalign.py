'''fix alignemnt by enforcing consistency over verious subregions'''
import sys
import re
from functools import lru_cache
from collections import defaultdict,Counter
import argparse

import readseq
import mutantx as mutant
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--output","-o",
        help="output inconsistent sequences")
    paa("--fix",
        help="fix sequences and write output to this fasta file")
    paa("--sitepartition","-s",nargs='+',
        help="partition full sequence into subsequences")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args


dedash = re.compile("-")

@lru_cache(maxsize=None)
def de_gap(seq):
    '''remove '-'s from sequence'''
    return dedash.sub("",seq)

def choose_alignment(MM,gseqs):
    '''
    given a batch of gseqs whose de-gapped dseqs are identical;
    choose the gseqs with
    1/ smallest mutation length
    2/ most common form
    '''
    ## first, make sure that de-gapped sequences are identical
    dseqs = [de_gap(seq) for seq in gseqs]
    assert all(dseq == dseqs[0] for dseq in dseqs[1:])
    ## find forms with smallest mutation length
    if 1:
        ## makes sense for amino acids, maybe not so much for nucleotides
        mutlen = {gseq: len(MM.get_mutation(gseq)) for gseq in gseqs}
        #mutlen = {gseq: MM.get_hamming(gseq) for gseq in gseqs}
        minmutlen = min(mutlen.values())
        candidateseqs = [gseq for gseq in gseqs if mutlen[gseq] == minmutlen]
    else:
        candidateseqs = gseqs
    ## find most common form
    cnt = Counter(candidateseqs)
    [(cseq,_)] = cnt.most_common(1)
    return cseq

def siteadjust(mstring,site_offset=0):
    '''adjust site numbers in ssm's in mutstring'''
    if not site_offset:
        return mstring
    mut = mutant.Mutation(str(mstring))
    for ssm in mut:
        ssm.adjust_site(site_offset)
    return str(mut)

def align_subsequences(subseqs,site_offset=0):
    '''return a list of subsequences that are aligned'''

    first,subseqs = subseqs[0],subseqs[1:]
    MM = mutant.MutationMaker(first)

    ## two forms of every sequence: gseq (gapped), dseq (de-gapped)

    ## gseqs_with_dseq[dseq] is a list of gseqs whose de-gapped seqstring is dseq
    gseqs_with_dseq = defaultdict(list)
    for gseq in subseqs:
        gseqs_with_dseq[de_gap(gseq)].append(gseq)

    ## find the short dseqs with inconsistent long gseqs
    inconsistent_dseqs = set(dseq
                             for dseq,gseqs in gseqs_with_dseq.items()
                             if any(seq != gseqs[0] for seq in gseqs[1:]))

    countbad = 0
    fix_table = dict()
    for dseq in inconsistent_dseqs:   ## check if dseq == de_gap(first), then goodmut = first
        gseqs = gseqs_with_dseq[dseq]
        goodseq = choose_alignment(MM,gseqs)
        goodmut = MM.get_mutation(goodseq)
        badseqs = [seq for seq in gseqs if seq != goodseq]
        countbad += len(badseqs)
        badcntr = Counter(seq for seq in badseqs)
        vprint(f"{first} ref")
        vprint(f"{goodseq} good {siteadjust(goodmut,site_offset)} "
               f"count={len(gseqs)-len(badseqs)}")
        for badseq in badcntr:
            vprint(f"{badseq} bad  {siteadjust(MM.get_mutation(badseq),site_offset)} "
                   f"count={badcntr[badseq]}")
            fix_table[badseq] = goodseq

    if countbad:
        vprint("Inconsistent sequences:",countbad,
               "at range start site:",site_offset+1)

    return [first] + [fix_table.get(gseq,gseq) for gseq in subseqs]



def _main(args):
    '''fixalignment main'''

    ## Read full sequences
    args.keepdashcols = True
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs)

    changed_sequences=[]

    ## note, should check that sitepartition is always increasing, and check early

    T = mutant.SiteIndexTranslator(first.seq)
    for lo,hi in zip(args.sitepartition[:-2],
                     args.sitepartition[2:]):

        if lo == 'x' or hi == 'x':
            continue

        site_lo = lo = 1                           if lo=="." else int(lo)
        hi = T.site_from_index(len(first.seq)-1)-1 if hi=="." else int(hi)

        if lo > hi:
            vprint("Out of order site range:",lo,hi)
            continue

        if lo == hi:
            vprint("Empty site range:",lo,hi)
            continue

        if lo >= T.topsite():
            vprint(f"lo = {lo} is larger than size of sequence")
            break

        if hi >= T.topsite():
            hi = T.topsite()-1

        ndxlo = min(T.indices_from_site(lo))
        ndxhi = max(T.indices_from_site(hi))

        vvprint(f"sites {lo}:{hi}, indices {ndxlo}:{ndxhi}")

        seqs=list(seqs)

        subseqs = [s.seq[ndxlo:ndxhi] for s in seqs]
        subseqs = align_subsequences(subseqs,site_offset=site_lo-1)

        for s,subseq in zip(seqs,subseqs):
            if s.seq[ndxlo:ndxhi] != subseq:
                if len(s.seq[ndxlo:ndxhi]) != len(subseq):
                    raise RuntimeError(f"lengths disagree! "
                                       f"{len(s.seq[ndxlo:ndxhi])} != {len(subseq)}")
                s.seq = s.seq[:ndxlo] + subseq + s.seq[ndxhi:]
                changed_sequences.append(s)

    print("changed:",
          len(changed_sequences),"total changes in",
          len(set(changed_sequences)),"distict sequences")

    if args.output:
        changed_set = set(changed_sequences)
        oseqs = [first] + [s for s in seqs if s in changed_set]
        readseq.write_seqfile(args.output,oseqs)

    if args.fix:
        readseq.write_seqfile(args.fix,seqs)

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
