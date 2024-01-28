'''fix alignemnt by enforcing consistency over various subregions'''

import re
from functools import cache
from collections import defaultdict,Counter
import argparse

import verbose as v
import readseq
import sequtil
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input sequence file")
    paa("--badseqs","-b",
        help="output inconsistent sequences")
    paa("--mstringpairs",
        help="write file of mstring pairs")
    paa("--fix",
        help="fix sequences and write output to this fasta file")
    paa("--sitepartition","-s",nargs='+',
        help="partition full sequence into subsequences")
    paa("--windowsize","-w",type=int,default=0,
        help="subsequence window size")
    paa("--nsweeps",type=int,default=2,
        help="number of sweeps through the sequences, to get overlaps")
    paa("--phases",type=int,nargs='+',
        help="Which sweeps to do in this run: subset of 0,...,nsweeps-1")
    paa("--na",action="store_true",
        help="set for nucleotide alignment (default is amino acid alignment)")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args


dedash = re.compile("-")

@cache
def de_gap(seq):
    '''remove '-'s from sequence'''
    return dedash.sub("",seq)

def choose_alignment(mut_mgr,gseqs,nuc_align=False):
    '''
    given a batch of gseqs whose de-gapped dseqs are identical;
    choose the gseqs with
    1/ smallest mutation length
    2/ most common form
    '''

    ## take care of special case right away
    if len(set(gseqs))==1:
        return gseqs[0]

    ## first, make sure that de-gapped sequences are identical
    assert len(set(de_gap(seq) for seq in gseqs)) == 1
    ## find forms with smallest mutation length


    if not nuc_align:
        ## makes sense for amino acids, maybe not so much for nucleotides
        mutlen = {gseq: len(mut_mgr.get_mutation(gseq)) for gseq in gseqs}
        #mutlen = {gseq: mut_mgr.get_hamming(gseq) for gseq in gseqs}
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
        ssm.site += site_offset
        #ssm.adjust_site(site_offset)
    return str(mut)

def deletion_only_solution(refseq,dseq):
    '''
    return a string that matches refseq using dseq characters (in order) and dashes
    return None if not all of dseq characters are used
    '''
    ## "A--BC--D--E--FG-H","BCDEG" --> '---BC--D--E------'
    ## "A--BC--D--E--FG-H","BDCEG" --> None
    dseq = iter(dseq)
    nextd = next(dseq)
    dos = ""
    for r in refseq:
        if r == nextd:
            dos += nextd
            nextd = next(dseq,None)
        else:
            dos += "-"
    return None if nextd else dos

def check_subsequences(subseqset):
    '''returns True if aligned and otherwise False'''
    ## strategy is to ensure that if two "gapped" subseqs (gseqs)
    ## have the same de-gapped string (dseq), then they are identical

    v.vvvprint('distinct subsequences:',len(subseqset))

    gseq_with_dseq = dict()
    retval=True
    for gseq in subseqset:
        dseq = de_gap(gseq)
        v.vvvprint_only(5,'gseq/dseq:',f'{gseq=} {dseq=}')
        if dseq in gseq_with_dseq:
            if gseq != gseq_with_dseq[dseq]:
                v.print(f'{dseq=}: {gseq} != {gseq_with_dseq[dseq]}')
                retval=False
                return False ## comment out if you want to see all inconsistencies
        else:
            gseq_with_dseq[dseq]=gseq
    #return True
    return retval

def align_subsequences(subseqs,site_offset=0,nuc_align=False,badgoodmuts=None):
    '''return a list of subsequences that are aligned'''

    firstseq,subseqs = subseqs[0],subseqs[1:]
    mut_mgr = mutant.MutationManager(firstseq)

    ## two forms of every sequence: gseq (gapped), dseq (de-gapped)

    ## gseqs_with_dseq[dseq] is a list of gseqs whose de-gapped seqstring is dseq
    gseqs_with_dseq = defaultdict(list)
    for gseq in subseqs:
        gseqs_with_dseq[de_gap(gseq)].append(gseq)

    fix_table = dict()
    dfirstseq = de_gap(firstseq)
    for dseq,gseqs in gseqs_with_dseq.items():
        if nuc_align and len(set(gseqs))==1:
            ## skip fixing of consistently bad alignments
            ## but don't trust dos for nucleotides
            ## (for nuc_align, might prefer something that preserved triplets)
            continue

        goodseq = firstseq if dseq==dfirstseq \
            else choose_alignment(mut_mgr,gseqs,nuc_align=nuc_align)
        goodmut = mut_mgr.get_mutation(goodseq)
        good_mstring = str(siteadjust(goodmut,site_offset))

        ## check for a special case; doesn't happen very often
        if not nuc_align and any(ssm.mut != '-' for ssm in goodmut):
            ## if not already a deletion-only solution,
            ## then see if one exists
            dos = deletion_only_solution(firstseq,dseq)
            if dos:
                ## if so, then use it instead
                goodseq = dos
                goodmut = mut_mgr.get_mutation(goodseq)
                good_mstring = str(siteadjust(goodmut,site_offset))

        badseqs = [seq for seq in gseqs if seq != goodseq]
        if not badseqs:
            continue
        badseqs_counter = Counter(seq for seq in badseqs)
        v.vprint(f"{firstseq} ref")
        v.vprint(f"{goodseq} good {siteadjust(goodmut,site_offset)} "
               f"count={len(gseqs)-len(badseqs)}")
        for badseq in badseqs_counter:
            badmut = mut_mgr.get_mutation(badseq)
            bad_mstring = str(siteadjust(badmut,
                                         site_offset))
            v.vprint(f"{badseq} bad {bad_mstring} "
                     f"count={badseqs_counter[badseq]}")
            v.vprint(f'ffix-mstring: {bad_mstring} {good_mstring}') ## full fix
            ## remove common ssm's
            commonmut = set(badmut) & set(goodmut)
            xbadmut = set(badmut) - commonmut
            xgoodmut = set(goodmut) - commonmut
            xbad_mstring = str(siteadjust(mutant.Mutation(xbadmut),site_offset))
            xgood_mstring = str(siteadjust(mutant.Mutation(xgoodmut),site_offset))
            v.vprint(f'sfix-mstring: {bad_mstring} {good_mstring}') ## short fix (too short?)
            ## add back common ssm's that are in range
            minsite = min(min(mut.site for mut in xbadmut),
                          min(mut.site for mut in xgoodmut))
            maxsite = max(max(mut.site for mut in xbadmut),
                          max(mut.site for mut in xgoodmut))
            for mut in commonmut:
                if minsite < mut.site < maxsite:
                    xbadmut.add(mut)
                    xgoodmut.add(mut)
            xbadmut = mutant.Mutation(xbadmut)
            xgoodmut = mutant.Mutation(xgoodmut)
            xbad_mstring = str(siteadjust(xbadmut,site_offset))
            xgood_mstring = str(siteadjust(xgoodmut,site_offset))
            v.vprint(f'xfix-mstring: {xbad_mstring} {xgood_mstring}') ## corrected fix
            ## add back explicit identity ssm's
            for site in range(minsite,maxsite+1):
                ref = mut_mgr.refval(site)
                if site not in [ssm.site for ssm in xbadmut]:
                    xbadmut.append(mutant.SingleSiteMutation(f'{ref}{site}{ref}'))
                if site not in [ssm.site for ssm in xgoodmut]:
                    xgoodmut.append(mutant.SingleSiteMutation(f'{ref}{site}{ref}'))
            xbad_mstring = str(siteadjust(xbadmut,site_offset))
            xgood_mstring = str(siteadjust(xgoodmut,site_offset))
            v.vprint(f'zfix-mstring: {xbad_mstring} {xgood_mstring}') ## corrected fix
            fix_table[badseq] = goodseq
            if badgoodmuts is not None:
                badmut = mut_mgr.get_mutation(badseq)
                badgoodmuts.append( (str(siteadjust(badmut,site_offset)),
                                     str(siteadjust(goodmut,site_offset))) )

    return [firstseq] + [fix_table.get(gseq,gseq) for gseq in subseqs]

def mk_subseq_site_ranges(top,window,nsweeps=2,phases=None):
    ''' return a list of site lo,hi pairs with lo starting at 1,
and hi ending at top, with hi-lo=window (or smaller than window
on the edges), and with nsweeps separate overlapping sweeps'''
    assert nsweeps <= window
    pairs = []
    phases = phases or range(nsweeps)
    offsetlist = [1-n*window//nsweeps for n in phases]
    for offset in offsetlist:
        r = range(offset,top+window,window)
        for lo,hi in zip(r[:-1],r[1:]):
            pairs.append( (max([1,lo]),min([top,hi])) )
    return pairs

def ndx_ranges(xlator,subseq_ranges):
    '''convert ranges of site numbers to ranges of indices;
    with conversion done by SiteIndexTranslator 'xlator'
    '''
    for lo,hi in subseq_ranges:

        if lo == 'x' or hi == 'x':
            continue

        lo = 1                if lo=="." else int(lo)
        hi = xlator.topsite+1 if hi=="." else int(hi)

        lo = max(lo,1)

        if lo > hi:
            v.vprint("Out of order site range:",lo,hi)
            continue

        if lo >= xlator.topsite:
            # we're done here
            break

        if hi > xlator.topsite:
            hi = xlator.topsite+1

        ndxlo = min(xlator.indices_from_site(lo))
        ndxhi = max(xlator.indices_from_site(hi-1))+1

        v.vvvprint(f"sites {lo}:{hi}, indices {ndxlo}:{ndxhi}")

        yield lo,hi,ndxlo,ndxhi

def _main(args):
    '''fixalign main'''
    v.vprint(args)

    seqs = sequtil.read_seqfile(args.input)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)

    xlator = mutant.SiteIndexTranslator(first.seq)

    changed_sequences=[]
    bad_good_mstrings=[]

    ## partition seq into overlapping windows of width 2*stepsize
    winsize = 60 if args.na else 20
    if args.windowsize:
        winsize = args.windowsize
        if args.na:
            winsize = 3*(winsize//3)

    subseq_ranges = mk_subseq_site_ranges(xlator.topsite+1,winsize,
                                          nsweeps=args.nsweeps,
                                          phases=args.phases)

    if args.sitepartition:
        ## over-rides windowsize/nsweep
        subseq_ranges = list(zip(args.sitepartition[:-2],
                                 args.sitepartition[2:]))

    v.vvvprint('ranges:',subseq_ranges)

    ## save til now to make generator a list
    ## that way, can debug subseq_ranges
    seqs=list(seqs)

    for lo,hi,ndxlo,ndxhi in ndx_ranges(xlator,subseq_ranges):

        v.vvprint(f'Checking sites {lo}:{hi} / indices {ndxlo}:{ndxhi}')
        subseqset = set(s.seq[ndxlo:ndxhi] for s in seqs)
        check = check_subsequences(subseqset)
        if check:
            continue

        v.vprint(f'Fixing sites {lo}:{hi} / indices {ndxlo}:{ndxhi}')
        subseqs = [s.seq[ndxlo:ndxhi] for s in seqs]
        subseqs = align_subsequences(subseqs,
                                     site_offset=lo-1,
                                     nuc_align=args.na,
                                     badgoodmuts=bad_good_mstrings)

        assert len(seqs) == len(subseqs)
        countbad=0
        for s,subseq in zip(seqs,subseqs):
            if s.seq[ndxlo:ndxhi] != subseq:
                countbad += 1
                assert len(s.seq[ndxlo:ndxhi]) == len(subseq)
                v.vvprint_only(1,f'{s.seq[ndxlo:ndxhi]} -> {subseq}:',s.name)
                s.seq = s.seq[:ndxlo] + subseq + s.seq[ndxhi:]
                changed_sequences.append(s)

        if countbad:
            v.vprint("Changed",countbad,"inconsistent sequences in range:",lo,hi-1)

    v.vprint("Total:",
             len(changed_sequences),"changes in",
             len(set(changed_sequences)),"distict sequences",
             f"{args.windowsize=} {args.phases=}")
    print("Total:",
          len(changed_sequences),"changes in",
          len(set(changed_sequences)),"distict sequences",
          f"{args.windowsize=} {args.phases=}")

    if args.mstringpairs:
        with open(args.mstringpairs,"w") as fileptr:
            for bad,good in bad_good_mstrings:
                fileptr.write(f"{bad} {good}\n")

    if args.badseqs:
        changed_set = set(changed_sequences)
        oseqs = [first] + [s for s in seqs if s in changed_set]
        readseq.write_seqfile(args.badseqs,oseqs)

    if args.fix:
        readseq.write_seqfile(args.fix,seqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
