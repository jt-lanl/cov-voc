'''Tweak alignment according to specifed mutant strings'''
import re
import itertools as it
from collections import Counter,defaultdict
import random
import argparse

import warnings

import verbose as v
import covid
import readseq
import sequtil
import mutant
import mstringfix

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--mstrings","-m",nargs=2,action="append",
        help="pair of from/to mstrings")
    paa("--mfile","-M",action="append",
        help="read from/to mstring pairs from a file")
    paa("--checktweaks",action="store_true",
        help="Check for dependencies among the tweaks [EXPERIMENTAL!]")
    paa("--shuffletweaks",action="store_true",
        help="apply the mstring-pair tweaks in random order")
    paa("--output","-o",
        help="output tweaked fasta file")
    paa("--jobno",type=int,default=1,
        help="job number if using tweakalign in parallel")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def de_gap(seq):
    '''remove '-'s from sequence'''
    return re.sub('-','',seq)

def mstrings_to_ndx_seqs(mut_mgr,mstring_a,mstring_b):
    '''
    convert mstrings into short sequence-alignment fragments
    along with the indices of where those sequences start/end
    '''
    mut_a = mutant.Mutation.from_mstring(mstring_a)
    mut_b = mutant.Mutation.from_mstring(mstring_b)

    sites = sorted(set(ssm.site for ssm in it.chain(mut_a,mut_b)))
    lo,hi = sites[0],sites[-1]+1
    ndxlo = min(mut_mgr.indices_from_site(lo))
    ndxhi = max(mut_mgr.indices_from_site(hi-1))+1

    seq_r = mut_mgr.refseq
    seq_a = mut_mgr.seq_from_mutation(mut_a)
    seq_b = mut_mgr.seq_from_mutation(mut_b)

    seq_r = seq_r[ndxlo:ndxhi]
    seq_a = seq_a[ndxlo:ndxhi]
    seq_b = seq_b[ndxlo:ndxhi]

    v.vvprint(f'{mstring_a}->{mstring_b}:')
    v.vvprint(f'   r: {seq_r}')
    v.vvprint(f'   a: {seq_a}')
    v.vvprint(f'   b: {seq_b}')

    ## None of these should happen
    if seq_r == seq_a:
        raise RuntimeError(f"Edit {mstring_a} will be inconsistent "
                           "with ref sequence!")
    if seq_a == seq_b:
        v.vprint(f"Edit {mstring_a}->{mstring_b} will do nothing!")
    if de_gap(seq_a) != de_gap(seq_b):
        v.print(".".join(a+b for a,b in zip(de_gap(seq_a),de_gap(seq_b))))
        v.print(f'   r: {seq_r}')
        v.print(f'   a: {seq_a}')
        v.print(f'   b: {seq_b}')
        warnings.warn(f"Edit {mstring_a}->{mstring_b} "
                      "will change actual sequence!"
                      " not just the alignment\n"
                      f"  {seq_a}->{seq_b}")
        ## this shouldn't happen...but if it does, then don't do any replacing
        seq_b = seq_a

    return ndxlo,ndxhi,seq_a,seq_b

def read_mstring_pairs(filename):
    '''
    wrapper of the mstringfix.read_mstring_pairs
    that returns a list instead of an iterator
    '''
    return list(mstringfix.read_mstring_pairs(filename))

def add_xtra_dashes(seqs,xxtras):
    '''xxtras is dict keyed by indices of s.seq strings;
       for each of those strings, we expand by extra dashes
    '''
    for s in seqs:
        sseq = list(s.seq)
        for ndx in xxtras:
            sseq[ndx] += "-"*xxtras[ndx]
        s.seq = "".join(sseq)
        yield s


def gen_seq_clumps(allseqs):
    '''yields lists of sequences (aka clumps);
    a clump is a list that begins with a reference sequence
    and has no further reference sequences in it'''
    first,seqs = sequtil.get_first_item(allseqs,keepfirst=False)
    clump = [first]
    for s in seqs:
        if s.name == first.name:
            yield clump
            clump = [s]
        else:
            clump.append(s)
    yield clump


def get_mstring_pairs(args):
    '''
    Return a list of 2-tuples of m-strings;
    Each pair has a from_mstring and a to_mstring
    '''
    mstringpairs = []

    ## first: read the file(s) indicated by '-M file'
    ## A file named '.' indicates that one should use the defaults
    ## (which are currently listed in mstringfix.py)
    ## Note that '-M' can be invoked more than once
    ##   That is: '-M file1 -M file2' is okay
    ##   Invalid: '-M file1 file2'
    for mfile in args.mfile or []:
        if mfile == '.':
            mfile = None
        mstringpairs.extend( read_mstring_pairs(mfile) )

    ## second: read from the strings on the command line
    ## these are of the from -m from_mstring to_mstring
    ## and multiple invocations of '-m' are allowed
    for mspair in args.mstrings or []:
        assert len(mspair)==2
        mstringpairs.append(mspair)

    ## finally: if after all that, there are no mstring pairs
    ## in the mstringpairs list, then use the defaults
    if not mstringpairs:
        mstringpairs.extend( read_mstring_pairs(None) )

    ## ensure mstrings have brackets around them
    mstringpairs = [ (mstringfix.mstring_brackets(a),
                      mstringfix.mstring_brackets(b))
                     for a,b in mstringpairs ]

    ## are there any duplicates?
    if len(mstringpairs) != len(set(mstringpairs)):
        mspairs_sofar = set()
        for mspair in mstringpairs:
            if mspair in mspairs_sofar:
                v.print('Duplicated pair:',mspair[0],mspair[1])
            mspairs_sofar.add(mspair)
        raise RuntimeError('Duplicated pair(s) of mstrings specified')

    ## Now we have all of our mstrings
    v.vprint("mstring pairs:",len(mstringpairs))
    for mspair in mstringpairs:
        v.vvprint("%s -> %s" % mspair)

    return mstringpairs

def check_tweaks(mstringpairs):
    '''Given a list of mstringpair tweaks, check to see if there
    are dependencies that could lead to different alignments
    if the tweaks were applied in different order'''

    ## From list of tweaks
    ## Build a dictionary of lists, keyed on site number
    ## List is of changes at that site
    ##   eg, A19C->A19D is C19D or CD
    ## As new items are added to the list
    ##   1/ check if on the list
    ##   2/ see if it over-rides something on the list
    ##     eg, CD,BE are ok, but
    ##     a/ CD,CE would make result depend on order
    ##     b/ CD,DE would override (also depend on order)
    ##     c/ CD,BD would be okay I think


    sitetweaks = defaultdict(set)
    warnlist = list()
    for ma,mb in mstringpairs:
        ssms_a = mutant.Mutation(ma)
        ssms_b = mutant.Mutation(mb)
        dtmp_a = {(ssm.ref,ssm.site): ssm.mut
                  for ssm in ssms_a}
        dtmp_b = {(ssm.ref,ssm.site): ssm.mut
                  for ssm in ssms_b}
        keys = sorted(set(dtmp_a)|set(dtmp_b))
        v.vprint(f'{ma}->{mb}: keys={keys}')
        for k in keys:
            amut,bmut = (dtmp_a.get(k,''),dtmp_b.get(k,''))
            if (amut,bmut) in sitetweaks[k]:
                continue
            for (a,b) in sitetweaks[k]:
                if b==bmut:
                    continue
                if (a==amut and b!=bmut) or (a == bmut) or (b == amut):
                    warnlist.append(f'{k[0]}{k[1]}: {amut}{bmut} vs {a}{b}')
            sitetweaks[k].add((amut,bmut))
    if warnlist:
        v.print("Warnings:")
        for warn in warnlist:
            v.print(warn)
    return warnlist

def get_extra_chars(mstringpairs):
    '''
    Characterize extra chars:
    xtras[site] = number of extra chars after site
    Add xtra chars for all the +nnnABC mstring component
    '''
    xtras = dict()
    for mspair in mstringpairs:
        for mstring in mspair:
            for ssm in mutant.Mutation.from_mstring(mstring):
                if ssm.ref == "+":
                    xtras[ssm.site] = max( [xtras.get(ssm.site,0),
                                            len(ssm.mut)] )
    return xtras

def _main(args):
    '''tweakalign main'''
    v.vprint(args)

    ## Get all the mstring pair
    mstringpairs = get_mstring_pairs(args)

    if args.checktweaks:
        check_tweaks(mstringpairs)
        return

    if args.shuffletweaks:
        mstringpairs = random.sample(mstringpairs,k=len(mstringpairs))
    ## Identify where extra columns are needed
    xtras = get_extra_chars(mstringpairs)
    v.vprint("xtras:",xtras)

    ## Read full sequences
    allseqs = covid.read_filter_seqfile(args)
    veryfirst = None
    outseqs = []
    for seqs in gen_seq_clumps(allseqs):
        first,seqs = sequtil.get_first_item(seqs)
        if veryfirst is None:
            veryfirst = first
        if args.jobno == 1 or veryfirst != first:
            outseqs.append(first)
        v.vprint('chunk: len=',len(seqs),'first=',first.name)

        mut_mgr = mutant.MutationManager(first.seq)
        v.vprint("len(first):",len(first.seq))

        xxtras = dict() ## xx is ndx based instead of site based
        for site in sorted(xtras):
            mm_xtras = len(mut_mgr.indices_from_site(site))-1
            if mm_xtras < xtras[site]:
                v.print("xtras:",site,mm_xtras,"->",xtras[site])
                ndx = mut_mgr.index_from_site(site)
                xxtras[ndx] = xtras[site] - mm_xtras

        if len(xxtras):
            seqs = add_xtra_dashes(seqs,xxtras)
        first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
        mut_mgr = mutant.MutationManager(first.seq)
        v.vprint("len(first):",len(first.seq))

        ndx_seqs_tuples = []
        for ma,mb in mstringpairs:
            ndxlo,ndxhi,seq_a,seq_b = mstrings_to_ndx_seqs(mut_mgr,ma,mb)
            ndx_seqs_tuples.append( (ma,mb,ndxlo,ndxhi,seq_a,seq_b) )
            v.vvprint(f"{ma} -> {mb}:   {seq_a} -> {seq_b}")

        changed_sequences=[]
        seqs = list(seqs)
        changes_by_mstring = Counter()
        for s in seqs:
            ## Normalize sequence (seq -> mut -> seq)
            nmut = mut_mgr.get_mutation(s.seq)
            s.seq = mut_mgr.seq_from_mutation(nmut)
            for ma,mb,ndxlo,ndxhi,seq_a,seq_b in ndx_seqs_tuples:
                if s.seq[ndxlo:ndxhi] == seq_a:
                    s.seq = s.seq[:ndxlo] + seq_b + s.seq[ndxhi:]
                    changed_sequences.append(s)
                    changes_by_mstring[(ma,mb)] += 1

        v.vprint("Made",len(changed_sequences),"changes, affecting",
                len(set(changed_sequences)),"distinct sequences")

        v.vprint("Changes by mstring:")
        for (ma,mb),cnt in changes_by_mstring.items():
            v.vprint(f'{cnt:8d} {ma}->{mb}')

        outseqs.extend(seqs)
        v.vprint('outseqs:',len(outseqs))

    if args.output:
        readseq.write_seqfile(args.output,outseqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
