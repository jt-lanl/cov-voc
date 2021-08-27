'''Tweak alignment according to specifed mutant strings'''
import sys
import re
import itertools as it
import argparse

import covid
import readseq
import mutanty as mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--mstrings","-m",nargs=2,
        help="pair of bad/good mstrings")
    paa("--mfile","-M",
        help="read bad/good mstrings from a file")
    paa("--output","-o",
        help="output tweaked fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def de_gap(seq):
    '''remove '-'s from sequence'''
    return re.sub('-','',seq)

def mstrings_to_ndx_seqs(MM,mstring_a,mstring_b):

    mut_a = mutant.Mutation.from_mstring(mstring_a)
    mut_b = mutant.Mutation.from_mstring(mstring_b)

    sites = sorted(set(ssm.site for ssm in it.chain(mut_a,mut_b)))
    lo,hi = sites[0],sites[-1]+1
    ndxlo = min(MM.indices_from_site(lo))
    ndxhi = max(MM.indices_from_site(hi-1))+1

    seq_r = MM.refseq
    seq_a = MM.seq_from_mutation(mut_a)
    seq_b = MM.seq_from_mutation(mut_b)

    seq_r = seq_r[ndxlo:ndxhi]
    seq_a = seq_a[ndxlo:ndxhi]
    seq_b = seq_b[ndxlo:ndxhi]

    ## None of these should happen
    if seq_r == seq_a:
        raise RuntimeError(f"Edit {mstring_a} will be inconsistent with ref sequence!")
    if seq_a == seq_b:
        vprint(f"Edit {mstring_a}->{mstring_b} will do nothing!")
    if de_gap(seq_a) != de_gap(seq_b):
        raise RuntimeError(f"Edit {mstring_a}->{mstring_b} will change actual sequence!"
                           " not just the alignment")

    return ndxlo,ndxhi,seq_a,seq_b

def read_mstring_pairs(filename):
    mstringpairs=[]
    with open(filename) as f_in:
        for line in f_in:
            line = line.strip()
            line = re.sub('#.*','',line)
            m = re.match(r'.*(\[.*\]).+(\[.*\]).*',line)
            if not m:
                if line:
                    vprint("Could not read line:",line)
                continue
            mstringpairs.append( (m[1],m[2]) )
    return mstringpairs

def _main(args):
    '''tweakalign main'''

    ## Read full sequences
    args.keepdashcols = True
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,putitemback=False)

    MM = mutant.MutationManager(first.seq)

    ## Get bad and good mutations
    ndx_seqs_tuples = []
    if args.mfile:
        mstringpairs = read_mstring_pairs(args.mfile)
        for ma,mb in mstringpairs:
           ndxlo,ndxhi,seq_a,seq_b = mstrings_to_ndx_seqs(MM,ma,mb)
           ndx_seqs_tuples.append( (ndxlo,ndxhi,seq_a,seq_b) )
    if args.mstrings:
        assert len(args.mstrings)==2
        ma,mb = args.mstrings
        ndxlo,ndxhi,seq_a,seq_b = mstrings_to_ndx_seqs(MM,ma,mb)
        ndx_seqs_tuples.append( (ndxlo,ndxhi,seq_a,seq_b) )

    changed_sequences=[]
    seqs = list(seqs)
    for s in seqs:
        for ndxlo,ndxhi,seq_a,seq_b in ndx_seqs_tuples:
            if s.seq[ndxlo:ndxhi] == seq_a:
                s.seq = s.seq[:ndxlo] + seq_b + s.seq[ndxhi:]
                changed_sequences.append(s)

    print("Made",len(changed_sequences),"changes, affecting",len(set(changed_sequences)),"distinct sequences")
    if args.output:
        seqs = [first] + seqs
        readseq.write_seqfile(args.output,seqs)

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
