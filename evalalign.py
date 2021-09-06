'''Evaluate alignment by computing various summary statistics'''
import sys
import itertools as it
import argparse

import covid
import sequtil
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("-N",type=int,default=0,
        help="Only read N sequences")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def count_dashes(seq,dash='-'):
    '''count how many dashes, and how many runs of dashes in sequence'''
    ndxlist = sequtil.str_indexes(seq,dash)
    dcount = len(ndxlist)
    if dcount == 0:
        return 0,0
    rcount = 1
    for n,m in zip(ndxlist[:-1],ndxlist[1:]):
        if m > n+1:
            rcount += 1
    return dcount,rcount

def _main(args):
    '''main'''
    vprint(args)
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs)
    m_mgr = mutant.MutationManager(first.seq)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    dcnt,rcnt = count_dashes(first.seq)
    print("Ref:",dcnt,rcnt)
    stripseqs,seqs = it.tee(seqs)
    stripseqs = (s.copy() for s in stripseqs)
    stripseqs = sequtil.stripdashcols(first.seq,stripseqs)

    dcnt_sum = rcnt_sum = 0
    ssm_sum = ins_sum = inschr_sum = 0
    sub_sum = 0
    nseq = 0
    for s,ss in zip(seqs,stripseqs):
        nseq += 1
        d,r = count_dashes(ss.seq)
        dcnt_sum += d
        rcnt_sum += r
        mut = m_mgr.get_mutation(s.seq)
        ssm_sum += len(mut)
        ins_sum += sum(1 for ssm in mut if ssm.ref=='+')
        inschr_sum += sum(len(ssm.mut) for ssm in mut if ssm.ref=='+')
        sub_sum += sum(1 for ssm in mut if ssm.mut != '-')

    print("Nseq:",nseq)
    print("d,r:",dcnt_sum/nseq,rcnt_sum/nseq)
    print("ssm:",ssm_sum/nseq,ins_sum/nseq,inschr_sum/nseq,sub_sum/nseq)
    print("d,r:",dcnt_sum,rcnt_sum)
    print("ssm:",ssm_sum,ins_sum,inschr_sum,sub_sum)

    print(dcnt,rcnt,dcnt_sum,rcnt_sum,sub_sum,ins_sum,inschr_sum)
    

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
