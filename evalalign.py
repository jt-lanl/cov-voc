'''Evaluate alignment by computing various summary statistics'''

import itertools as it
import argparse
import verbose as v

import covid
import sequtil
import mutant

HEADER_KEY='''
 D: number of deletions in ref sequence
 R: number of runs of deletions in ref sequence
SD: sum of deletions in all sequence
SR: sum of runs of deletions in all sequence
SS: sum of substitutions in all sequence
SI: sum of insertions in all sequence
Es: entropy over the sites that are not "-" in the ref sequence
Ea: entropy over all site
'''

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--keepsites",
        help="integer list of sites to include")
    paa("--jobno",type=int,default=1,
        help="Use --jobno={#} if using GNU parallel")
    paa("--key",action="store_true",
        help="Print out the key to the header abbreviations")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
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

def empty_names(seqs):
    '''return the SequenceSample's with empty names to save memory'''
    for s in seqs:
        s.name=""
        yield s

def _main(args):
    '''main'''
    v.vprint(args)
    seqs = covid.read_filter_seqfile(args)
    if args.keepsites:
        first,seqs = covid.get_first_item(seqs,keepfirst=True)
        m_mgr = mutant.MutationManager(first.seq)
        seqs = covid.keepsites(m_mgr,seqs,args.keepsites)
        seqs = empty_names(seqs)

    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)

    dcnt,rcnt = count_dashes(first.seq)
    v.vprint("Ref:",dcnt,rcnt)

    seqs = list(seqs)
    E = sequtil.chunked_entropy(seqs,keepx=False)

    sitelist = range(1,1+m_mgr.topsite)
    ndxsites = [m_mgr.index_from_site(site) for site in sitelist]
    esites = sum(E[n] for n in ndxsites)
    eall = sum(E)
    v.vprint("Entropy:",esites,eall)

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
        mut = m_mgr.seq_to_mutation(s.seq)
        ssm_sum += len(mut)
        ins_sum += sum(1 for ssm in mut if ssm.ref=='+')
        inschr_sum += sum(len(ssm.mut) for ssm in mut if ssm.ref=='+')
        sub_sum += sum(1 for ssm in mut if ssm.mut != '-')

    v.vprint("Nseq:",nseq)
    v.vprint("d,r:",dcnt_sum/nseq,rcnt_sum/nseq)
    v.vprint("ssm:",ssm_sum/nseq,ins_sum/nseq,inschr_sum/nseq,sub_sum/nseq)
    v.vprint("d,r:",dcnt_sum,rcnt_sum)
    v.vprint("ssm:",ssm_sum,ins_sum,inschr_sum,sub_sum)

    if args.jobno==1:
        if args.key:
            v.print(HEADER_KEY)
        else:
            v.vprint(HEADER_KEY)
        print("\t".join("D R SD SR SS SI SC Es Ea".split()))
        print("\t".join("%d %d %d %d %d %d %d %.6f %.6f".split()) %
              (dcnt,rcnt,dcnt_sum,rcnt_sum,sub_sum,
               ins_sum,inschr_sum,esites,eall))

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
