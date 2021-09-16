'''Try re-aligning sequences in a fasta file using skbio package'''
import sys
import re
from collections import Counter,defaultdict
from functools import lru_cache
import itertools as it
import argparse
from tqdm import tqdm

import skbio
#from blosum50 import blosum50

from seqsample import SequenceSample
import readseq
import sequtil
import covid
import mutant
import nw_align as nw

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(argparser)
    paa = argparser.add_argument
    paa("-N",type=int,default=0,
        help="only read in the first N sequences")
    paa("--consensus","-c",action="store_true",
        help="align to consensus instead of to first sequence")
    paa("--dedash",action="store_true",
        help="remove dashes from sequences")
    paa("--output","-o",
        help="output file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def dedash(seq):
    return re.sub("-","",seq)

def dedashiter(seqs):
    for s in seqs:
        s.seq = dedash(s.seq)
        yield s

alphabet='ABCDEFGHIJKLMNOPQRSTUVWXYZ-*#'
idmat = skbio.alignment.make_identity_substitution_matrix(2,-3,
                                                          alphabet=alphabet)

@lru_cache(maxsize=None)
def nwalign(fseq,cseq):
    '''possibly slow alignment scheme -- used for short substrings'''
    xfseq,xcseq,_ = nw.needleman_wunsch(fseq,cseq)
    return xfseq,xcseq
    
def align_seqs(fseq,cseq):
    '''main alignemnt: uses skbio's fast ssw algorithm, plus some fixes'''

    alignment,_,se = \
        skbio.alignment.local_pairwise_align_ssw(skbio.Protein(fseq),
                                                 skbio.Protein(cseq),
                                                 substitution_matrix=idmat)
                                                 #substitution_matrix=blosum50)
    xfseq = str(alignment[0])
    xcseq = str(alignment[1])
    [(fs,fe),(cs,ce)] = se

    if 0:

        ## local ssw does not align head or tail of the sequences
        ## so need to add those back ... with an appropriate number
        ## of -'s that they all line up  
    
        xfseq = fseq[:fs] + "-"*(cs-fs) + xfseq + fseq[1+fe:]
        xcseq = "-"*(fs-cs) + cseq[:cs] + xcseq + cseq[1+ce:]

    else:
        ## IDEA
        ## use NW algorithm to re-align the "head" and "tail"
        ## ok that NW is slow because these are short fragments
        ## and since short, NW can probably be cached
        fpre,cpre = '',''
        if fs or cs:
            fpre,cpre = nwalign(fseq[:fs],cseq[:cs])
        fpost,cpost = '',''
        if 1+fe < len(fseq) or 1+ce < len(cseq):
            fpost,cpost = nwalign(fseq[1+fe:],cseq[1+ce:])

        xfseq = fpre + xfseq + fpost
        xcseq = cpre + xcseq + cpost
        
        if len(xfseq) != len(xcseq):
            ## pad with dashes
            xfseq = xfseq + "-"*(len(xcseq)-len(xfseq))
            xcseq = xcseq + "-"*(len(xfseq)-len(xcseq))
        if len(xfseq) != len(xcseq):
            print("f:",len(xfseq),xfseq[:3],xfseq[-3:])
            print("c:",len(xcseq),xcseq[:3],xcseq[-3:])
            assert 0
    return xfseq,xcseq

@lru_cache(maxsize=None)
def get_manager(refseq):
    return mutant.MutationManager(refseq)

@lru_cache(maxsize=None)
def seq_to_mut(fseq,cseq):
    '''align cseq to fseq, then express cseq as an m-string, relative to aligned xfseq'''
    xfseq,xcseq = align_seqs(fseq,cseq)
    mut = get_manager(xfseq).get_mutation(xcseq)
    if len(dedash(xfseq)) > len(fseq):
        print(xfseq)
        print(dedash(xfseq))
        print(fseq)
        assert 0
    #if len(mut):
    #    print("max site:",max(ssm.site for ssm in mut),"vs",len(dedash(xfseq)),len(fseq))
    return mut

def update_xtrachars(mut,xtras):
    '''xtras dict keeps track of how many extra spaces we need at each site'''
    for ssm in mut:
        if ssm.ref == '+':
            xtras[ssm.site] = max(xtras.get(ssm.site,0), len(ssm.mut))
            
def newref(ref,xtras):
    '''use xtras dict to pad the ref sequence with dashes'''
    nref = ""
    for site,rchar in enumerate(ref,start=1):
        nref += rchar
        if xtras.get(site,None):
            nref += "-"*xtras[site]
    return nref
        
def _main(args):
    '''realign main'''

    seqs = covid.read_filter_seqfile(args)
    if args.dedash:
        seqs = dedashiter(seqs)
    if args.N:
        seqs = it.islice(seqs,args.N+1)
        
    first, seqs = sequtil.get_first_item(seqs)
    refseq = first.seq

    outputseqs = []
    def ssappend(name,seq):
        if args.output:
            outputseqs.append(SequenceSample(name,seq))


    if args.consensus:
        seqs = list(seqs)
        vprint("getting consensus...")
        refseq = sequtil.consensus(seqs[1:])
        seqs[0] = SequenceSample("cons",refseq)
        vprint("Consensus defined...")

    xtras = dict()
            
    for s in seqs: #tqdm(seqs):
        mut = seq_to_mut(refseq,s.seq)
        ssappend(s.name,mut) ## abusing SequenceSample structure
        update_xtrachars(mut,xtras)

    vprint("xtras:",xtras)
    refseq = newref(refseq,xtras)
    vprint("refseq:",refseq)

    m_mgr = mutant.MutationManager(refseq)
    for s in outputseqs:
        oldmut = s.seq
        try:
            ## convert 'abused' s.seq as mut into an actual sequence
            s.seq = m_mgr.seq_from_mutation(s.seq)
        except (IndexError,TypeError) as e:
            ## possible error: if mut includes [+0ABC]
            print("name:",s.name)
            print("oldmut:",oldmut)
            raise RuntimeError(e)
        ## this check appears to be cheap (20s in an 8m run)
        ## and makes sure things are consistent
        newmut = m_mgr.get_mutation(s.seq)
        if newmut != oldmut:
            raise RuntimeError(f"Error {s.name}\n {oldmut} {newmut}")

    if args.output:
        readseq.write_seqfile(args.output,outputseqs)

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
