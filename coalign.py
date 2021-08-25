'''Template doc string'''
import sys
import itertools as it
from pathlib import Path
import argparse

import readseq
import wrapgen

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--aa",type=Path,
        help="amino acid sequence file")
    paa("--nu",type=Path,
        help="nucleic acid sequence file")
    paa("-N",type=int,default=0,
        help="for debugging purposes, just do the first N sequences")
    paa("--output","-o",
        help="write new nucleic acid file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

class seqdict():
    '''roughly dict-like, but it fills itself from iterator seqs only as needed'''
    def __init__(self,seqs):
        self.seqs = iter(seqs)
        self.D = dict()

    def get(self,key,default=None):
        while key not in self.D:
            s = next(self.seqs,None)
            if s is not None:
                self.D[s.name] = s.seq
            else:
                return default
        return self.D[key]

    ## NOTE: dunder-Methods below make seqsdict dict-like
    ## but in our algorithm only use get() method above    
    def __getitem__(self, key):
        return self.get(key)

    def __setitem__(self, key, value):
        self.D[key] = value

    def __delitem__(self, key):
        del self.D[key]

    def __contains__(self, key):
        if key in self.D:
            return True
        v = self.get(key,None)
        return v is not None

def grouper(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    return it.zip_longest(*args, fillvalue=fillvalue)

def grouper_str(iterable, n, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx"
    args = [iter(iterable)] * n
    for item in it.zip_longest(*args, fillvalue=fillvalue):
        ## each item is a tuple of n values
        yield "".join(item)

DASH = '-'
THREE_DASHES = ('-','-','-')
def nuseq_realign(nuseq,aaseq):
    nuseq_triples = filter(lambda x: x!=THREE_DASHES,
                           grouper(nuseq,3,fillvalue=DASH))
    aaseq_indices = filter(lambda n: aaseq[n] != DASH, range(len(aaseq)))
    nuseq_realigned = [THREE_DASHES] * len(aaseq)
    for n,triple in zip(aaseq_indices,nuseq_triples):
        nuseq_realigned[n] = triple
    flatlist = it.chain.from_iterable(nuseq_realigned) # flatten list of tuples
    return "".join(flatlist)
    
def nuseq_realign_filter(seqs,aaseqs,aanotfound=None,upfront=True):
    ''' realign the nuseqs to match the aaseqs '''

    ## Use aadict to get aa sequence that matches nu sequence
    if upfront:
        ## convert the aaseqs iterable into a dictionary, up front
        aadict = {s.name: s.seq for s in aaseqs}
    else:
        ## use seqdict class, which is dict-like, but loads up only as needed
        aadict = seqdict(aaseqs)
        
    for s in seqs:
        aa = aadict.get(s.name,None)
        if aa is not None:
            s.seq = nuseq_realign(s.seq,aa)
            del aadict[s.name] ## don't need this any more
        else:
            if aanotfound is not None:
                aanotfound.append(s.name)
        yield s

def nuseq_realign_filter2(seqs,aaseqs,aanotfound=None):
    ''' realign the nuseqs to match the aaseqs '''
    ## use 'seqdict' which is like dict, but reads as necessary
    aadict = seqdict(aaseqs)  # {aa.name: aa.seq for aa in aaseqs}
    for s in seqs:
        if s.name in aadict:
            s.seq = nuseq_realign(s.seq,aadict[s.name])
        else:
            if aanotfound is not None:
                aanotfound.append(s.name)
        yield s

def _main(args):
    '''main'''
    vprint(args)
    aaseqs = readseq.read_seqfile(args.aa)
    if args.verbose:
        aaseqs = wrapgen.keepcount(aaseqs,"AA sequences:")
    
    nuseqs = readseq.read_seqfile(args.nu)
    if args.verbose:
        nuseqs = wrapgen.keepcount(nuseqs,"NU sequences:")

    aanotfound=[]
    upfront = True if args.N==0 else False
    seqs = nuseq_realign_filter(nuseqs,aaseqs,aanotfound=aanotfound,upfront=upfront)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"realigned sequences")

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    if args.output:
        readseq.write_seqfile(args.output,seqs)

    if aanotfound:
        print("Did not find",len(aanotfound),"aa sequences")
        print("Including: ",aanotfound[:5])

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
