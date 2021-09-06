'''Compare two alignments'''
import sys
import itertools as it
from pathlib import Path
import argparse

import readseq
import wrapgen
import covid
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--alignments","-a",nargs=2,type=Path,
        help="two sequence alignment files")
    paa("-N",type=int,default=0,
        help="for debugging purposes, just do the first N sequences")
    paa("--output","-o",nargs=2,type=Path,
        help="write two fasta files with the differing sequences")
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
                self.D[s.name] = s
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

def _main(args):
    '''main'''
    vprint(args)
    afile,bfile = args.alignments
    aseqs = readseq.read_seqfile(afile)
    if args.verbose:
        aseqs = wrapgen.keepcount(aseqs,"A sequences:")
    afirst,aseqs = covid.get_first_item(aseqs)
    a_mgr = mutant.MutationManager(afirst.seq)
    
    bseqs = readseq.read_seqfile(bfile)
    if args.verbose:
        bseqs = wrapgen.keepcount(bseqs,"B sequences:")
    bfirst,bseqs = covid.get_first_item(bseqs)
    b_mgr = mutant.MutationManager(bfirst.seq)

    if args.N:
        aseqs = it.islice(aseqs,args.N+1)

    a_del_seqs=[afirst]
    b_del_seqs=[bfirst]
        
    bdict = seqdict(bseqs)

    set_of_diffs = set()
    
    for a in aseqs:
        b = bdict[a.name]
        amut = a_mgr.get_mutation(a.seq)
        bmut = b_mgr.get_mutation(b.seq)
        if amut == bmut:
            continue
        if (str(amut),str(bmut)) not in set_of_diffs:
            if len(amut) != len(bmut):
                vprint(a.name)
                vprint("   A:",amut)
                vprint("   B:",bmut)
            a_del_seqs.append(a)
            b_del_seqs.append(b)
        set_of_diffs.add( (str(amut),str(bmut)) )

    #for diff in set_of_diffs:
    #    vprint("A:",diff[0])
    #    vprint("B:",diff[1])

    if args.output:
        adel,bdel = args.output
        readseq.write_seqfile(adel,a_del_seqs)
        readseq.write_seqfile(bdel,b_del_seqs)
        
        
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
