'''
USITES: union over all continents of high-entropy sites 
'''

import sys
from collections import Counter
from pathlib import Path
import numpy as np
import itertools as it
import scipy.stats as sst

import argparse

import verbose as v
import readseq
import sequtil
import intlist
import covid

def _getargs():
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--sites",type=int,default=30,
        help="Number of highest-entropy sites to use")
    paa("--restrictsites",
        help=("Consider only these sites "
              "(RBD, NTD, NTD+RBD, or comma-separated list)"))
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    return args

def xentropy(clist,keepx=False):
    cnt = Counter(clist)
    if not keepx:
        cnt.pop('X',None)
    return sst.entropy(list(cnt.values()))

def stripxs(alist,blist,badchar='X'):
    '''given a pair of lists (or strings), strip element n where 
    either alist[n]=='X' or blist[n]=='X'
    return pair of tuples
    '''
    ab = [ (a,b) for a,b in zip(alist,blist)
           if a!=badchar and b!=badchar ]
    alist,blist = zip(*ab)
    return alist,blist

def _main(args):

    allseqs = covid.read_filter_seqfile(args)
    allseqs = list(allseqs)

    firstseq = allseqs[0].seq

    esites=dict()
    
    for cx,c,x in covid.parse_continents(withglobal=True):

        seqs = []
        for s in allseqs[1:]:
            if x and x in s.name:
                continue
            if c=="Global" or c in s.name:
                seqs.append(s)

        if len(seqs)==0:
            v.vprint(f"No sequences in {cx}")
            esites[cx]=[]
            continue

        N = len(seqs)
        M = len(firstseq)

        #### SINGLE-SITE ENTROPY
        v.vprint(f"Single-site entropy {cx}...",end="")
        E = [xentropy(clist,keepx=args.keepx)
             for clist in sequtil.gen_columns_seqlist(seqs)]
        v.vprint("ok")

        etop = [e+1 for e in np.argsort(E)[::-1]]
        if args.restrictsites:
            rsites = covid.spike_sites(args.restrictsites)
            etop = [e for e in etop if e in rsites]

        esites[cx] = etop

    for n in range(1,args.sites):
        esitelist=[]
        for cx,c,x in covid.parse_continents(withglobal=True):
            esitelist.extend( esites[cx][:n] )

        esitelist = sorted(set(esitelist))
        v.vprint(n,len(esitelist))
        if len(esitelist) >= args.sites:
            break
        
    print(",".join(str(e) for e in esitelist))

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)

  

