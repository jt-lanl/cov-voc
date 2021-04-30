DESCRIPTION='''
USITES: union over all continents of high-entropy sites 
'''

import sys
from collections import Counter
from pathlib import Path
import numpy as np
import itertools as it
import scipy.stats as sst

import argparse


import readseq
import sequtil
import intlist
import covid

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
#                                 conflict_handler='resolve')
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--sites",type=int,default=30,
        help="Number of highest-entropy sites to use")
    paa("--keepx",action="store_true",
        help="Keep sequences that have X's in them")
    paa("--stripdashcols",
        help="Strip dashes from master, and align rest, saving to specified file")
    paa("--restrictsites",
        help="Consider only these sites (RBD, NTD, NTD+RBD, or comma-separated list)")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
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

def main(args):

    allseqs = covid.read_seqfile(args)
    allseqs = covid.filter_seqs_by_date(allseqs,args)

    firstseq = allseqs[0].seq
    if "-" in firstseq:
        vprint("Stripping sites with dashes in first sequence...",end="")
        sequtil.stripdashcols(firstseq,allseqs)
        vprint("ok")

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
            vprint(f"No sequences in {cx}")
            esites[cx]=[]
            continue

        N = len(seqs)
        M = len(firstseq)

        #### SINGLE-SITE ENTROPY
        vprint(f"Single-site entropy {cx}...",end="")
        E = [xentropy(clist,keepx=args.keepx)
             for clist in sequtil.gen_columns_seqlist(seqs)]
        vprint("ok")

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
        vprint(n,len(esitelist))
        if len(esitelist) >= args.sites:
            break
        
    print(",".join(str(e) for e in esitelist))

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose and args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)

  

