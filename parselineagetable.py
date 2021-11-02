'''
create output text with three columns: sequence name, lineage designation, color
'''

import sys
import itertools as it
import re
import argparse

import warnings

import sequtil
import covid
import colornames

from lineagetablev1 import LineageTable

def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--usehex",action="store_true",
        help="use six-character hex-codes instead of X11 colornames")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def main(args):
    '''main parsecolormut'''

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    first,seqs = sequtil.get_first_item(seqs)

    table = [
        ('Gray', 'other', 'other'),
        ] + LineageTable

    patterns = []
    colors = dict()
    fullnames = dict()

    for color,name,patt in table:
        colors[patt] = colornames.tohex(color) if args.usehex else color
        patterns.append(patt)
        fullnames[patt] = name
    
    for s in seqs:

        vocmatch = [voc for voc in patterns[1:]
                    if re.search(r'\.'+voc+'$',s.name)]

        voc = vocmatch[0] if vocmatch else 'other'
        col = colors[voc]
        nom = fullnames[voc]

        print(s.name,nom,col)

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

    main(_args)
