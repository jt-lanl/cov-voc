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
import verbose as v

import lineagetable

def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--skipnone",action="store_true",
        help="sequences labeled 'None' are totally ignored, not put in OTHER")
    paa("--skipother",action="store_true",
        help="skip all sequences in OTHER category")
    paa("--usehex",action="store_true",
        help="use six-character hex-codes instead of X11 colornames")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    return args

def main(args):
    '''main parselineagetable'''

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)

    first,seqs = sequtil.get_first_item(seqs)
 
    T = lineagetable.get_lineage_table(args.lineagetable)
    
    for s in seqs:
        lineage = covid.get_lineage(s)
        if args.skipnone:
            if lineage in [None,"None","Unassigned",""]:
                v.vprint_only(5,"skip None:",f'[{s.name}]')
                continue
        if lineage is not None:
            voc = T.last_match("."+lineage)
        else:
            voc = lineagetable.OTHER
        if args.skipother and voc == lineagetable.OTHER:
            v.vprint_only(5,"skip Other:",f'[{s.name}]')
            continue
            
        print(s.name,T.names[voc],T.colors[voc])

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
