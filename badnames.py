'''Find bad names'''
import sys
import re
import itertools as it
import argparse

import sequtil
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(argparser)
    paa = argparser.add_argument
    paa("-N",type=int,default=0,
        help="only check first N sequences")
    paa("--checkdates",action="store_true",
        help="check dates, too")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def lin_regexp(name):
    return re.sub(r'.*EPI_ISL_\d+\.','',name)

def lin_bydots(name):
    try:
        tokens = name.split('.',6)
        lintok = tokens[6]
    except IndexError:
        lintok = None
    return lintok

def date_regexp(name):
    m = re.search(r'\.(20\d\d(-\d\d)?(-\d\d)?)',name)
    return m[1] if m else None

def date_bydots(name):
    try:
        tokens = name.split('.')
        datestring = tokens[4]
    except IndexError:
        datestring = None
    return datestring

def _main(args):
    '''badnames main'''
    vprint(args)

    seqs = covid.read_filter_seqfile(args)
    _,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    if args.N:
        seqs = it.islice(seqs,args.N)

    for s in seqs:
        ldot = lin_bydots(s.name)
        lreg = lin_regexp(s.name)
        if ldot != lreg:
            print(s.name,"Lineage:",ldot,lreg)
            continue

        if args.checkdates:
            ddot = date_bydots(s.name)
            dreg = date_regexp(s.name)
            if ddot != dreg:
                print(s.name,"Dates:",ddot,dreg)

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
