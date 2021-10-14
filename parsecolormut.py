'''
parse the color mutation table to provide a table for Werner
that has three columns: sequence name, mutant designation, color
'''

import sys
import itertools as it
import argparse

import warnings

from spikevariants import SpikeVariants
import wrapgen
import sequtil
import covid
import colornames

def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--colormut","-c",required=True,
        help="name of color mutation file")
    paa("-N",type=int,default=0,
        help="show at most this many sequences")
    paa("--usehex",action="store_true",
        help="use six-character hex-codes instead of X11 colornames")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def mk_table(svar,seqs,usehex=False):
    '''yield lines of the table: seq_name, variant_name, color'''
    for s in seqs:
        vocs = svar.vocmatch(s.seq)
        if len(vocs) == 1:
            voc_name = vocs[0].name
            voc_color = vocs[0].color
        elif len(vocs) == 0:
            voc_name = "other"
            voc_color = "Gray"
        else:
            voc_names = [voc.name for voc in vocs]
            if any(vname != voc_names[0] for vname in voc_names):
                for voc in vocs:
                    print(voc,voc.name,file=sys.stderr)
                warnings.warn("multiple patterns matched!")
            voc_name = vocs[0].name
            voc_color = vocs[0].color

        if usehex:
            voc_color = colornames.tohex(voc_color) or voc_color
            if voc_name == "other":
                voc_color = "#DDDDDD"  #a slightly different shade of Gray

        yield (s.name,voc_name,voc_color)

def main(args):
    '''main parsecolormut'''

    seqs = covid.read_seqfile(args)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences read:")
    seqs = covid.filter_seqs(seqs,args)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences after filtering:")

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    first,seqs = sequtil.get_first_item(seqs)

    svar = SpikeVariants.from_colormut(args.colormut,refseq=first.seq)

    for seqname,n,c in mk_table(svar,seqs,usehex=args.usehex):
        print(seqname,n,c)

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
