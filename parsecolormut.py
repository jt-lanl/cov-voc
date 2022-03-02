'''
parse the color mutation table to provide a table for Werner
that has three columns: sequence name, mutant designation, color
'''

import itertools as it
import argparse

from verbose import verbose as v
from spikevariants import SpikeVariants
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
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args


def mk_table(svar,seqs):
    '''yield lines of the table: seq_name, variant_name, color'''
    for s in seqs:
        vocs = svar.vocmatch(s.seq)
        if len(vocs) == 1:
            voc_name = vocs[0].name
            voc_color = vocs[0].color
        elif len(vocs) == 0:
            voc_name = "other"
            voc_color = "#808080" #corresponds to Gray om X11
        else:
            voc_names = [voc.name for voc in vocs]
            if any(vname != voc_names[0] for vname in voc_names):
                v.vprint_only(10,'multiple matches:',voc_names)
            voc_name = vocs[0].name
            voc_color = vocs[0].color

        ## color names should already be hex strings,
        ## but just to make sure:
        voc_color = colornames.tohex(voc_color)

        yield (s.name,voc_name,voc_color)

    v.vprint_only_summary('multiple matches:',
                          'sequences with multiple matches')

def main(args):
    '''main parsecolormut'''

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)

    if args.N:
        seqs = it.islice(seqs,args.N+1)

    first,seqs = sequtil.get_first_item(seqs)

    svar = SpikeVariants.from_colormut(args.colormut,refseq=first.seq)

    for seqname,name,color in mk_table(svar,seqs):
        print(seqname,name,color)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
