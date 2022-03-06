'''Make a keyfile that is used by mkkey to make a pdf key table'''

import argparse
import matplotlib.pyplot as plt

from verbose import verbose as v
from spikevariants import SpikeVariants
import sequtil
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--colormut","-c",
        help="color mutation table")
    paa("--view",type=int,default=1,choices=[1,2,3],
        help="view: 1, 2, or 3")
    paa("--pdf",
        help="write key as a pdf file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def mk_key_figure(colors,names,output_pdf=None):
    '''make a pdf figure'''
    maxnamelen = max(len(name) for name in names)

    plt.figure(figsize=((7+maxnamelen)/11,len(colors)/4.5))
    for color,name in zip(colors,names):
        plt.bar([0],[0],color=color,label=name)
        plt.bar([0],[1],bottom=[-1],color='white') #blank it out!

    plt.axis('off')
    plt.legend(bbox_to_anchor=(1.02, 1),
               #handlelength=3,
               #markerfirst=False,
               frameon=False,
               handletextpad=1,
               labelspacing=0.45,
               prop={'family' : 'monospace'})
    plt.tight_layout()
    if output_pdf:
        if not output_pdf.endswith(".pdf"):
            output_pdf += ".pdf"
        plt.savefig(output_pdf)
    else:
        plt.show()

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs)

    if args.colormut:
        svar = SpikeVariants.from_colormut(args.colormut,refseq=first.seq)
    else:
        svar = SpikeVariants.default(refseq=first.seq)

    if not args.pdf:
        svar.key_print(args.view,seqs=seqs)
    else:
        keylines = svar.key_view(args.view,seqs=seqs)
        colors_and_labels = [line.split(' ',1) for line in keylines]
        colors = [cn[0] for cn in colors_and_labels]
        labels = [cn[1] for cn in colors_and_labels]
        mk_key_figure(colors,labels,args.pdf)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
