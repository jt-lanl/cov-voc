'''Create a dataframe with mutations as columns and sequences as rows;
each individual cell will have the codon from the row-sequence
associated iwth the column-mutation'''

import argparse
from operator import attrgetter
from pathlib import Path
import pandas as pd
import verbose as v
import sequtil
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="Input DNA sequences, with filnames indicating mutations")
    paa("--output","-o",
        help="tsv file with sequence names, mutations, and codons")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    if not Path(args.input).is_file():
        v.print(f'File {args.input} not found; will not produce {args.output}')
        return

    ## first, read that DNA file
    seqs = sequtil.read_seqfile(args.input)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    seqs=list(seqs)

    ## Get mutations from the sequence names
    mutations = set()
    for s in seqs:
        tokens = s.name.split("_")
        #lin,parent = tokens[0:2]
        #isl = "_".join(tokens[2:5])
        mutations.update(tokens[5:])
    mutations.discard("MCF")
    v.vprint("Mutations:",len(mutations))
    mutations = map(mutant.SingleSiteMutation,mutations)
    mutations = sorted(mutations,key=attrgetter("site","ref","mut"))
    v.vprint("Mutations:",len(mutations),
            " ".join(str(m) for m in mutations[:5]),"...")

    sitelist = sorted(set(map(attrgetter("site"),mutations)))
    v.vprint('Sites:',sitelist[:5],'...')
    sndxtrans = mutant.SiteIndexTranslator(first.seq)
    ndxdict = {str(mut): sndxtrans.index_from_site(3*mut.site-2)
               for mut in mutations}

    v.vvprint(f'ndxdict={ndxdict}')

    SEQNAME="Lineage_Parent_EPI_ISL_number_[MCF]?_Mutations"

    df = pd.DataFrame( columns = mutations + [SEQNAME],
                       #index = range(len(seqs)),
                       dtype=str)
    v.vprint(f'{df}')
    for nr,s in enumerate(seqs):
        for mut in mutations:
            ndx = ndxdict[str(mut)]
            asterisk="*" if str(mut) in s.name else ""
            df.loc[nr,mut] = s.seq[ndx:ndx+3]+asterisk
        df.loc[nr,SEQNAME] = s.name

    v.vprint(f'{df}')

    df.to_csv(args.output,sep='\t',index=False)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
