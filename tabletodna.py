'''
Input is tsv file that is crated by mktable and a DNA fasta file.
Output is DNA fasta file with the ISL numbers indicated in the table,
using names also provided in the table
'''


import re
from collections import Counter
import argparse
import pandas as pd
from dataclasses import dataclass

import warnings

import verbose as v
import breakpipe
import xopen
import sequtil
import covid
import numu

def _getargs():
    '''get arguments from command line'''
    aparser = argparse.ArgumentParser(description=__doc__)
    paa = aparser.add_argument
    paa("--tablefile", "-t",required=True,
        help="tsv file created by mktable")
    paa("--dnainputfile","-d",required=True,
        help="sequence file with reference (eg, DNA) sequences")
    paa("--output","-o",
        help="output fasta file with reference sequences")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = aparser.parse_args()
    return args

@dataclass
class SeqInfo:
    who: str
    pango: str
    isl: str


def get_row_item(row,shortkey):
    fullkey = numu.match_column_name(row.keys(),shortkey)
    if not fullkey:
        raise RuntimeError(f'Cannot match {shortkey} among {row.keys()}')
    item = row[fullkey]
    try:
        item = item.strip()
    except AttributeError:
        pass
    return item

def read_tablefile(filename):
    '''read mutant file and return mstring list'''
    ## assume tab separated columns
    df = pd.read_csv(filename,sep='\t')
    for _,row in df.iterrows():
        who = get_row_item(row,"WHO")
        pango = get_row_item(row,"Pango lineage")
        isl = get_row_item(row,"ISL")
        yield SeqInfo(who,pango,isl)

def get_refseqdict(dnainputfile,seqinfolist):
    '''read DNA file, make a dict of sequences indexed by ISL number,
    only including those in the isl_matches dictionary'''
    isl_setofall = set()
    for seqinfo in seqinfolist:
        isl_setofall.add(seqinfo.isl)
    v.vprint('ISL set:',isl_setofall)

    refseqs = sequtil.read_seqfile(dnainputfile)
    firstref,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    refseqdict = dict()
    for s in refseqs:
        isl_name = covid.get_isl(s)
        if isl_name in isl_setofall:
            refseqdict[isl_name] = s
    v.vprint("Read",len(refseqdict),"reference sequences")
    return firstref,refseqdict

@breakpipe.no_broken_pipe
def _main(args):
    '''mutlineage main'''

    seqinfolist = list(read_tablefile(args.tablefile))
    v.vprint(seqinfolist)
    firstref,refseqdict = get_refseqdict(args.dnainputfile,seqinfolist)

    dna_seqs = []
    if firstref:
        dna_seqs.append(firstref)

    for seqinfo in seqinfolist:
        if seqinfo.isl not in refseqdict:
            warnings.warn(f"Did not find seq for ISL={seqinfo.isl}")
            continue
        dnaseq = refseqdict[seqinfo.isl]
        sample_date = covid.get_date(dnaseq.name,as_string=True,check_date=False)
        new_name = f'{seqinfo.who}_{seqinfo.pango}_{seqinfo.isl}_{sample_date}'
        dna_seqs.append( sequtil.SequenceSample(new_name,dnaseq.seq) )

    if args.output:
        sequtil.write_seqfile(args.output,dna_seqs)



if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
