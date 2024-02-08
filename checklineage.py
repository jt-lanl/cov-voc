'''from a lineage table and a names file, make a list of full pango names associated with each lineage'''
from collections import Counter,defaultdict
import argparse
import re

import covid
import verbose as v
import embersutil as emu
import lineagetable
import owid

OTHER = lineagetable.OTHER
EMPTY_LINEAGE_REGEX = re.compile(r'EPI_ISL_\d+\.\s*$')
DEFAULTNAMESFILE="Latest-names.nm"

def _getargs():
    ap = argparse.ArgumentParser(description=__doc__,
                                 conflict_handler='resolve')
    covid.corona_args(ap)
    ap.set_defaults(input=covid.find_seqfile(DEFAULTNAMESFILE))
    emu.embers_args(ap) ## between this and corona_args, total overkill!
    paa = ap.add_argument
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    args = ap.parse_args()
    return args

def main(args):
    '''sparks main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args)
    ## for names file, the reference seq has already been stripped (so firstisref=False)
    seqs = covid.filter_seqs_by_pattern(seqs,args,firstisref=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args,firstisref=False)
    v.vvprint(args)

    T = lineagetable.get_lineage_table(args.lineagetable)

    v.vvprint('patterns',T.patterns)
    v.vvprint('names',list(T.names.values()))

    pango_names_by_lineage = defaultdict(set)
    pango_counts = Counter()
    other_lineages = Counter()
    for s in seqs:

        pango = covid.get_lineage(s)
        voc = T.last_match("."+pango)
        pango_names_by_lineage[voc].add(pango)
        pango_counts[pango] += 1

    for patt in T.patterns:
        print("Name:",T.names[patt])
        print("Patt:",patt)
        for pango in sorted(pango_names_by_lineage[patt],
                            key=pango_counts.get,
                            reverse=True):
            print("      %6d %s" % (pango_counts[pango],pango))

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
