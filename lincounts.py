"""Provide counts for all sublineages of a specified lineage"""

from collections.abc import Iterator
import argparse
import verbose as v
import breakpipe

import covid
from lineagenotes import LineageNotes
from commonforms import LineagePartition


def _getargs():
    """parse options from command line"""
    argparser = argparse.ArgumentParser(description=__doc__)
    generic_paa = argparser.add_argument
    generic_paa("--verbose", "-v", action="count", default=0, help="verbose")
    covid.corona_args(argparser)
    paa = argparser.add_argument_group("Program Options").add_argument
    paa("--clade", help="include all sublineages of this clade")
    paa("--xclade", help="exclude all sublineages of this clade")
    paa("--notesfile", help="lineage_notes.txt")

    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    if not args.notesfile:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)
        v.vprint(f"Using notes file: {args.notesfile}")

    return args


def order_lineages(lin_notes: LineageNotes, clade: str, restrict_to: set = None) -> Iterator[str]:
    """yield lineages in the order of a depth-first tree"""
    yield clade
    for lin in sorted(lin_notes.children_of(clade)):
        if restrict_to and lin not in restrict_to:
            continue
        yield from order_lineages(lin_notes, lin, restrict_to=restrict_to)


@breakpipe.no_broken_pipe
def _main(args: argparse.Namespace):
    """main"""
    v.vprint("args:", args)
    lin_notes = LineageNotes.from_file(args.notesfile)
    lineage_set = lin_notes.get_lineage_set(args.clade)
    if args.xclade:
        xlineage_set = lin_notes.get_lineage_set(args.xclade)
        lineage_set = lineage_set - xlineage_set
    v.vprint("Lineages in clade:", len(lineage_set))

    seqs = covid.read_filter_seqfile(args)
    seqs = list(seqs)
    lin_partition = LineagePartition(seqs, restrict_to=lineage_set)

    args.days = 60
    seqs60 = covid.filter_seqs(seqs, args)
    v.vprint("Range of dates", args.dates)
    lin_part60 = LineagePartition(seqs60, restrict_to=lineage_set)

    init_offset = lin_notes.get_fullname(args.clade).count(".")
    for lin in order_lineages(lin_notes, args.clade, restrict_to=lineage_set):
        full = lin_notes.get_fullname(lin)
        offset = full.count(".") - init_offset
        count = lin_partition.counts.get(lin, 0)
        count60 = lin_part60.counts.get(lin, 0)
        xlin = " " * offset + lin
        print(f"{xlin:20s} {count:6d} {count60:6d}  {lin_notes.describe.get(lin)}")


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
