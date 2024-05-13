"""
Find parent/child lineage pairs for which new mutations appear

1. From sequence file, obtain all the pango lineages
2. For each lineage, get most common sequence string
3.     Get mutation list (vis-a-vis Wuhan) for that string
4. For each lineage, get parent, and parent's mutation list
5.     Identify new mutations vis-a-vis the parent
6. While doing all this, be collecting a list of mutations
7. For each mutation, identify child/parent pairs, including counts
"""

from operator import attrgetter
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
import argparse

import verbose as v
import breakpipe

import covid
from xopen import xopen
from mutant import SingleSiteMutation, MutationManager
from lineagenotes import LineageNotes
from commonforms import LineagePartition


def getargs():
    """get arguments from command line"""
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    paa("--verbose", "-v", action="count", default=0, help="verbose")
    covid.corona_args(ap)
    paa = ap.add_argument_group("New Mutations Options").add_argument
    paa("--notesfile", help="lineage_notes.txt")
    paa(
        "--cutoff",
        "-c",
        type=int,
        default=0,
        help="Do not count any lineage mutations that occur fewer then c times",
    )
    paa("--clade", help="Restrict attention to sequences in clade's lineages")
    paa("--xclade", help="Avoid sequences in xclade's lineages")
    paa(
        "--skiprecomb",
        action="store_true",
        help="skip recombinant lineages (whose names begin with X)",
    )
    paa("--mutationsfile", "-M", help="Write mutations to a tsv file")  # default="-",
    paa("--reversionsfile", "-R", help="Write reversions to a tsv file")
    paa("--output", "-o", help="Write two tsv files: equivealent to -M new-OUTPUT -R rev-OUTPUT ")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    if args.notesfile is None:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)
    if args.output:
        if not args.mutationsfile:
            args.mutationsfile = covid.filename_prepend("new-", args.output)
        if not args.reversionsfile:
            args.reversionsfile = covid.filename_prepend("rev-", args.output)
    if args.clade and "-" in args.clade:
        if args.xclade:
            raise RuntimeError("Cannot use minus (-) in clade and also have an xclade")
        args.clade, args.xclade = args.clade.split("-")
    if args.clade and args.clade.startswith("X"):
        if args.skiprecomb:
            v.print(f"Setting --skiprecomb=False since --clade={args.clade} is a recombinant")
        args.skiprecomb = False

    return args


def most_common_forms(mmgr: MutationManager, lin_partition: LineagePartition) -> dict[str, set]:
    """return a dict with most common form (mcf) for each lineage,
    with mcf expressed as a set of SingleSiteMutation's"""
    mcf_dict = dict()  ## most common form for each lineage
    for lin, seqlin in lin_partition.sequences.items():
        ## Get most common form
        cntr = Counter(s.seq for s in seqlin)
        mcf, _ = cntr.most_common(1)[0]
        ## Express as mutations
        mcf_muts = mmgr.seq_to_mutation(mcf)
        mcf_dict[lin] = set(mcf_muts)
    return mcf_dict


@dataclass
class Transition:
    """Transition is a single-site mutation between two lineages"""

    parent: str
    child: str
    count: int

    def __str__(self):
        """__str__ for Transiton type"""
        if self.parent == "Wuhan":
            return f"{self.child}({self.count})"
        return f"{self.parent}->{self.child}({self.count})"


@dataclass
class CountAppearances:
    """count appearances of various mutations/reversions"""

    seq_child: Counter = field(default_factory=Counter)  ## how many sequences have mutation
    lin_child: Counter = field(default_factory=Counter)  ## how many disctinct lineages
    lin_parent: Counter = field(default_factory=Counter)  ## how many distinct parents
    mcf_small: Counter = field(default_factory=Counter)  ## how many below-cutoff MCF mutations


def count_appearances(
    transdict: dict[str, list[Transition]], cutoff: int
) -> tuple[Counter, Counter, Counter, Counter]:
    """parse mutsdict dict to count apppearances"""
    ca = CountAppearances()
    for ssm in transdict:
        pset = set()
        for transit in transdict[ssm]:
            if transit.count < cutoff:
                if transit.parent != transit.child:
                    ca.mcf_small[ssm] += 1
                continue
            pset.add(transit.parent)
            ca.seq_child[ssm] += transit.count
            ca.lin_child[ssm] += 1
        ca.lin_parent[ssm] = len(pset)
    return ca


def tab_join_map_str(*pargs) -> str:
    """convert positional arguments (pargs) into a tab separated string"""
    return "\t".join(map(str, pargs))


def write_mutations_summary(
    filename: str,
    transdict: dict[str, list[Transition]],
    count_seq_total: dict[str, int],
    denominator: int,
    cutoff: int,
):
    """write tsv file summarizing the mutations"""
    if not filename:
        return
    ca = count_appearances(transdict, cutoff)

    header = tab_join_map_str(
        "site",
        "mutation",
        "parent_lineages",
        "child_lineages",
        "transitory_lineages",
        "lineage_sequences",
        "total_sequences",
        "denominator",
        "lineage_transitions",
    )

    with xopen(filename, "w") as fout:
        print(header, file=fout)
        v.vprint(f"Writing {len(transdict)} mutations to file={filename}")
        for ssm in sorted(transdict):
            if count_seq_total[ssm] < cutoff:
                continue
            lineage_transitions = ", ".join(
                str(transit)
                for transit in sorted(transdict[ssm], key=attrgetter("count"), reverse=True)
            )
            print(
                tab_join_map_str(
                    ssm.site,
                    ssm,
                    ca.lin_parent[ssm],
                    ca.lin_child[ssm],
                    ca.mcf_small[ssm],
                    ca.seq_child[ssm],
                    count_seq_total[ssm],
                    denominator,
                    lineage_transitions,
                ),
                file=fout,
            )


def get_lineages(
    notesfile: str, clade: str, xclade: str = None, skiprecomb: bool = False
) -> tuple[LineageNotes, set]:
    """return set of lineages in the given clade (minus those in xclade)"""
    lin_notes = LineageNotes.from_file(notesfile)
    lineage_set = lin_notes.get_lineage_set(clade)
    if xclade:
        xlineage_set = lin_notes.get_lineage_set(xclade)
        lineage_set = lineage_set - xlineage_set
    if skiprecomb:
        lineage_set = set(
            lin for lin in lineage_set if not lin_notes.get_fullname(lin).startswith("X")
        )
    v.vprint("Total lineages:", len(lineage_set))
    return lin_notes, lineage_set


def set_filter_out_x(ssms: Iterable[SingleSiteMutation]) -> set[SingleSiteMutation]:
    """return set of ssm's, with X's filtered out"""
    return set(ssm for ssm in ssms if "X" not in str(ssm))


@breakpipe.no_broken_pipe
def main(args: argparse.Namespace):
    """newmuts main"""
    v.vprint(args)

    lin_notes, lineage_set = get_lineages(
        args.notesfile, args.clade, args.xclade, skiprecomb=args.skiprecomb
    )

    seqs = covid.read_filter_seqfile(args)
    first, seqs = covid.get_first_item(seqs)
    mut_manager = MutationManager(first.seq)

    lin_partition = LineagePartition(seqs, restrict_to=lineage_set)
    denominator = sum(lin_partition.counts.values())
    v.vprint(f"denominator={denominator}")
    v.vprint(f"lineages in partition: {len(lin_partition.lineages)}")

    mcf_dict = most_common_forms(mut_manager, lin_partition)
    siteset = set()
    for mcf in mcf_dict.values():
        siteset.update(ssm.site for ssm in mcf if ssm.ref != "+")

    ## lists of lineage Transition's for each ssm
    ssm_appearances = defaultdict(list)
    ssm_reversions = defaultdict(list)
    ## count seqs with ssm, regardless of lin
    ssm_total = Counter()
    site_total_ref = Counter()  ## eg: site_total_ref[373] counts seqs with "S373S"
    for lin in lineage_set:  ## only consider lineages in clade
        if lin not in mcf_dict:
            v.vprint(f"{lin=} not in mcf_dict, skipping...")
            continue
        parent = lin_notes.parent_of(lin)
        if parent not in lineage_set:
            if args.clade not in lineage_set:
                parent = "Wuhan"
            elif lin == args.clade:
                parent = lin
            else:
                raise RuntimeError(f"lineage {lin} has {parent=} not in lineage set")
        parmutset = set_filter_out_x(mcf_dict.get(parent, set()))
        linmutset = set_filter_out_x(mcf_dict[lin])
        seqlin = lin_partition.sequences[lin]
        seq_counter = Counter(s.seq for s in seqlin)
        ## ssm_appear: for seqs in child(lin) population that have ssm...
        ssm_appear_int = Counter()  ## but child MCF does not have ssm
        ssm_appear_ext = Counter()  ## and child MCF also has ssm, but parent MCF does not have ssm
        ## ssm_revert: for seqs inn child(lin) population that do not have ssm...
        ##              AND mut[site] = ref[site]
        ssm_revert_int = Counter()  ## but child MCF does have ssm
        ssm_revert_ext = Counter()  ## and child MCF also lacks ssm, but parent MCF does have ssm
        for seq, cnt in seq_counter.items():
            mut = mut_manager.seq_to_mutation(seq)
            # mut_ref_sites: set of sites NOT in mut; so have Wuhan ref values
            mut_ref_sites = siteset - set(ssm.site for ssm in mut if ssm.ref != "+")
            for site in mut_ref_sites:
                site_total_ref[site] += cnt
            mutset = set_filter_out_x(mut)
            for ssm in mutset:
                ssm_total[ssm] += cnt

            for ssm in mutset - linmutset:
                ssm_appear_int[ssm] += cnt
            for ssm in (mutset & linmutset) - parmutset:
                ssm_appear_ext[ssm] += cnt
            for ssm in linmutset - mutset:
                if ssm.site in mut_ref_sites:
                    ssm_revert_int[ssm] += cnt
            for ssm in parmutset - (mutset | linmutset):
                if ssm.site in mut_ref_sites:
                    ssm_revert_ext[ssm] += cnt

        for ssm, mcnt in ssm_appear_int.items():
            if mcnt >= args.cutoff:
                ssm_appearances[ssm].append(Transition(lin, lin, mcnt))
        for ssm, mcnt in ssm_appear_ext.items():
            ## no cutoff for external transitions
            ssm_appearances[ssm].append(Transition(parent, lin, mcnt))
        for ssm, mcnt in ssm_revert_int.items():
            if mcnt >= args.cutoff:
                ssm_reversions[ssm].append(Transition(lin, lin, mcnt))
        for ssm, mcnt in ssm_revert_ext.items():
            ssm_reversions[ssm].append(Transition(parent, lin, mcnt))

    ssm_total_reversions = {ssm: site_total_ref[ssm.site] for ssm in ssm_reversions}

    write_mutations_summary(
        args.mutationsfile, ssm_appearances, ssm_total, denominator, cutoff=args.cutoff
    )
    write_mutations_summary(
        args.reversionsfile, ssm_reversions, ssm_total_reversions, denominator, cutoff=args.cutoff
    )


if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
