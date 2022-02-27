'''embers-style analysis/plotting based only on names of sequences'''
from collections import Counter
import argparse

import covid
from verbose import verbose as v
import embersutil as emu
import lineagetable


OTHER = lineagetable.OTHER
DEFAULTNAMESFILE="Latest-names.nm"

def _getargs():
    ap = argparse.ArgumentParser(description=__doc__,
                                 conflict_handler='resolve')
    covid.corona_args(ap)
    ap.set_defaults(input=covid.default_seqfile(DEFAULTNAMESFILE))
    emu.embers_args(ap)
    paa = ap.add_argument
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--other",action="store_true",
        help="write out the lineages in the 'other' class")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    args = ap.parse_args()
    return args

def print_other_lineages(other_lineages):
    '''if other lineages found, print out a summary'''
    otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
    for otherlin in otherlist:
        print("%6d %s" % (other_lineages[otherlin],otherlin))

def main(args):
    '''sparks main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs_by_pattern(seqs,args,keepfirst=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args)
    v.vprint(args)

    T = lineagetable.get_lineage_table(args.lineagetable)

    date_counter = {m: Counter() for m in T.patterns}
    other_lineages = Counter()
    for s in seqs:

        seqdate = emu.date_from_seqname(s.name)
        if not seqdate:
            v.vprint_only(5,"No seqdate:",s.name)
            continue

        if (seqdate.year,seqdate.month) < (2019,11):
            v.vprint_only(5,"Bad seqdate:",seqdate)
            continue

        voc = T.first_match(s.name)
        date_counter[voc][seqdate] += 1

        if args.other and voc == OTHER:
            lineage = covid.get_lineage_from_name(s.name)
            other_lineages[lineage] += 1
            if lineage != 'None':
                v.vvprint_only(10,'Other:',s.name)

    v.vprint_only_summary("No seqdate:","warnings triggered")
    v.vprint_only_summary("Bad seqdate:","warnings triggered")
    v.vprint_only_summary("Other:","non-None sequences in OTHER category")

    if args.other:
        print_other_lineages(other_lineages)

    nmatches = sum(sum(date_counter[p].values()) for p in T.patterns)
    v.vprint("matched sequences:",nmatches)
    if nmatches==0:
        raise RuntimeError("No sequences for: " + " ".join(args.filterbyname))

    onsets=dict()
    if args.onsets:
        ## Don't include OTHER or T.patterns that don't appear in sequence set
        onsets.update( {m: min(date_counter[m])
                        for m in T.patterns[1:]
                        if date_counter[m] and m != OTHER} )

    for line in emu.mk_counts_table(date_counter,T.names):
        v.vprint(line)

    ord_range, ord_plot_range = emu.get_ord_daterange(date_counter,args.dates)
    v.vprint("ordinal range:",ord_range,ord_plot_range)
    cum_counts = emu.get_cumulative_counts(date_counter,ord_range,
                                           daysperweek=args.daily)

    emu.make_emberstyle_plots(args,'bynames',cum_counts,T.names,T.colors,ord_range[0],
                              ordplotrange = ord_plot_range,
                              title=": ".join([covid.get_title(args),
                                               f"{nmatches} sequences"]),
                              onsets=onsets)


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
