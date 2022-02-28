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
    paa("--other",
        help="write out the lineages in the 'other' class to this file")
    paa("--cnpo",nargs=3,
        metavar=('COLOR','NAME','POSN'),
        help="color, name, and position (0=bottom) to be used for OTHER")
    ## typical use: --cnpo Yellow D614G 1
    ## treats all 'other' as D614G and then moves it up one notch in the plot,
    ## so that it comes after the 'Ancestral'
    ## or if 'Black None None' is in the lineage table as well,
    ## you might do --cnpo Yellow D614G 2, to put it above None and Ancestral
    paa("--skipnone",action="store_true",
        help="sequences labeled 'None' are totally ignored, not put in OTHER")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    args = ap.parse_args()

    ## should we check that args.dates is in correct order?
    if args.dates:
        if not any(bool("." in date) for date in args.dates):
            if args.dates[0] > args.dates[1]:
                raise RuntimeError(f'out of order --dates {args.dates}')

    return args

def print_other_lineages(filename,other_lineages):
    '''if other lineages found, print out a summary'''
    if not filename or not other_lineages:
        return
    otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
    with open(filename,'w') as fout:
        for otherlin in otherlist:
            print("%6d %s" % (other_lineages[otherlin],otherlin),file=fout)

def main(args):
    '''sparks main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs_by_pattern(seqs,args,keepfirst=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args)
    v.vprint(args)

    T = lineagetable.get_lineage_table(args.lineagetable,
                                       other=args.cnpo)

    v.vprint('T',T)
    v.vprint('patterns',T.patterns)
    v.vprint('names',list(T.names.values()))

    date_counter = {m: Counter() for m in T.patterns}
    other_lineages = Counter()
    for s in seqs:

        seqdate = emu.date_from_seqname(s.name)
        if not seqdate:
            v.vprint_only(5,"No seqdate:",s.name)
            continue

        if (seqdate.year,seqdate.month) < (2019,11):
            v.vprint_only(5,"Bad seqdate:",seqdate,s.name)
            continue

        if args.skipnone and "None" in s.name:
            v.vprint_only(5,"None:",f'[{s.name}]')
            continue

        voc = T.first_match(s.name)
        date_counter[voc][seqdate] += 1

        if args.other and voc == OTHER:
            lineage = covid.get_lineage_from_name(s.name)
            other_lineages[lineage] += 1
            v.vvprint_only(10,'Other:',s.name)

    v.vprint_only_summary("No seqdate:","warnings triggered")
    v.vprint_only_summary("Bad seqdate:","warnings triggered")
    v.vprint_only_summary("Other:","sequences in OTHER category")
    v.vprint_only_summary("None:","sequences skipped")

    v.vprint('Other lineages:',sum(other_lineages.values()))
    v.vprint('OTHER lineages:',sum(date_counter[OTHER].values()))

    print_other_lineages(args.other,other_lineages)

    nmatches = sum(sum(date_counter[p].values()) for p in T.patterns)
    v.vprint("matched sequences:",nmatches)
    if nmatches==0:
        ifilters = " ".join(args.filterbyname) if args.filterbyname else "Global"
        xfilters = " w/o " + " ".join(args.xfilterbyname) if args.xfilterbyname else ""
        raise RuntimeError(f"No sequences for: {ifilters}{xfilters}")

    onsets=dict()
    if args.onsets:
        ## Don't include OTHER or T.patterns that don't appear in sequence set
        onsets.update( {m: min(date_counter[m])
                        for m in T.patterns
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
