'''embers-style analysis/plotting based only on names of sequences'''
from collections import Counter
import argparse
import re
import datetime

import covid
from verbose import verbose as v
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
    ap.set_defaults(input=covid.default_seqfile(DEFAULTNAMESFILE))
    emu.embers_args(ap)
    paa = ap.add_argument
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--writeother",
        help="write out the lineages in the 'other' class to this file")
    paa("--skipnone",action="store_true",
        help="sequences labeled 'None' are totally ignored, not put in OTHER")
    paa("--skipother",action="store_true",
        help="skip all sequences in OTHER category")
    paa("--cases",#default="data/owid-harmonized.csv",
        help="csv file with case counts (harmonized, from OWID)")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    args = ap.parse_args()
    return args


def print_other_lineages(filename,other_lineages):
    '''if other lineages found, print out a summary'''
    if not filename:
        return
    otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
    with open(filename,'w') as fout:
        for lin in otherlist:
            print("%6d %s" % (other_lineages[lin],lin),file=fout)

def main(args):
    '''sparks main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs_by_pattern(seqs,args,keepfirst=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args)
    v.vvprint(args)

    if args.lineagetable and args.lineagetable.upper() == "FROMSEQS":
        ## special case in which lineage table is "created on the fly"
        ## with the most common pango lineages (with no grouping of
        ## similar lineages as is done inthe lineage table files)
        seqs = list(seqs)
        T = lineagetable.get_lineage_table_from_seqs(seqs,
                                                     num_lineages=24,
                                                     skipnone=args.skipnone)
    else:
        ## build lineage table from the file
        T = lineagetable.get_lineage_table(args.lineagetable)

    v.vvprint('patterns',T.patterns)
    v.vvprint('names',list(T.names.values()))

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

        lineage = covid.get_lineage_from_name(s.name)

        if args.skipnone:
            if lineage in ["None","Unassigned",""]:
                v.vprint_only(5,"skip None:",f'[{s.name}]')
                continue

        voc = T.last_match("."+lineage)

        if voc == OTHER:
            if args.writeother:
                other_lineages[lineage] += 1
                v.vvprint_only(10,'Other:',s.name)
            if args.skipother:
                v.vprint_only(5,"skip other:",voc,s.name)
                continue

        date_counter[voc][seqdate] += 1

    v.vprint_only_summary("No seqdate:","warnings triggered")
    v.vprint_only_summary("Bad seqdate:","warnings triggered")
    v.vprint_only_summary("Other:","sequences in OTHER category")
    v.vprint_only_summary("skip None:","sequences skipped")
    v.vprint_only_summary("skip other:","sequences skipped")

    v.vprint('Other lineages:',sum(other_lineages.values()))
    v.vprint('OTHER lineages:',sum(date_counter[OTHER].values()))

    print_other_lineages(args.writeother,other_lineages)

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
    v.vprint("ordinal to date:",args.dates,
             datetime.date.fromordinal(ord_range[0]),
             datetime.date.fromordinal(ord_range[1]))
    cum_counts = emu.get_cumulative_counts(date_counter,ord_range,
                                           daysperweek=args.daily)

    num_cases=None
    ncases=0
    if args.cases:
        df = owid.read_dataframe(args.cases)
        if df is None:
            raise RuntimeError(f'Cannot read OWID datafile: {args.cases}')
        df = owid.filter_cases(df,args.filterbyname,args.xfilterbyname)
        num_cases,ncases = owid.case_counts(df,ord_plot_range,
                                            daysperweek=args.daily)
        if max(num_cases) == 0:
            v.vprint(f"No OWID case data for {args.filterbyname}")
            num_cases = None

    if args.skipother:
        del cum_counts[OTHER]
        T.del_pattern(OTHER)

    if args.skipnone:
        for nonestring in ['None','Unssigned',
                           '(None|Unassigned)',
                           '(None|Unassigned|)']:
            try:
                del cum_counts[nonestring]
                T.del_pattern(nonestring)
                v.vvprint('Removing:',nonestring)
            except:
                v.vvprint('Cannot remove:',nonestring)

    emu.make_emberstyle_plots(args,'bynames',cum_counts,
                              T.names,T.colors,ord_range[0],
                              ordplotrange = ord_plot_range,
                              num_cases = num_cases,
                              title=covid.get_title(args),
                              nmatches=nmatches,
                              ncases=ncases,
                              daily=args.daily,
                              onsets=onsets)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
