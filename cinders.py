'''variant of embers that considers both lineage AND sequence'''
from collections import Counter
import argparse

import covid
import colornames
import sequtil
import mutant
from verbose import verbose as v
import embersutil as emu

import lineagetable #_cinders as lineagetable

OTHER = lineagetable.OTHER
MUTATED = 'mutated'

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__, conflict_handler='resolve')
    paa = ap.add_argument
    covid.corona_args(ap)  ## except want default to be names, not sequences!
    emu.embers_args(ap)
    ## Back to std options
    paa("--lineagetable","-l",
        help="read lineage table from file")
    paa("--mutant","-m",default='A222V',
        help="mutant to consider as variant of lineage")
    paa("--mincount",type=int,default=0,
        help="minimum count needed for inclusion in plot/legend")
    paa("--other",action="store_true",
        help="write out the lineages in the 'other' class")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    ap.set_defaults(keepx=True) ## no reason to exclude X's here
    args = ap.parse_args()
    emu.check_dates_order(args.dates)
    return args

def main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    seqs = covid.filter_seqs_by_pattern(seqs,args,keepfirst=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args)
    v.vprint(args)
    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)
    seqs = sequtil.checkseqlengths(seqs)

    T = lineagetable.get_lineage_table(args.lineagetable)

    for m in T.patterns:
        T.names[m+MUTATED] = T.names[m]+'+'+args.mutant
        T.colors[m+MUTATED] = colornames.darker(T.colors[m])

    ssm = mutant.SingleSiteMutation(args.mutant) ## eg '[A222V]'
    ndx = mutant.SiteIndexTranslator(first.seq).index_from_site(ssm.site)
    def seq_mut_consistent(seq):
        return seq[ndx]==ssm.mut

    date_counter = dict()
    for m in T.patterns:
        date_counter[m] = Counter()
        date_counter[m+MUTATED] = Counter()
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
        if seq_mut_consistent(s.seq):
            date_counter[voc+MUTATED][seqdate] += 1
        else:
            date_counter[voc][seqdate] += 1

        if args.other and voc == OTHER:
            lineage = covid.get_lineage_from_name(s.name)
            other_lineages[lineage] += 1
            v.vprint_only(5,'other:',s.name)

    v.vprint_only_summary("No seqdate:")
    v.vprint_only_summary("Bad seqdate:")
    v.vprint_only_summary("other:")

    if args.other:
        otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
        for otherlin in otherlist:
            print("%6d %s" % (other_lineages[otherlin],otherlin))

    nmatches = sum(sum(date_counter[p].values()) for p in T.patterns)
    v.vprint("matched plain  sequences:",nmatches)
    mut_nmatches = sum(sum(date_counter[p+MUTATED].values()) for p in T.patterns)
    v.vprint("matched mutant sequences:",mut_nmatches)
    if nmatches + mut_nmatches == 0:
        raise RuntimeError("No sequences for: " + " ".join(args.filterbyname))

    onsets=dict()
    if args.onsets:
        ## Don't include OTHER or patterns that don't appear in sequence set
        onsets.update( {m: min(date_counter[m])
                        for m in T.patterns[1:] if date_counter[m]} )

    total_counts=dict()
    for p,counter in date_counter.items():
        total_counts[p]=sum(counter.values())
        v.vvprint('Total counts: %5d %s' % (total_counts[p],p))

    if args.verbose:
        ## make a little table of counts and onsets
        for line in emu.mk_counts_table(date_counter,T.names):
            v.vprint(line)


    ord_range, ord_plot_range = emu.get_ord_daterange(date_counter,args.dates)
    v.vprint("ordinal range:",ord_range,ord_plot_range)
    cum_counts = emu.get_cumulative_counts(date_counter,ord_range,
                                           daysperweek=args.daily)

    for k,tcounts in total_counts.items():
        if MUTATED in k and tcounts <= args.mincount:
            v.vvprint('pop key:',k)
            cum_counts.pop(k)

    emu.make_emberstyle_plots(args,args.mutant,cum_counts,T.names,T.colors,ord_range[0],
                              ordplotrange = ord_plot_range,
                              title = ": ".join([covid.get_title(args),
                                                 f"{nmatches} sequences"]),
                              onsets=onsets)

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
