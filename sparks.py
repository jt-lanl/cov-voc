'''embers-style analysis/plotting based only on names of sequences'''
from collections import Counter
import argparse

import covid
import verbose as v
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

def main(args):
    '''sparks main'''
    v.vprint(args)

    seqs = covid.read_seqfile(args) ## <- should keepfirst be here??
    ## filter by country, then by date
    seqs = covid.filter_seqs_by_pattern(seqs,args)
    ## keepfirst=False because names files don't begin with a reference name
    seqs = covid.filter_seqs_by_date(seqs,args,keepfirst=False)
    seqs = covid.fix_seqs(seqs,args)  ## make any sense for names???

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

        ## voc is the first pattern that matches
        voc = T.match_name(s.name)
        date_counter[voc][seqdate] += 1

        if args.other and voc == OTHER:
            lineage = covid.get_lineage_from_name(s.name)
            other_lineages[lineage] += 1
            if lineage != 'None':
                v.vvprint('other:',s.name)

    if args.other:
        otherlist = sorted(other_lineages,key=other_lineages.get,reverse=True)
        for otherlin in otherlist:
            print("%6d %s" % (other_lineages[otherlin],otherlin))

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

    if args.verbose:
        ## make a little table of counts and onsets
        maxnamelen = max(len(T.names[p]) for p in date_counter)
        fmt = f'%{maxnamelen}s %6s %10s %s'
        v.vprint(fmt % ('Name','Count','Onset','Pattern'))
        for p in date_counter:
            name = T.names[p]
            count = sum(date_counter[p].values())
            onset = min(date_counter[p]) if date_counter[p] else ''
            v.vprint(fmt % (name,str(count),onset,p))

    ordmin, ordmax, ordplotmin, ordplotmax = emu.get_daterange(date_counter,args.dates)
    cum_counts = emu.get_cumulative_counts(date_counter,(ordmin, ordmax),
                                       daysperweek=args.daily)

    ## Only keep data that is within the specified date range
    ## That way, automatic scaling on the y-axis will be based on available data
    ## Are we sure this is really necessary!?
    if args.dates:
        for m in cum_counts:
            ztmp = []
            for i,cnt in enumerate(cum_counts[m]):
                if ordplotmin <= ordmin+i <= ordplotmax:
                    ztmp.append(cnt)
            cum_counts[m] = ztmp
        ordmin = ordplotmin

    emu.make_emberstyle_plots(args,'bynames',cum_counts,T.names,T.colors,ordmin,
                              ordplotrange = (ordplotmin,ordplotmax),
                              title = covid.get_title(args) + ": %d sequences" % (nmatches,),
                              onsets=onsets)


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
