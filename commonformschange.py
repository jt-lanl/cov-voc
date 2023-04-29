'''
Show how the most common sequence forms change over two halves of an interval;
can be used to look for rapidly increasing variants
'''
## note, this is based on, indeed derived from, pangocommonforms
## 1/ it shares a lot of the same code.  ideally, that common code should
##    be abstacted into a module that this and pangocommonforms could import
## 2/ there may be vestigal features of pangocommonforms that are still here,
##    but which don't make sense in the context of how common forms change

from collections import Counter,defaultdict
import datetime
import argparse
import numpy as np
from scipy import stats

import sequtil
import covid
import mutant
import verbose as v
from hamming import hamming

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--npatterns","-n",type=int,default=0,
        help="How many of the most common patterns per lineage (0=all)")
    paa("--mincount","-m",type=int,default=10,
        help="Show only patterns with at least this many counts")
    paa("--consensusalways","-c",action="store_true",
        help="Always show consensus, even if not common")
    paa("--consensusnever",action="store_true",
        help="Do not compute consensus for each form [faster compute]")
    paa("--protein",default="Spike",
        help="Protein name to be used in the header")
    paa("--baseline",default=None,
        choices=tuple(covid.BASELINE_MSTRINGS),
        help="Use this sequence as basline for mutation strings")
    paa("--lineagebaseline",action="store_true",
        help="Use each lineage most common form as mstring baseline for that lineage")
    paa("--nopango",action="store_true",
        help="Don't divide up sequences by pango name")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    if args.consensusalways and args.consensusnever:
        raise RuntimeError("Cannot have both --consensusalways and --consensusnever")
    return args

def mostcommonchar(clist):
    '''return the most common item in the list'''
    [(c,_)] = Counter(clist).most_common(1)
    return c

def consensus(seqlist):
    '''create a consesnsus sequence from the sequence list'''
    return "".join(mostcommonchar(clist)
                   for clist in sequtil.gen_columns_seqlist(seqlist))

def print_header(args):
    '''print the header before the table itself'''
    print(f"COMMON FORMS OF {args.protein.upper()} "
          f"WITH A GIVEN PANGO LINEAGE DESIGNATION")
    print()
    print(f"For each lineage, we show the most common forms of {args.protein}, "
          "as well as their counts within (and percentages of) the lineage.  "
          "[Note that if a lineage contains several divergent forms, "
          "the consensus form might not be found among these common forms.] "
          "Also shown is the Hamming distance (HD) between each form "
          "and the most common form in that lineage. Deletions relative to "
          "the baseline reference strain are indicated with a dash "
          "(e.g. the two amino acid deletion at positions 156-157 is "
          "indicated with 'E156-,F157-'), "
          "and insertions are denoted by a plus sign "
          "(e.g. an extra T at position 143 is written '+143T'). ")
    if args.baseline:
        print(f"[Note: Mutation strings are relative to baseline {args.baseline}].")
    if args.lineagebaseline:
        print("[Note: Mutation strings are relative to the most common variant in each lineage.]")
    print()

    count_forms = f"the {args.npatterns} most common" if args.npatterns \
        else "all the"
    min_count = \
        f" that have at least {args.mincount} counts " \
        "(but we always show the most common form)" \
        if args.mincount>1 else ""
    consensus_always = \
        "And we always show the consensus form; " \
        "this form has the most common amino acid at each position. " \
        if args.consensusalways else ""
    print(f"We show {count_forms} forms{min_count}. {consensus_always}")

def split_date_range(date_range):
    '''Split a date range into early and later halves'''
    ## input date_range tuple can be datetime.date objects or iso-strings
    ## output is two tuples of datetime.date objects
    start_date, end_date = tuple(map(covid.date_fromiso,date_range))
    mid_date = start_date + (end_date - start_date) // 2
    early_range = (start_date, mid_date)
    later_range = (mid_date + datetime.timedelta(days=1), end_date)
    return early_range, later_range

def main(args):
    '''pangocommonforms main'''

    if args.baseline and args.lineagebaseline:
        raise RuntimeError('Cannot have both --baseline and --lineagebaseline')

    print_header(args)

    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqlist = sequtil.get_first_item(seqs,keepfirst=False)
    firstseq = first.seq


    ## baseline mutation for mstrings
    base_mut = mutant.Mutation(covid.get_baseline_mstring(args.baseline))
    if args.baseline:
        print()
        print("Baseline:",str(base_mut))
        print()

    seqlist = list(seqlist)
    v.vprint_only_summary('Invalid date:','skipped sequences')

    last_days = f" in the last {args.days} days from our last update,"
    last_days = last_days if args.days else ""
    try:
        (f_date,t_date) = covid.range_of_dates(seqlist)
    except ValueError:
        (f_date,t_date) = ('Unknown','Unknown')

    early,later = split_date_range((f_date,t_date))
    n_early = len(list(covid.filter_by_date(seqlist,*early)))
    n_later = len(list(covid.filter_by_date(seqlist,*later)))

    print(f"This output is based on sequences sampled{last_days} "
          "from %s to %s." % (f_date,t_date))
    print(f"This interval is split into two intervals.")
    print(f"Early: {early[0].isoformat()} to {early[1].isoformat()} ({n_early})")
    print(f"Later: {later[0].isoformat()} to {later[1].isoformat()} ({n_later})")

    ## Partition seqlist by lineages, separate list for each lineage
    seqlist_by_lineage=defaultdict(list)
    for s in seqlist:
        lin = covid.get_lineage_from_name(s.name) if not args.nopango else "N/A"
        lin = lin or "None"
        seqlist_by_lineage[lin].append(s)

    cnt_lin = {lin: len(seqlist_by_lineage[lin]) for lin in seqlist_by_lineage}
    lineages = sorted(cnt_lin,key=cnt_lin.get,reverse=True)

    ## format lineage strings so they line up
    maxlinlen = max(len(lin) for lin in lineages+["Lineage"])
    fmt = "%%-%ds" % (maxlinlen,)
    fmt_lin = {lin: fmt%lin for lin in lineages}
    fmt_lin["EMPTY"] = fmt % ("",)

    ## print header table
    print()
    print(fmt % "Pango",
          "Lineage    Form   Form    Counts        Fractions       Differences -log10")
    print(fmt % "Lineage",
          "  Count   Count    Pct  Early/Later    Early/Later      Abs Relative pval  HD [Form as mutation string]")

    mut_manager = mutant.MutationManager(firstseq)

    for lin in lineages:

        seqlin = seqlist_by_lineage[lin]

        ## First get consensus form
        cons = consensus(seqlin) if not args.consensusnever else None

        seqlin_early = list(covid.filter_by_date(seqlin,*early))
        seqlin_later = list(covid.filter_by_date(seqlin,*later))
        ne,nl = len(seqlin_early),len(seqlin_later)
        v.vvprint('early:',ne)
        v.vvprint('later:',nl)

        if nl+ne < args.mincount:
            continue

        _,pval = stats.fisher_exact([[ne,nl],[n_early-ne,n_later-nl]])
        print()
        print("%s %7d %7d   100%% %6d/%-6d"
              "                                 "
              "%4.1f     Full lineage" %
              (fmt_lin[lin],ne+nl,ne+nl,ne,nl,-np.log10(pval+1.25e-100)))
        #-np.log10(pval)))

        ## Now get most common forms
        cntr_both = Counter(s.seq for s in seqlin)
        top_comm = sorted(cntr_both,key=cntr_both.get,reverse=True)[0]
        lineage_baseline = mut_manager.get_mutation(top_comm)
        cntr_early = Counter(s.seq for s in seqlin_early)
        cntr_later = Counter(s.seq for s in seqlin_later)
        cflag = False
        def relative_diff(comm):
            ce,cl = cntr_early[comm],cntr_later[comm]
            #den = cl/nl if cl*ne > ce*nl else ce/ne
            return 100*(cl*ne-ce*nl)/max([cl*ne,ce*nl,1])
        def neg_log_pval(comm):
            ce,cl = cntr_early[comm],cntr_later[comm]
            oddrat,pval = stats.fisher_exact([[ce,cl],[n_early-ce,n_later-cl]])
            return -np.log10(pval+1.26e-100) ## caps neglog at 99.9
        cntrlist = sorted(cntr_both,key=neg_log_pval,reverse=True) #key=cntr_later.get
        if args.npatterns:
            cntrlist = cntrlist[:args.npatterns]
        for n,comm in enumerate(cntrlist):
            ce,cl = cntr_early[comm],cntr_later[comm]
            if cl+ce < args.mincount:
                continue
            cons_string = ""
            if comm == cons:
                cflag = True
                cons_string = "(consensus)"
            m = mut_manager.get_mutation(comm)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)

            if args.lineagebaseline:
                mstring = m.relative_to(lineage_baseline)

            h = hamming(top_comm,comm)
            cene = ce/ne if ne>0 else np.inf
            clnl = cl/nl if nl>0 else np.inf
            oddrat,pval = stats.fisher_exact([[ce,cl],[n_early-ce,n_later-cl]])
            print("%s %7d %7d %5.1f%% %6d/%-6d %7.5f/%7.5f %+5.2f%% %+7.1f%% %4.1f %3d %s %s" %
                  (fmt_lin[lin],cnt_lin[lin],ce+cl,
                   100*(ce+cl)/cnt_lin[lin],
                   ce,cl,
                   cene,clnl,
                   100*(clnl-cene),
                   relative_diff(comm),
                   neg_log_pval(comm),
                   h,mstring,cons_string))

        ## print statement is out of date, so use if 0 for now
        if 0 and args.consensusalways and not cflag:
            h = hamming(top_comm,cons)
            m = mut_manager.get_mutation(cons)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)
            cnt = cntr[cons]
            print("%s %7d %6d/%6d %5.1f%% %3d %s %s" %
                  (fmt_lin[lin],cnt_lin[lin],ce,cl,
                       100*(ce+cl)/cnt_lin[lin],
                   h,mstring,"(consensus)"))

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
