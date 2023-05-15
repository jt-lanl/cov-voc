'''
Show how the most common sequence forms change over two halves of an interval;
can be used to look for rapidly increasing variants
'''
## note, this is based on, indeed derived from, pangocommonforms
## 1/ it shares a lot of the same code.  ideally, that common code should
##    be abstacted into a module that this and pangocommonforms could import
## 2/ there may be vestigal features of pangocommonforms that are still here,
##    but which don't make sense in the context of how common forms change

import re
from collections import Counter
import datetime
import argparse
import numpy as np
from scipy import stats

import verbose as v
from hamming import hamming

import covid
import mutant
import commonforms as cf

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--npatterns","-n",type=int,default=0,
        help="How many of the most common patterns per lineage (0=all)")
    paa("--mincount","-m",type=int,default=10,
        help="Show only patterns with at least this many counts")
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
    return args

def print_header(args):
    '''print the header before the table itself'''
    print(f"COMMON FORMS CHANGES FOR {args.protein.upper()} "
          f"WITH A GIVEN PANGO LINEAGE DESIGNATION")
    print()
    print(f"For each lineage, we show the most common forms of {args.protein}, "
          "and we show the forms that are most significantly increasing "
          "or decreasing. "
          "Also shown is the Hamming distance (HD) between each form "
          "and the most common form in that lineage. Deletions relative to "
          "the baseline reference strain are indicated with a dash "
          "(e.g. the two amino acid deletion at positions 156-157 is "
          "indicated with 'E156-,F157-'), "
          "and insertions are denoted by a plus sign "
          "(e.g. an extra T at position 143 is written '+143T'). ")
    print(f"Abs Differences in later vs early fractions, expressed as percent.")
    print(f"Relative Differences range between -100% and +100%, using formula: "
          "rel = 100*(later-early)/max(later,early) percent.")
    print(f"The p-value associated with increase or decrease is meant to be used only as a rough guide to the siginficcance of the change. Since the computation is based on assumptions (such as independent and unbiased sampling) that may not hold in practice, the p-value should not be taken too literally.")
    print(f"The first line of each lineage section indicates counts "
          f"for the full lineage relative to all the sequences. "
          f"'Lineage Count' in this first line is actually the full sequence count. "
          f"'Form Count' in this first line is actually the full lineage count, "
          f"and 'Form Pct' refers to lineage count relative to full sequence set.")
    print()

    count_forms = f"the {args.npatterns} most common" if args.npatterns \
        else "all the"
    min_count = \
        f" that have at least {args.mincount} counts " \
        "(but we always show the most common form)" \
        if args.mincount>1 else ""
    print(f"We show {count_forms} forms{min_count}. ")
    print()
    if args.baseline:
        print(f"[Note: Mutation strings are relative to baseline {args.baseline}].")
    if args.lineagebaseline:
        print("[Note: Mutation strings are relative to the most common variant in each lineage.]")

def split_date_range(date_range):
    '''Split a date range into early and later halves'''
    ## input date_range tuple can be datetime.date objects or iso-strings
    ## output is two tuples of datetime.date objects
    start_date, end_date = tuple(map(covid.date_fromiso,date_range))
    mid_date = start_date + (end_date - start_date) // 2
    early_range = (start_date, mid_date)
    later_range = (mid_date + datetime.timedelta(days=1), end_date)
    return early_range, later_range

def split_date_range_bycounts(seqlist):
    '''Produce two adjacent date ranges, one early and one later,
       that cover the full range of dates in the input sequence list,
       with the split chosen so that there is a
       roughly equal number of sequences in each range
    '''
    datelist = [covid.date_from_seqname(s.name) for s in seqlist]
    datelist = [d for d in datelist if d is not None]
    datelist = sorted(datelist)
    start_date = datelist[0]
    mid_date = datelist[len(datelist)//2]
    end_date = datelist[-1]
    early_range = (start_date, mid_date + datetime.timedelta(days=-1))
    later_range = (mid_date, end_date)
    return early_range, later_range

PVALMIN = 1e-9
def strpval(pval):
    lessmin = "<" if pval < PVALMIN else " "
    pval = max([pval,PVALMIN])
    pval = "%5.0e" % (pval,)
    mantissa,exponent = pval.split('e')
    pval = '%1de%+1d' % (int(mantissa),int(exponent))
    return lessmin + pval

def main(args):
    '''commonformschange main'''

    if args.baseline and args.lineagebaseline:
        raise RuntimeError('Cannot have both --baseline and --lineagebaseline')
    if args.nopango and args.lineagebaseline:
        v.print('Warning: --nopango and --lineagebaseline not recommended.')

    print_header(args)

    firstseq,seqlist = cf.get_input_sequences(args)
    mut_manager = mutant.MutationManager(firstseq)

    ## baseline mutation for mstrings (assumes protein==Spike)
    base_mut = mutant.Mutation(covid.get_baseline_mstring(args.baseline))
    if args.baseline:
        print("Baseline:",str(base_mut))
        print()

    early,later = split_date_range_bycounts(seqlist)
    n_early = len(list(covid.filter_by_date(seqlist,*early)))
    n_later = len(list(covid.filter_by_date(seqlist,*later)))

    last_days = f" in the last {args.days} days from our last update,"
    last_days = last_days if args.days else ""
    (f_date,_),(_,t_date) = early,later
    print(f"This output is based on {n_early+n_later} sequences sampled{last_days} "
          f"from {f_date} to {t_date}")
    print(f"This total interval is split into two sub-intervals.")
    print(f"Total: {early[0]} to {later[1]} ({n_early+n_later})")
    print(f"Early: {early[0]} to {early[1]} ({n_early})")
    print(f"Later: {later[0]} to {later[1]} ({n_later})")

    ## Partition seqlist by lineages, separate list for each lineage
    lp = cf.LineagePartition(seqlist,nopango=args.nopango)

    ## print header for table
    print()
    print(lp.format("Pango"),
          "Lineage    Form   Form    Counts        Fractions     "
          " Differences")
    print(lp.format("Lineage"),
          "  Count   Count    Pct  Early/Later    Early/Later    "
          "Abs   Relative  pval   HD [Form as mutation string]")

    table_format = \
        "%s %7d %7d %5.1f%% %6d/%-6d %7.5f/%7.5f "\
        "%+6.2f%% %+6.1f%% %s %3d %s %s"

    for lin in lp.lineages:

        seqlin = lp.sequences[lin]
        countlin = lp.counts[lin]
        fmtlin = lp.format(lin)

        ## First get consensus form
        cons = cf.consensus(seqlin)

        ## Partition sequences into early and later
        seqlin_early = list(covid.filter_by_date(seqlin,*early))
        seqlin_later = list(covid.filter_by_date(seqlin,*later))
        ne,nl = len(seqlin_early),len(seqlin_later)
        v.vvprint('early:',ne)
        v.vvprint('later:',nl)

        if nl+ne < args.mincount:
            continue

        ## Compute pval for the lineage
        if lin != "N/A":
            _,pval = stats.fisher_exact([[ne,nl],[n_early-ne,n_later-nl]])
            print()
            table_line = table_format % \
                (lp.format(''),n_early+n_later,ne+nl,
                 100*(ne+nl)/(n_early+n_later),
                 ne,nl,ne/n_early,nl/n_later,
                 100*(nl/n_later-ne/n_early),
                 100*(nl*n_early-ne*n_later)/max([nl*n_early,ne*n_later]),
                 strpval(pval),
                 0,f" Full {lin}","lineage")
            table_line = re.sub(" ","_",table_line)
            print(table_line)


        ## Now get most common forms
        cntr_both = Counter(s.seq for s in seqlin)
        top_comm = sorted(cntr_both,key=cntr_both.get,reverse=True)[0]
        lineage_baseline = mut_manager.get_mutation(top_comm)
        cntr_early = Counter(s.seq for s in seqlin_early)
        cntr_later = Counter(s.seq for s in seqlin_later)
        def relative_diff(comm):
            ce,cl = cntr_early[comm],cntr_later[comm]
            #den = cl/nl if cl*ne > ce*nl else ce/ne
            return 100*(cl*ne-ce*nl)/max([cl*ne,ce*nl,1])
        def neg_log_pval(comm,cap=False):
            ce,cl = cntr_early[comm],cntr_later[comm]
            _,pval = stats.fisher_exact([[ce,cl],[n_early-ce,n_later-cl]])
            if cap:
                pval = max([pval,PVALMIN])
            return np.log10(1/pval)
        def form_count(comm):
            ce,cl = cntr_early[comm],cntr_later[comm]
            return ce+cl
        cntrlist = sorted(cntr_both,key=form_count,reverse=True)

        if args.npatterns:
            cntrlist = cntrlist[:args.npatterns]

        for comm in cntrlist:
            ce,cl = cntr_early[comm],cntr_later[comm]
            if cl+ce < args.mincount:
                continue
            cons_string = ""
            if comm == cons:
                cons_string = "(consensus)"
            m = mut_manager.get_mutation(comm)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)

            if args.lineagebaseline:
                mstring = m.relative_to(lineage_baseline)

            h = hamming(top_comm,comm)
            cene = ce/ne if ne>0 else np.inf
            clnl = cl/nl if nl>0 else np.inf
            _,pval = stats.fisher_exact([[ce,cl],[n_early-ce,n_later-cl]])
            print(table_format %
                  (fmtlin,countlin,ce+cl,
                   100*(ce+cl)/countlin,
                   ce,cl,
                   cene,clnl,
                   100*(clnl-cene),
                   relative_diff(comm),
                   strpval(pval),
                   h,mstring,cons_string))

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
