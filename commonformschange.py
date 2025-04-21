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
    cf.commonforms_args(ap)
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    ap.set_defaults(bylineage=False)
    args = ap.parse_args()
    args = cf.commonforms_fixargs(args)
    return args

def print_header(args):
    '''print the header before the table itself'''
    print(f"COMMON FORMS CHANGES FOR {args.protein.upper()} "
          "WITH A GIVEN PANGO LINEAGE DESIGNATION")
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
    print("Abs Differences in later vs early fractions, expressed as percent.")
    print("Relative Differences range between -100% and +100%, using formula: "
          "rel = 100*(later-early)/max(later,early) percent, applied to the "
          "early and later fractions.")
    print("The p-value associated with increase or decrease is computed by "
          "Fisher's exact test, and is meant to be used only as a rough guide "
          "to the significance of the change. Since the computation is based on "
          "assumptions (such as independent and unbiased sampling) that may not "
          "hold in practice, the p-value should not be taken too literally.")
    print("The first line of each lineage section indicates counts "
          "for the full lineage relative to all the sequences. "
          "'Lineage Count' in this first line is actually the full sequence count. "
          "'Form Count' in this first line is actually the full lineage count, "
          "and 'Form Pct' refers to lineage count relative to full sequence set.")
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
    datelist = [covid.get_date(s.name) for s in seqlist]
    datelist = [d for d in datelist if d is not None]
    datelist = sorted(datelist)
    start_date = datelist[0]
    mid_date = datelist[len(datelist)//2]
    end_date = datelist[-1]
    early_range = (start_date, mid_date + datetime.timedelta(days=-1))
    later_range = (mid_date, end_date)
    return early_range, later_range

def relative_diff(ne,nl,ce,cl):
    '''informally: relative difference between ce/ne and cl/nl;
    but bounded in range (-100,100) and ok when ne=0 or nl=0'''
    return 100*(cl*ne-ce*nl)/max([cl*ne,ce*nl,1])

PVALMIN = 1.2589e-10 ## caps neglog at 9.9
def neg_log_pval(pval):
    '''compute neg log of pvalue, capping value at 9.9'''
    pval = max([pval,PVALMIN])
    return np.log10(1/pval)

def strpval(pval):
    '''convert p-value into a string of the form, eg " 2e-4", or else "<1e-9"
    '''
    lessmin = "<" if pval < PVALMIN else " "
    pval = max([pval,PVALMIN])
    pval = "%5.0e" % (pval,)
    mantissa,exponent = pval.split('e')
    pval = '%1de%+1d' % (int(mantissa),int(exponent))
    return lessmin + pval

def main(args):
    '''commonformschange main'''

    firstseq,seqlist = cf.get_input_sequences(args,minseqs=2)
    mut_manager = mutant.MutationManager(firstseq)

    print_header(args)

    early,later = split_date_range_bycounts(seqlist)
    total_ne = len(list(covid.filter_by_date(seqlist,*early)))
    total_nl = len(list(covid.filter_by_date(seqlist,*later)))

    last_days = f" in the last {args.days} days from our last update,"
    last_days = last_days if args.days else ""
    (f_date,_),(_,t_date) = early,later
    print(f"This output is based on {total_ne+total_nl} sequences sampled{last_days} "
          f"from {f_date} to {t_date}")
    print("This total interval is split into two sub-intervals.")
    print(f"Total: {early[0]} to {later[1]} ({total_ne+total_nl})")
    print(f"Early: {early[0]} to {early[1]} ({total_ne})")
    print(f"Later: {later[0]} to {later[1]} ({total_nl})")

    ## Partition seqlist by lineages, separate list for each lineage
    lp = cf.LineagePartition(seqlist)
    base_mut = cf.get_baseline_mutation(args.baseline,mut_manager,lp,args.protein)
    if not args.bylineage:
        ## Re-partition seqlist into one big partition, not by lineage
        lp = cf.LineagePartition(seqlist,bylineage=False)

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
            _,pval = stats.fisher_exact([[ne,nl],[total_ne-ne,total_nl-nl]])
            print()
            rne = ne/total_ne if total_ne>0 else 0
            rnl = nl/total_nl if total_nl>0 else 0
            if max([nl*total_ne,ne*total_nl]) == 0:
                vprint("ZERO:",nl,total_ne,ne,total_nl)
            table_line = table_format % \
                (lp.format(''),total_ne+total_nl,ne+nl,
                 100*(ne+nl)/(total_ne+total_nl),
                 ne,nl,rne,rnl,
                 100*(rnl-rne),
                 100*(nl*total_ne-ne*total_nl)/max([nl*total_ne,ne*total_nl,1]),
                 strpval(pval),
                 0,f" Full {lin}","lineage")
            table_line = re.sub(" ","_",table_line)
            print(table_line)


        ## Now get most common forms
        cntr_both = Counter(s.seq for s in seqlin)
        top_comm = sorted(cntr_both,key=cntr_both.get,reverse=True)[0]
        lineage_baseline = mut_manager.seq_to_mutation(top_comm)
        cntr_early = Counter(s.seq for s in seqlin_early)
        cntr_later = Counter(s.seq for s in seqlin_later)
        cntrlist = sorted(cntr_both,key=cntr_both.get,reverse=True)

        if args.npatterns:
            cntrlist = cntrlist[:args.npatterns]

        for comm in cntrlist:
            ce,cl = cntr_early[comm],cntr_later[comm]
            if cl+ce < args.mincount:
                continue
            cons_string = ""
            if comm == cons:
                cons_string = "(consensus)"
            mut = mut_manager.seq_to_mutation(comm)
            mstring = mut.relative_to(base_mut) if args.baseline else str(mut)

            if args.lineagebaseline:
                mstring = mut.relative_to(lineage_baseline)

            hdist = hamming(top_comm,comm)
            cene = ce/ne if ne>0 else (np.inf if ne>0 else 0)
            clnl = cl/nl if nl>0 else (np.inf if nl>0 else 0)
            _,pval = stats.fisher_exact([[ce,cl],[total_ne-ce,total_nl-cl]])
            print(table_format %
                  (fmtlin,countlin,ce+cl,
                   100*(ce+cl)/countlin,
                   ce,cl,
                   cene,clnl,
                   100*(clnl-cene),
                   relative_diff(ne,nl,ce,cl),
                   strpval(pval),
                   hdist,mstring,cons_string))

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)
    main(_args)
