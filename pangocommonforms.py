'''
From sequence file, obtain all the pango lineages, and for each one,
find the mutation string that corresponds to the consensus sequence
for that lineage, and show the most commmon forms.
'''
## note, consensus is the most expensive part of the computation
## use --consensusnever to avoid that computation

from collections import Counter
import argparse

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
        help="Show only patters with at least this many counts")
    paa("--consensusalways","-c",action="store_true",
        help="Always show consensus, even if not common")
    paa("--consensusnever",action="store_true",
        help="Do not compute consensus for each form [faster compute]")
    paa("--protein",default="Spike",
        help="Protein name to be used in the header")
    paa("--baseline",default="XBB.1.5",
        choices=tuple(covid.BASELINE_MSTRINGS),
        help="Use this sequence as basline for mutation strings")
    paa("--bylineage",action="store_true",default=True,
        help="Partition sequences by pango lineage")
    paa("--notbylineage",action="store_false",dest='bylineage',
        help="Do not partition sequences by pango lineges")
    paa("--lineagebaseline",action="store_true",
        help="Use each lineage most common form as mstring baseline for that lineage")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    if args.consensusalways and args.consensusnever:
        raise RuntimeError("Cannot have both --consensusalways and --consensusnever")
    if args.baseline == 'Wuhan':
        args.baseline = None
    return args

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

def main(args):
    '''pangocommonforms main'''

    if args.baseline and args.lineagebaseline:
        raise RuntimeError('Cannot have both --baseline and --lineagebaseline')
    if not args.bylineage and args.lineagebaseline:
        v.print('Warning: not recommended to use --notbylineage AND --lineagebaseline.')

    firstseq,seqlist = cf.get_input_sequences(args)
    mut_manager = mutant.MutationManager(firstseq)

    if len(seqlist) < 1:
        v.print(args)
        v.print(f'Only {len(seqlist)} sequences -- aborting.')
        return

    print_header(args)

    last_days = f" in the last {args.days} days from our last update,"
    last_days = last_days if args.days else ""
    try:
        (f_date,t_date) = covid.range_of_dates(seqlist)
    except ValueError:
        (f_date,t_date) = ('Unknown','Unknown')

    print(f"This output is based on sequences sampled{last_days} "
          "from %s to %s." % (f_date,t_date))

    ## Partition seqlist by lineages, separate list for each lineage
    lin_partition = cf.LineagePartition(seqlist)

    base_mut = cf.get_baseline_mutation(args.baseline,mut_manager,
                                        lin_partition,args.protein)

    if not args.bylineage:
        lin_partition = cf.LineagePartition(seqlist,bylineage=False)

    ## Print header for table:
    print()
    print(lin_partition.format("Pango"),
          "Lineage   Form   Form")
    print(lin_partition.format("Lineage"),
          "  Count  Count    Pct  HD [Form as mutation string]")

    for lin in lin_partition.lineages:

        seqlin = lin_partition.sequences[lin]
        countlin = lin_partition.counts[lin]
        fmtlin = lin_partition.format(lin)

        print()

        ## First get consensus form
        cons = cf.consensus(seqlin) if not args.consensusnever else None

        ## Now get most common forms
        cntr = Counter(s.seq for s in seqlin)
        cflag = False
        cntrlist = sorted(cntr,key=cntr.get,reverse=True)
        if args.npatterns:
            cntrlist = cntrlist[:args.npatterns]
        for n,comm in enumerate(cntrlist):
            cnt = cntr[comm]
            if n>0 and cnt < args.mincount:
                break
            cons_string = ""
            if comm == cons:
                cflag = True
                cons_string = "(consensus)"
            m = mut_manager.get_mutation(comm)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)

            if args.lineagebaseline:
                if n==0:
                    lineage_baseline = m
                else:
                    mstring = m.relative_to(lineage_baseline)

            if n == 0:
                top_comm = comm
            h = hamming(top_comm,comm)
            print("%s %7d %6d %5.1f%% %3d %s %s" %
                  (fmtlin,countlin,cnt,100*cnt/countlin,
                   h,mstring,cons_string))
        if args.consensusalways and not cflag:
            h = hamming(top_comm,cons)
            m = mut_manager.get_mutation(cons)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)
            cnt = cntr[cons]
            print("%s %7d %6d %5.1f%% %3d %s %s" %
                  (fmtlin,countlin,cnt,100*cnt/countlin,
                   h,mstring,"(consensus)"))

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
