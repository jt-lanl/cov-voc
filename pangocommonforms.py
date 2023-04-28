'''
From sequence file, obtain all the pango lineages, and for each one,
find the mutation string that corresponds to the consensus sequence
for that lineage, and show the most commmon forms.
'''
## note, consensus is the most expensive part of the computation
## use --consensusnever to avoid that computation

from collections import Counter,defaultdict
import argparse

import sequtil
import covid
import mutant
from verbose import verbose as v
from hamming import hamming

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
    paa("--baseline",default=None,
        choices=tuple(covid.BASELINE_MSTRINGS),
        help="Use this sequence as basline for mutation strings")
    paa("--lineagebaseline",action="store_true",
        help="Use each lineage most common form as mstring baseline for that lineage")
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

def main(args):
    '''pangocommonforms main'''

    if args.baseline and args.lineagebaseline:
        raise RuntimeError('Cannot have both --baseline and --lineagebaseline')

    print_header(args)

    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqlist = sequtil.get_first_item(seqs,keepfirst=False)
    firstseq = first.seq
    
    seqlist = list(seqlist)
    v.vprint_only_summary('Invalid date:','skipped sequences')

    last_days = f" in the last {args.days} days from our last update,"
    last_days = last_days if args.days else ""
    try:
        (f_date,t_date) = covid.range_of_dates(seqlist)
    except ValueError:
        (f_date,t_date) = ('Unknown','Unknown')
        
    print(f"This output is based on sequences sampled{last_days} "
          "from %s to %s." % (f_date,t_date))

    seqlist_by_lineage=defaultdict(list)
    for s in seqlist:
        lin = covid.get_lineage_from_name(s.name)
        lin = lin or "None"
        seqlist_by_lineage[lin].append(s)

    cnt_lin = {lin: len(seqlist_by_lineage[lin]) for lin in seqlist_by_lineage}
    lineages = sorted(cnt_lin,key=cnt_lin.get,reverse=True)

    maxlinlen = max(len(lin) for lin in lineages)
    fmt = "%%-%ds" % (maxlinlen,)
    fmt_lin = {lin: fmt%lin for lin in lineages}
    fmt_lin["EMPTY"] = fmt % ("",)

    mut_manager = mutant.MutationManager(firstseq)

    ## Get baseline:
    ##     For Spike, use hardcoded baseline sequences
    ##     For Other proteins, use most common form of the given baseline pango type
    if not args.baseline:
        base_mut = mutant.Mutation("[]")
    elif args.protein == 'Spike':
        base_mut = mutant.Mutation(covid.get_baseline_mstring(args.baseline))
    else:
        v.vprint('Will obtain baseline from most common',args.baseline)
        if args.baseline not in seqlist_by_lineage:
            ## should this be fatal?
            v.vprint(f'Baseline {args.baseline} not in data!')
        else:
            cntr = Counter(s.seq for s in seqlist_by_lineage[args.baseline])
            base_seq = cntr.most_common(1)[0][0]
            base_mut = mut_manager.get_mutation(base_seq)
    if args.baseline:
        print()
        print("Baseline:",str(base_mut))
        print()

    print()
    print(fmt % "Pango","Lineage   Form   Form")
    print(fmt % "Lineage","  Count  Count    Pct  HD [Form as mutation string]")

    for lin in lineages:

        #if lin == "None": continue
        #if not lin: continue

        seqlin = seqlist_by_lineage[lin]

        ## First get consensus form
        cons = consensus(seqlin) if not args.consensusnever else None

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
                h = 0
                print("%s %7d %6d %5.1f%% %3d %s %s" %
                      (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],
                       h,mstring,cons_string))
            else:
                h = hamming(top_comm,comm)
                print("%s %7d %6d %5.1f%% %3d %s %s" %
                      (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],
                       h,mstring,cons_string))
        if args.consensusalways and not cflag:
            h = hamming(top_comm,cons)
            m = mut_manager.get_mutation(cons)
            mstring = m.relative_to(base_mut) if args.baseline else str(m)
            cnt = cntr[cons]
            print("%s %7d %6d %5.1f%% %3d %s %s" %
                  (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],
                   h,mstring,"(consensus)"))
        print()

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
