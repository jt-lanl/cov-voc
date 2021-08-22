'''
From sequence file, obtain all the pango lineages, and for each one,
find the mutation string that corresponds to the consensus sequence
for that lineage, and show the most commmon forms.
'''
import sys
import re
from collections import Counter,defaultdict
import argparse

import sequtil
import covid
import mutanty as mutant
import wrapgen
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
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def get_matches(mutantstring,seqlist,fullmatch=False):
    '''return sequences in seqlist that match the mutantstring'''

    firstseq = seqlist[0].seq
    muts = mutant.Mutation().init_from_line(mutantstring)
    assert muts.checkref(firstseq,verbose=True)
    matchpatt = muts.regex_pattern(firstseq,exact=bool(fullmatch))
    rematchpatt = re.compile(matchpatt)

    seqmatches = [s for s in seqlist[1:]
                  if rematchpatt.match(s.seq)]

    return seqmatches

def get_lineage_from_name(name):
    return re.sub(r".*EPI_ISL_\d+\.","",name)

def count_lineages(seqlist):
    '''return a counter with number of sequences for each lineage'''
    lineages = [get_lineage_from_name(s.name) for s in seqlist]
    return Counter(lineages)

def format_counter(cnt):
    '''make a nicely formatted string summarizing contents of counter'''
    total = sum(cnt.values())
    scnt = sorted(cnt,key=cnt.get,reverse=True)
    s="%6d %3d " % (total,len(cnt))
    for lin in scnt[:3]:
        s += "%6d %5.1f%% %10s; " % (cnt[lin],100*cnt[lin]/total,lin)
    return s

def mostcommonchar(clist):
    '''return the most common item in the list'''
    [(c,_)] = Counter(clist).most_common(1)
    return c

def consensus(seqlist):
    '''create a consesnsus sequence from the sequence list'''
    return "".join(mostcommonchar(clist)
                   for clist in sequtil.gen_columns_seqlist(seqlist))

def main(args):

    print("COMMON FORMS OF SPIKE WITH A GIVEN PANGO LINEAGE DESIGNATION")
    print()
    print("For each lineage, we show the most common forms of Spike, "
          "as well as their counts within (and percentages of) the lineage.  "
          "[Note that if a lineage contains several divergent forms of Spike, "
          "the consensus form might not be found among these common forms.] "
          "Also shown is the Hamming distance (HD) between each form "
          "and the most common form in that lineage. Deletions relative to "
          "the ancestral reference strain are indicated with a dash "
          "(e.g. the two amino acid deletion in Delta variant are indicated as: "
          "'E156-,F157-'), and insertions are denoted by a plus sign "
          "(e.g. B.1.621 usually has an extra T at position 143, "
          "it is written '+143T')."
          "")
    print()

    count_forms = f" {args.npatterns}" if args.npatterns else ""
    min_count = f" that have at least {args.mincount} counts (but we always show the most common form)" if args.mincount else ""
    consensus_always = f"And we always show the consensus form. " if args.consensusalways else ""
    print(f"We show the{count_forms} most common forms{min_count}. {consensus_always}")

    seqs = covid.read_filter_seqfile(args)
    seqlist = list(seqs)

    last_days = f" in the last {args.days} days from our last update," if args.days else ""
    print(f"This output is based on sequences sampled{last_days} from %s to %s." \
          % sequtil.range_of_dates(seqlist))

    firstseq = seqlist[0].seq

    cnt_lin = count_lineages(seqlist[1:])
    lineages = sorted(cnt_lin,key=cnt_lin.get,reverse=True)

    seqlist_by_lineage=defaultdict(list)
    for s in seqlist[1:]:
        lin = get_lineage_from_name(s.name)
        seqlist_by_lineage[lin].append(s)

    maxlinlen = max(len(lin) for lin in lineages)
    fmt = "%%-%ds" % (maxlinlen,)
    fmt_lin = {lin: fmt%lin for lin in lineages}
    fmt_lin["EMPTY"] = fmt % ("",)

    print()
    print(fmt % "Pango","Lineage   Form   Form")
    print(fmt % "Lineage","  Count  Count    Pct  HD [Spike form as mutation string]")

    mut_manager = mutant.MutationManager(firstseq)

    for lin in lineages:

        if lin == "None":
            continue
        if not lin:
            continue
        
        seqlin = seqlist_by_lineage[lin]

        ## First get consensus form
        cons = consensus(seqlin)
        cntcons = sum(1 for s in seqlin if s.seq == cons)
        mcons = mut_manager.get_mutation(cons) #mutant.Mutation((firstseq,cons))

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
            m = mut_manager.get_mutation(comm) #ant.Mutation((firstseq,comm))
            if n == 0:
                top_comm = comm
                h = 0
                print("%s %7d %6d %5.1f%% %3d %s %s" %
                      (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],
                       h,str(m),cons_string))
            else:
                h = hamming(top_comm,comm)
                print("%s %7d %6d %5.1f%% %3d %s %s" %
                      (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],
                       h,str(m),cons_string))
        if args.consensusalways and not cflag:
            h = hamming(top_comm,cons)
            print("%s %7d %6d %5.1f%% %3d %s (consensus)" %
                  (fmt_lin[lin],cnt_lin[lin],cntcons,100*cntcons/cnt_lin[lin],
                   h,str(mcons)))
        print()

if __name__ == "__main__":

    _args = getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    def vcount(seqs,*p,**kw):
        '''if verbose, count seqs as they go by'''
        return wrapgen.keepcount(seqs,*p,**kw) if _args.verbose else seqs

    main(_args)
