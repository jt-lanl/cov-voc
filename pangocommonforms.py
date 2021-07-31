DESCRIPTION='''
from sequence file, obtain all the pango lineages, and for each one,
find the mutation string that corresponds to the consensus sequence
for that lineage, and show the most commmon forms
'''
import sys
import re
from pathlib import Path
from collections import Counter,defaultdict
import argparse

import warnings

import readseq
import sequtil
import intlist
import covid
import mutant
import wrapgen



def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--mostcommon","-m",type=int,default=5,
        help="How many of the most common patterns per lineage")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def get_matches(mutantstring,seqlist,fullmatch=False):

    firstseq = seqlist[0].seq
    
    muts = mutant.Mutation().init_from_line(mutantstring)
    assert( muts.checkref(firstseq,verbose=True) )
    matchpatt = muts.regex_pattern(firstseq,exact=bool(fullmatch))
    rematchpatt = re.compile(matchpatt)

    seqmatches = [s for s in seqlist[1:]
                  if rematchpatt.match(s.seq)]

    return seqmatches

def get_lineage_from_name(name):
    return re.sub(".*EPI_ISL_\d+\.","",name)

def count_lineages(seqlist):
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
    [(c,cnt)] = Counter(clist).most_common(1)
    return c

def consensus(seqlist):
    return "".join(mostcommonchar(clist) for clist in sequtil.gen_columns_seqlist(seqlist))

def main(args):

    seqs = covid.read_filter_seqfile(args)
    seqlist = list(seqs)

    print("COMMON FORMS OF SPIKE WITH A GIVEN PANGO LINEAGE DESIGNATION")
    print()
    print("For each lineage, we show the consensus form of Spike, "
          "as well as its count within (and percentage of) the lineage.  "
          f"We also show the {args.mostcommon} most common forms, "
          "which may or may not include the consensus form.  "
          "[Note that if a lineage contains several divergent forms of Spike, "
          "the consensus form might not be found in nature.]  "
          "Also, note that insertions are not included in these patterns.")
    print()
    
    lastDays = f" in the last {args.days} days from our last update," if args.days else ""
    print(f"This output is based on sequences sampled{lastDays} from %s to %s." \
          % sequtil.range_of_dates(seqlist))

    firstseq = seqlist[0].seq

    cnt_lin = count_lineages(seqlist[1:])
    lineages = sorted(cnt_lin,key=cnt_lin.get,reverse=True)

    list_by_lineage=defaultdict(list)
    for s in seqlist[1:]:
        lin = get_lineage_from_name(s.name)
        list_by_lineage[lin].append(s)

    maxlinlen = max(len(lin) for lin in lineages)
    fmt = "%%-%ds" % (maxlinlen,)
    fmt_lin = {lin: fmt%lin for lin in lineages}
    fmt_lin["EMPTY"] = fmt % ("",)

    print()
    print(fmt % "Pango","Lineage   Form   Form")
    print(fmt % "Lineage","  Count  Count    Pct [Mutation string]")
        
    for lin in lineages:
        seqlin = list_by_lineage[lin]
        
        cons = consensus(seqlin)
        cnt = sum(1 for s in seqlin if s.seq == cons)
        m = mutant.Mutation((firstseq,cons))
        print("%s %7d %6d %5.1f%% %s (consensus)" % (fmt_lin[lin],cnt_lin[lin],cnt,100*cnt/cnt_lin[lin],str(m)))

        cntr = Counter(s.seq for s in seqlin)
        for comm,cnt in cntr.most_common(args.mostcommon):
            if comm == cons:
                continue
            m = mutant.Mutation((firstseq,comm))
            print("%s         %6d %5.1f%% %s" % (fmt_lin["EMPTY"],cnt,100*cnt/cnt_lin[lin],str(m)))
        
              
        

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    def vcount(seqs,*p,**kw):
        if args.verbose:
            return wrapgen.keepcount(seqs,*p,**kw)
        else:
            return seqs
        
    main(args)
    

