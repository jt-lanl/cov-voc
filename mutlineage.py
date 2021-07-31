DESCRIPTION='''
read mutant strings and linege assignments from a list
find sequences in a fasta file that match the mutant string 
and tabulate the lineages associated with them
'''
import sys
import re
from pathlib import Path
from collections import Counter
import argparse

import warnings

import readseq
import sequtil
import intlist
import covid
import mutant

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--mutant","-m",
        help="file with list of mutant strings and lineages")
    paa("--output","-o",type=Path,
        help="output fasta file")
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

    seqmatches = (s for s in seqlist[1:]
                  if rematchpatt.match(s.seq))

    return seqmatches

def get_lineage_from_name(name):
    return re.sub(".*EPI_ISL_\d+\.","",name)

def count_lineages(seqlist):
    lineages = (get_lineage_from_name(s.name) for s in seqlist)
    return Counter(lineages)

def format_counter(cnt):
    '''make a nicely formatted string summarizing contents of counter'''
    total = sum(cnt.values())
    scnt = sorted(cnt,key=cnt.get,reverse=True)
    s="%6d %3d " % (total,len(cnt))
    for lin in scnt[:3]:
        s += "%6d %5.1f%% %10s; " % (cnt[lin],100*cnt[lin]/total,lin)
    return s

def main(args):

    mutantfilepatt = re.compile(r".*(\[.*\])\s+(\(.*\))")
    mutlist=[]
    linlist=[]
    with open(args.mutant) as f:
        for line in f:
            line = line.strip()
            m = mutantfilepatt.match(line)
            if m:
                mstring,lineage = m[1],m[2]
                vprint(mstring,lineage)
                mutlist.append(mstring)
                linlist.append(lineage)
            else:
                vprint("Invalid line:",line)

    mutmaxlen = max(len(mut) for mut in mutlist)
    linmaxlen = max(len(lin) for lin in linlist)
    fmt = "%%%ds %%-%ds" % (mutmaxlen,linmaxlen)
    for mut,lin in zip(mutlist,linlist):
        vprint(fmt % (mut,lin),)

    seqs = covid.read_filter_seqfile(args)
    seqlist = list(seqs)

    for fullmatch in [False,True]:
        print("Exact matches:" if fullmatch else "Inclusive matches:")
    
        for mut,lin in zip(mutlist,linlist):
            matches = get_matches(mut,seqlist,fullmatch=fullmatch)
            cnt = count_lineages(matches)
            print(fmt % (mut,lin),format_counter(cnt))
        

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)


    main(args)
    

