'''
read mutant strings and linege assignments from a list
each line format: [MutantString] (Lineage)
find sequences in a fasta file that match the mutant string
and tabulate the lineages associated with them
'''
import sys
import re
from collections import Counter
import itertools as it
import argparse

import sequtil
import covid
import mutant

def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--mutant","-m",
        help="file with list of mutant strings and lineages")
    paa("-N",type=int,default=0,
        help="just look at the first N sequences in input file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def get_lineage_from_name(name):
    '''return lineage (eg, B.1.1.7) from full sequence name'''
    return re.sub(r".*EPI_ISL_\d+\.","",name)

def count_lineages(seqlist):
    '''make a counter of all the lineages in the squlist'''
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

def read_mutantfile(filename):
    '''read mutant file and return mstring list and lineage list'''
    mutantfilepatt = re.compile(r".*(\[.*\])\s+(\(.*\)?)")
    mutlist=[]
    linlist=[]
    with open(filename) as fptr:
        for line in fptr:
            line = line.strip()
            m = mutantfilepatt.match(line)
            if m:
                mstring,lineage = m[1],m[2]
                vvprint(mstring,lineage)
                mutlist.append(mstring)
                linlist.append(lineage)
            else:
                vprint("Invalid line:",line)
    return mutlist,linlist

def format_string_mutlin(mutlist,linlist):
    '''produce format string for neater printing of final table'''
    mutmaxlen = max(len(mut) for mut in mutlist)
    linmaxlen = max(len(lin) for lin in linlist)
    fmt = "%%%ds %%-%ds" % (mutmaxlen,linmaxlen)
    return fmt

def _main(args):
    '''mutlineage main'''

    mutlist,linlist = read_mutantfile(args.mutant)
    fmt = format_string_mutlin(mutlist,linlist)
    for mut,lin in zip(mutlist,linlist):
        vprint(fmt % (mut,lin),)

    seqs = covid.read_filter_seqfile(args)
    if args.N:
        ## just grab the first N (for debugging w/ shorter runs)
        seqs = it.islice(seqs,args.N+1)

    seqs = list(seqs)
    first,seqs = sequtil.get_first_item(seqs)
    MM = mutant.MutationManager(first.seq)
    for fullmatch in [True, False]:
        print("Exact matches:" if fullmatch else "Inclusive matches:")
        for mstring,lin in zip(mutlist,linlist):
            mpatt = mutant.Mutation.from_mstring(mstring,exact=fullmatch)
            matches = MM.filter_seqs_by_pattern(mpatt,seqs)
            cnt = count_lineages(matches)
            print(fmt % (mstring,lin),format_counter(cnt))


if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    _main(_args)
