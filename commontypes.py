'''
given an input fasta file, the first of which is a reference
find the most common sequences, and express them in terms of a mutation string
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
import mutant

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        help="input fasta file with all spike sequences")
    paa("--mutants","-m",
        help="mutant string, such as '[W152R,N439K,D614G,P681R]'")
    paa("--sites","-s",
        help="restrict attention to just these sites")
    paa("--filtername","-f",
        help="pattern for filtering by name")
    paa("--keepx",action="store_true",
        help="keep sequences with bad characters")
    paa("--dates",nargs=2,
        help="range of dates (two dates, yyyy-mm-dd format)")
    paa("-N",type=int,default=10,
        help="Maximum number of sequence patterns")
    paa("-C",type=int,default=3,
        help="Maximum number of countries per pattern")
    paa("--level","-l",type=int,default=2,
        help="Region Level: 1=Continent, 2=Country, 3=State, 4=City")
    paa("--stripdashcols",action="store_true",
        help="Strip dash columns")
    paa("--output","-o",
        help="Write ouptut to fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def get_country(fullname):
    m = fullname.split(".")
    if m:
        return m[2]
    else:
        return None

def get_region(level,fullname):
    m = fullname.split(".")
    if m:
        return m[level]
    else:
        return None
    


def filterseqs(args,seqlist):
    ''' pull out a sublist of the sequence list, based on 
    various options in the args structure '''

    firstseq = seqlist[0].seq

    if args.dates:
        seqlist = sequtil.filter_by_date(seqlist,
                                         args.dates[0],args.dates[1],
                                         keepfirst=True)
        vprint(len(seqlist),"sequences in date range:",args.dates)

    if args.filtername:
        seqlist = seqlist[:1] + [s for s in seqlist[1:]
                                 if args.filtername in s.name]
        vprint(len(seqlist),"sequences filtered by:",args.filtername)

    if args.stripdashcols and "-" in firstseq:
        vprint("Stripping sites with dashes in first sequence...",end="")
        sequtil.stripdashcols(firstseq,seqlist)
        vprint("ok")
        if "-" in seqlist[0].seq:
            raise RuntimeError("strip dash failed!")

    firstseq = seqlist[0].seq
    
    if "-" in firstseq:
        warnings.warn("dashes in reference sequence")

    if not args.keepx:
        seqlist = seqlist[:1] + [s for s in seqlist[1:]
                                 if "X" not in s.seq]
        vprint(len(seqlist),"sequences without X")

    return seqlist

def filtermutants(seqlist,mutstring):
    firstseq = seqlist[0].seq
    mutlist = mutant.parse_mutline(mutstring)
    for mut in mutlist:
        if mut.ref != firstseq[mut.site-1]:
            vprint("Warning: at site",mut.site,":",
                   mut.ref,"should be",firstseq[mut.site-1])
    seqpattern = "".join(mut.mut for mut in mutlist)
    sites = [mut.site for mut in mutlist]

    fullpatt = ["."] * len(firstseq)
    for c,s in zip(seqpattern,sites):
        fullpatt[s-1] = c
    fullpatt = "".join(fullpatt)

    seqlist = seqlist[:1] + [s for s in seqlist[1:]
                             if re.match(fullpatt,s.seq)]

    return seqlist

def main(args):

    seqlist = readseq.read_seqfile(args.input)
    vprint(len(seqlist),"sequences read")

    seqlist = filterseqs(args,seqlist)
    vprint(len(seqlist),"sequences after filtering")

    if args.mutants:
        seqlist = filtermutants(seqlist,args.mutants)
        vprint(len(seqlist),"sequences match mutant string:",args.mutants)

    firstseq = seqlist[0].seq
        
    if args.sites:
        sites = intlist.string_to_intlist(args.sites)
        for s in seqlist:
            s.fullseq = s.seq
            ## replace char with "." if not in sites array
            slist = ["."]*len(s.seq)
            for n in sites:
                slist[n-1] = s.fullseq[n-1]
            s.seq = "".join(slist)

    counts = Counter([s.seq for s in seqlist[1:]])
    commonseqs = sorted(counts, key=counts.get, reverse=True)

    for seq in commonseqs[:args.N]:
        mutlist=[]
        for n,(refchar,mutchar) in enumerate(zip(firstseq,seq)):
            if refchar != mutchar:
                mutlist.append( f"{refchar}{n+1}{mutchar}" )
        mutstr = "[" + ",".join(mutlist) + "]"
        print(f"{counts[seq]:6d} {mutstr}")

        regions = Counter(get_region(args.level,s.name)
                            for s in seqlist[1:] if s.seq == seq)
        topregions = sorted(regions,key=regions.get,reverse=True)
        for t in topregions[:3]:
            print(f"       {regions[t]:6d} {t}")

    if args.output:
        seqs = []
        seqs.append(seqlist[0])
        for s in seqlist:
            if args.N>0 and s.seq not in commonseqs:
                continue
            if args.sites:
                s.seq = s.fullseq
            seqs.append(s)
        vprint(len(seqs),"sequences written to output file:",args.output)
        readseq.write_seqfile(args.output,seqs)
        
        


if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

