'''
given an input fasta file, the first of which is a reference
find the most common sequences, and express them in terms of a mutation string
'''
import sys
from collections import Counter,defaultdict
import argparse

import readseq
import sequtil
import wrapgen
import intlist
import mutant
import covid

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--mutants","-m",
        help="mutant string, such as '[W152R,N439K,D614G,P681R]'")
    paa("--sites","-s",
        help="restrict attention to just these sites")
    paa("-N",type=int,default=10,
        help="Maximum number of sequence patterns")
    paa("-C",type=int,default=3,
        help="Maximum number of countries per pattern")
    paa("--level","-l",type=int,default=2,
        help="Region Level: 1=Continent, 2=Country, 3=State_City, 4=Date")
    paa("--isl",action="store_true",
        help="Include sample ISL number for each type")
    paa("--output","-o",
        help="Write ouptut to fasta file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def get_country(fullname):
    '''return the name of the country in the full sequence name'''
    return get_region(2,fullname)

def get_region(level,fullname):
    '''return the name of the region (based on level) in the full sequence name'''
    m = fullname.split(".")
    return m[level] if m else None

def filtermutants(seqs,mstring):
    '''
    yield seqs that match mutant string;
    always yielding the first sequence
    '''
    first,seqs = sequtil.get_first_item(seqs)
    yield first
    mpatt = mutant.Mutation(mstring)
    MM = mutant.MutationManager(first.seq)
    for s in seqs:
        if MM.seq_fits_pattern(mpatt,s.seq):
            yield s

def main(args):
    '''commontypes main'''

    seqs = covid.read_filter_seqfile(args)

    if args.mutants:
        seqs = filtermutants(seqs,args.mutants)
        seqs = wrapgen.keepcount(seqs,"Sequences match mutation string:")

    seqlist = list(seqs)
    firstseq = seqlist[0].seq
    MM = mutant.MutationManager(firstseq)

    if args.sites:
        ## restrict attention to just these sites
        sites = intlist.string_to_intlist(args.sites)
        ndxlist = MM.indices_from_sitelist(sites)
        DOT='.' ## could be any character that doesn't occur in normal sequences
        for s in seqlist:
            if args.output:
                ## keep around the full seq if you are
                ## ultimately making an output fasta file
                s.fullseq = s.seq
            ## replace all non-essential sites with DOT
            slist = [DOT]*len(s.seq)
            for n in ndxlist:
                slist[n] = s.seq[n]
            s.seq = "".join(slist)


    isldict = defaultdict(list)
    if args.isl:
        for s in seqlist[1:]:
            isldict[s.seq].append( covid.get_isl(s.name) )

    counts = Counter([s.seq for s in seqlist[1:]])
    commonseqs = sorted(counts, key=counts.get, reverse=True)
    if args.N:
        commonseqs = commonseqs[:args.N]

    for seq in commonseqs:
        mut = MM.get_mutation(seq)
        if args.sites:
            mutx = [ssm for ssm in mut if DOT not in ssm.mut]
            mut = mutant.Mutation(mutx)
        mut.exact = False # just to avoid that '!'
        print(f"{counts[seq]:6d} {mut}",end="    ")
        if args.isl:
            print(", ".join(isldict[seq][:3]),end=" ")
            if len(isldict[seq]) > 3:
                print("...",end="")
        print()

        if args.C > 0:
            regions = Counter(get_region(args.level,s.name)
                              for s in seqlist[1:] if s.seq == seq)
            topregions = sorted(regions,key=regions.get,reverse=True)
            for t in topregions[:args.C]:
                print(f"       {regions[t]:6d} {t}")

    if args.output:
        seqs = []
        seqs.append(seqlist[0])
        for s in seqlist[1:]:
            if s.seq not in commonseqs:
                continue
            if args.sites:
                s.seq = s.fullseq
            seqs.append(s)
        vprint(len(seqs),"sequences written to output file:",args.output)
        readseq.write_seqfile(args.output,seqs)

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

    main(_args)
