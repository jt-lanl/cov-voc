import sys
import re
import argparse

import covid
import colornames

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--colormut",required=True,
        help="name of color mutation file (mutation_string,lineage_name) are 2nd,3rd columns")
    paa("--usehex",action="store_true",
        help="use six-character hex-codes instead of X11 colornames")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def main(args):

    seqlist = covid.read_seqfile(args)
    vprint(len(seqlist),"sequences read")
    seqlist = covid.filter_seqs(seqlist,args)
    vprint(len(seqlist),"sequences after filtering")

    firstseq = seqlist[0].seq

    lineages = covid.init_lineages(args.colormut,firstseq)
    Table = []
    for s in seqlist[1:]:
        n,c = covid.match_lineage_name_color(lineages,s.seq)
        if args.usehex:
            try:
                c = colornames.tohex(c)
            except:
                pass
            if n == "other":
                c = "#DDDDDD"
        Table.append((s.name,n,c))

    for seqname,n,c in Table:
        print(seqname,n,c)
        

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

