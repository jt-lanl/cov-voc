DESCRIPTION='''
Make a barpolot of the coverages in the output from shiver.py
'''
import sys
import argparse
import re

import numpy as np
import matplotlib.pyplot as plt

def getargs():
    ap = argparse.ArgumentParser(description=DESCRIPTION)
    paa = ap.add_argument
    paa("input",
        help="input file is output from shiver run")
    #paa("--set","-s",type=int,default=1,
    #    help="how many set vaccines")
    paa("--lack",action="store_true",
        help="plot lack of coverage (1-C) instead")
    paa("--title","-t",
        help="title on top of plot")
    paa("--output","-o",
        help="write plot to file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def main(args):
    lines=[]
    skipIntro=True ## Skip ahead until "Table of Coverages"
    with open(args.input) as f:
        for line in f:
            vvprint("line=",line.strip())
            if line.strip() == "Table of Coverages":
                ## just keep skipping until we get to Table of Variants
                vvprint("OK: START READING NOW")
                skipIntro=False
            if skipIntro:
                continue
            
            try:
                ## Relevant lines are of the form: Continent vacc-name fraction (possibly multiple values, only use first)
                cont,vacc,*flist = line.strip().split()
                fo = flist[0]
                cont = re.sub("-w/o.*","",cont)
                cont = re.sub("-"," ",cont)
                lines.append((cont,vacc,float(fo)))
            except:
                pass

    for line in lines:
        vprint("%16s %5s %f" % line )

    continents = []
    vnames = []
    for cont,name,f in lines:        
        if cont not in continents:
            continents.append(cont)
        if name not in vnames:
            vnames.append(name)

    Nv = len(vnames)
    Nc = len(continents)

    frac = {c: dict() for c in continents}
    for c,name,f in lines:
        frac[c][name] = 1-f if args.lack else f

    figsize=(0.35*(Nv*Nc),3)
    plt.figure(figsize=figsize)

    for n,name in enumerate(vnames):
        plt.bar(n + (Nv+1)*np.arange(Nc),
                [frac[c][name] for c in continents],
                zorder=3,
                label=name)

    plt.xticks(np.arange(Nv/2-1/2,Nc*(Nv+1),Nv+1),
               labels=continents)
    plt.grid(axis='y',linestyle='--',linewidth=0.5)
    if args.lack:
        plt.ylim([0.01,1])
        plt.ylabel("Lack of Coverage")
        plt.yscale("log")
    else:
        plt.ylim([0,1])
        plt.ylabel("Coverage")
    plt.title(args.title)
    plt.legend(bbox_to_anchor=(1, 1),
               loc="upper left")
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
    
        
    

    
            

if __name__ == "__main__":

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    main(args)
    

