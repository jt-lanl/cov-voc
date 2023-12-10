'''Plot output from newmut -- list of mutations as they emerge'''

import re
import argparse
from collections import defaultdict,Counter
import numpy as np
import matplotlib.pyplot as plt
import verbose as v

## Default color cycle: blue, orange, etc
_prop_cycle = plt.rcParams['axes.prop_cycle']
_colors = _prop_cycle.by_key()['color']
def nth_color(n):
    return _colors[n % len(_colors)]

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input file is output of newmut")
    paa("--top","-t",type=int,default=0,
        help="Only plot the top T mutations")
    paa("--stat","-s",choices=('par','lin','cnt'),default='lin',
        help="Which statistic to use (parent, lineage, counts)")
    paa("--rbd",action="store_true",
        help="Plot only RBD region (328-528) to plotfile")
    paa("--output","-o",
        help="Write output to file")
    paa("--plotfile","-p",
        help="write plot to file")
    paa("--title",
        help="plot title")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def parseline(line):
    '''each line of input file'''
    tokens = line.strip().split()
    if len(tokens) < 6:
        v.vprint('not six tokens in: ',tokens)
        return None
    ssm = tokens[0]
    site = int(re.sub('\D*(\d+)\D*',r'\1',ssm))
    npar = int(tokens[1])
    nlin = int(tokens[2])
    nseq = int(tokens[3])
    return site,ssm,npar,nlin,nseq

def bysite(parsedlines):
    prev_site = None
    bysitelines=[]
    parsedlines = list(parsedlines) + [(0,0,0,0,0)]
    for site,ssm,npar,nlin,nseq in parsedlines:
        if site != prev_site:
            if prev_site:
                bysitelines.append((prev_site,None,snpar,snlin,snseq))
            snpar = snlin = snseq = 0
            prev_site = site
        snpar += npar
        snlin += nlin
        snseq += nseq

    return bysitelines

def ssmlabel(ssmlist):
    xssmlist = [re.sub(r'[A-Z]*\d+','',ssm) for ssm in ssmlist]
    noins = [ssm for ssm in ssmlist if ssm[0]!='+']
    if noins:
        initletter = noins[0][0] + ": "
    else:
        initletter = ''
    label = initletter+",".join(xssmlist)
    if len(label) > 20:
        label = label[:20]
        label = re.sub(r'[^,*]$','',label) + "..."
    return label

def _main(args):
    '''main'''
    v.vprint(args)

    parsedlines = []
    with open(args.input) as fin:
        prev_ssm = None
        for line in fin:
            p = parseline(line)
            if not p:
                continue
            site,ssm,npar,nlin,nseq = p
            if ssm == prev_ssm:
                continue
            prev_ssm = ssm
            parsedlines.append((site,ssm,npar,nlin,nseq))

    if args.stat == 'par': VALNDX = 2
    if args.stat == 'lin': VALNDX = 3
    if args.stat == 'cnt': VALNDX = 4
    v.vprint(f'{VALNDX=}')
    ## get value as fcn of ssm
    val_byssm = {line[1]:line[VALNDX] for line in parsedlines}

    ## get list of ssm's or each site
    ssms_bysite = defaultdict(list)
    for line in parsedlines:
        site = line[0]
        ssm = line[1]
        ssms_bysite[site].append(ssm)

    val_bysite = Counter()
    for line in parsedlines:
        site = line[0]
        value = line[VALNDX]
        val_bysite[site] += value

    ## can get top sites using val_bysite.most_common(args.top)

    ## sum up all the values in each site (needed for getting top sites)
    parsedlines = bysite(parsedlines)

    if args.output:
        ## write consolidated data file (numbers only)
        with open(args.output,"w") as fout:
            for site,ssm,npar,nlin,nseq in parsedlines:
                print(site,npar,nlin,nseq,file=fout)
        return
                
    sitelist = []
    values = []
    for line in parsedlines:
        sitelist.append(line[0])
        values.append(line[VALNDX])

    if args.top:
        thresh=1
        if args.top < len(values):
            sv = sorted(values)
            thresh = max([thresh, sv[-args.top+1]])
        keep = [i for i in range(len(values)) if values[i] > thresh]
        v.vprint('Thresh',thresh,'reduces',len(values),"to",len(keep))
        sitelist = [sitelist[i] for i in keep]
        values = [values[i] for i in keep]

    if args.rbd:
        newsitelist = list(range(328,529))
        newvalues = [0]*len(newsitelist)
        for ndx,site in enumerate(newsitelist):
            try:
                k = sitelist.index(site);
                newvalues[ndx] = values[k]
            except ValueError:
                pass
        sitelist = newsitelist
        values = newvalues

    try:
        max_ssms_persite = max(len(ssms_bysite[site]) for site in sitelist)
        v.vprint('max ssms:',max_ssms_persite)
    except ValueError:
        v.print(_args)
        v.print(f'{len(sitelist)=}')
        raise ValueError
        
    plt.figure(figsize=(18,5))
    textspace = max(values)*0.02
    ssmstrlist = []
    for ndx,site in enumerate(sitelist):
        bot = 0
        ssmlines = sorted(ssms_bysite[site],
                          key=val_byssm.get,
                          reverse=True)
        ssmstr = ssmlabel(ssmlines)
        ssmstrlist.append(ssmstr)
        for n,ssm in enumerate(ssmlines):
            val = val_byssm[ssm]
            plt.bar(ndx,val,bottom=bot,color=nth_color(n))
            bot = bot+val
        plt.text(ndx,bot+textspace,ssmstr,
                 horizontalalignment='center',
                 rotation='vertical',
                 fontsize='x-small')

    plt.xlabel('site')
    if args.rbd:
        newsitelist = [str(site) for site in sitelist if site%2 == 0]
        newticks = [2*i for i in range(len(newsitelist))]
        plt.xticks(newticks,newsitelist,
                   rotation='vertical')
    else:
        plt.xticks(range(len(sitelist)),
                   [str(s) for s in sitelist],
                   rotation='vertical')
    plt.xlim([-1,len(sitelist)])
    plt.ylabel(f'recurrences ({args.stat})')
    plt.gca().spines[['right', 'top']].set_visible(False)

    if args.title:
        plt.title(args.title)
    plt.tight_layout()
    if args.plotfile:
        plt.savefig(args.plotfile)
    else:
        plt.show()    


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
