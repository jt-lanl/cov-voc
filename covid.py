'''covid-specific utilities and hardcoded values'''

import sys
import re
import datetime
import itertools
from pathlib import Path
from collections import Counter,namedtuple
import pickle
import warnings

import wrapgen
import readseq
import sequtil
import intlist

DEFAULTFASTA="Latest.ipkl.gz"

def corona_args(ap):
    ''' call this in the getargs() function, and these options will be added in '''
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        default=Path(DEFAULTFASTA),
        help="input file with aligned sequences (first is reference)")
    paa("--filterbyname","-f",nargs='+',
        help="Only use sequences whose name matches this pattern")
    paa("--xfilterbyname","-x",nargs='+',
        help="Do not use sequences whose name matches this pattern")
    paa("--dates","-d",nargs=2,
        help="Only use seqs in range of dates (two dates, yyyy-mm-dd format)")
    paa("--days",type=int,default=0,
        help="Consider date range of DAYS days ending on the last sampled date")
    paa("--stripdashcols",action="store_true",
        help="Strip columns with dash in reference sequence")
    paa("--keeplastchar",action="store_true",
        help="Do not strip final stop codon from end of sequences")
    paa("--keepx",action="store_true",
        help="Keep sequences that include bad characters, denoted X")

    return

def get_isl(fullname):
    epi_patt = re.compile(r"EPI_ISL_\d+")
    g = epi_patt.search(fullname)
    return g[0] if g else "X"


site_specifications = {
    "RBD"       : "330-521",
    "NTD"       : "14-292",
    "NTDSS"       : "13-20,140-158,242-264",
    "NTDSS-18"    : "13-17,19,20,140-158,242-264",
    "NTDSS-13-18" : "14-17,19,20,140-158,242-264",
    "NTDSS-13"    : "14-20,140-158,242-264",
    "NTDSS-TRUNC" : "140-158,242-264",
    "NTDSS+RBD"   : "13-20,140-158,242-264,330-521",
    "NTDSS-13-18+RBD": "14-17,19,20,140-158,242-264,330-521",
    "NTDSS-18+RBD": "13-17,19,20,140-158,242-264,330-521",
    "NTDSS+6970-18+RBD" : "13-17,19,20,69,70,140-158,242-264,330-521",
    "NTDSS+6970+RBD"    : "13-20,69,70,140-158,242-264,330-521",
    "NTDISH"    : "13,18,20,69,70,141,142,143,144,152,153,157,242,243,244,253,254,255,256,257,262",
    "RBDISH"    : "367,417,439,440,452,453,477,478,484,490,494,501,520,614",
    "NTDSS-18+RBDPLUS" : "13-17,19,20,140-158,242-264,330-521,655,675,679,681",
}
def spike_sites(sitespec,remove=None):
    '''return a list of site numbers based on spec'''
    sitespec = sitespec.upper()
    if sitespec in site_specifications:
        sitenums = intlist.string_to_intlist(site_specifications[sitespec])
    else:
        sitenums = intlist.string_to_intlist(sitespec)
    if remove:
        for r in remove:
            sitenums.remove(r)
            
    return sitenums

CONTINENTS = ["United-Kingdom",
              "Europe-w/o-United-Kingdom",
              "North-America",
              "Asia",
              "Africa",
              "South-America",
              "Oceania",
]
ABBREV_CONTINENTS = {"United-Kingdom"           : "UK",
                     "Europe-w/o-United-Kingdom": "Eu-UK",
                     "North-America"            : "NAmer",
                     "Asia"                     : "Asia",
                     "Africa"                   : "Africa",
                     "South-America"            : "SAmer",
                     "Oceania"                  : "Ocean",
}
    
def parse_continents(withglobal=False):
    ConExclude=[]
    if withglobal:
        ConExclude.append(("Global","Global",None))
    for cx in CONTINENTS:
        if "-w/o-" in cx:
            c,x = cx.split("-w/o-")
        else:
            c,x = cx,None
        ConExclude.append((cx,c,x))
    return ConExclude

def get_title(args):
    ## Get title for plots and tables
    title = "+".join(args.filterbyname) if args.filterbyname else "Global"
    if title==".":
        title = "Global"
    if args.xfilterbyname:
        title = title + " w/o " + "+".join(args.xfilterbyname)
    return title

def summarizeseqlengths(seqlist,args):
    #assert isinstance(seqlist,list)
    clen = Counter([len(s.seq) for s in seqlist])
    if args.verbose:
        for l in clen:
            print(clen[l],"sequences of length",l,file=sys.stderr)
    if len(clen)>1:
        warnings.warn("Not all sequences are the same length")
    

def get_first_item(items,putitemback=True):
    '''get first item in iterable, and and put it back'''
    if isinstance(items,list):
        first = items[0]
        if not putitemback:
            items = items[1:]
    else:
        first = next(items)
        if putitemback:
            items = itertools.chain([first],items)
    return first,items
        
def checkseqlengths(seqs):
    first,seqs = get_first_item(seqs)
    seqlen = len(first.seq)
    for s in seqs:
        if len(s.seq) != seqlen:
            raise RuntimeError(f"Sequence {s.name} has inconsistent length: {len(s.seq)} vs {seqlen}")
        yield s

def read_filter_seqfile(args,**kwargs):
    seqs = read_seqfile(args,**kwargs)
    seqs = filter_seqs(seqs,args)
    return seqs    
        
def read_seqfile(args,**kwargs):
    seqs = readseq.read_seqfile(args.input,badchar='X',**kwargs)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences read:")
    return seqs

def fix_seqs(seqs,args):

    firstseq,seqs = get_first_item(seqs)

    if "-" in firstseq.seq and args.stripdashcols:
        seqs = sequtil.stripdashcols(firstseq.seq,seqs)

    if not args.keeplastchar and firstseq.seq and firstseq.seq[-1] in "$X":
        seqs = striplastchar(seqs)

    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)

    return seqs

def striplastchar(seqs):
    for s in seqs:
        s.seq = s.seq[:-1]
        yield s

def filter_seqs(seqs,args):
    ## by date first so multiple runs will have the same date range
    ## with --days option
    seqs = filter_seqs_by_date(seqs,args)
    seqs = filter_seqs_by_pattern(seqs,args)
    seqs = fix_seqs(seqs,args)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences filtered:")

    return seqs

def xrepair(seqs,X='X'):
    '''replace all occurrences of X with the ancestral form in first sequence'''
    ## This seems dangerous!! Use with care ... or not at all.
    first,seqs = get_first_item(seqs)
    yield first
    for s in seqs:
        if X in s.seq:
            ss = list(s.seq)
            for n in range(len(ss)):
                if ss[n] == X:
                    ss[n] = ref[n]
            s.seq = "".join(ss)
        yield s

def filter_seqs_by_date(seqs,args):

    if not args.days and not args.dates:
        return seqs
    if args.days and args.dates:
        raise RuntimeError("Cannot specify both --days AND --dates")

    if args.days:
        seqs = list(seqs)
        _,lastdate = sequtil.range_of_dates(seqs)
        t = sequtil.date_fromiso(lastdate) ## not: datetime.date.today()
        f = t - datetime.timedelta(days=args.days)
        args.dates = f,t

    if args.dates:
        seqs = sequtil.filter_by_date(seqs,args.dates[0],args.dates[1],keepfirst=True)
        
    return seqs

def filter_seqs_by_pattern(seqs,args):
    '''input is iterable seqs; output is generator seqs'''

    keepers = []
    xcludes = []
    if args.filterbyname:
        for name in args.filterbyname:
            patt,wo,xpat = name.partition("-w/o-")
            if patt != "Global":
                keepers.append(patt)
            if xpat:
                xcludes.append(xpat)
    if args.xfilterbyname:
        for name in args.xfilterbyname:
            xcludes.append(name)

    ## Use r"\."+ to ensure that names have to be preceeded by a dot
    keepers = [r"\."+patt for patt in keepers]
    xcludes = [r"\."+xpat for xpat in xcludes]

    if keepers:
          seqs = sequtil.filter_by_patternlist(seqs,keepers,
                                               keepfirst=True)
    if xcludes:
          seqs = sequtil.filter_by_patternlist_exclude(seqs,xcludes,
                                                       keepfirst=True)
            
    return seqs

SARS_REGIONS = '''
SP     1   13
NTD   27  292
RBD  330  521
RBM  438  506
SD1  528  590
SD2  591  679
#S1   685  685
#S2   686  686
FP   816  833
HR1  908  985
CH   986 1035
CD  1076 1140
HR2 1162 1204
TM  1214 1236
'''

MERS_REGIONS='''
SP      1       18
NTD     18      350
RBD     382     588
RBM     483     566
SD1     595     656
SD2     657     747
#S1/S2   748     749     The cleavage site
UH      816     851     
FP      880     900     
CR      901     991
HR1     992     1054   
CH      1055    1109
BH      1110    1149
SD3     1150    1206
HR2     1246    1294
TM      1295    1329
'''

def get_srlist(virus="SARS"):
    ''' provides a tuple of three lists: n,x,y
        n = name of region
        x = start site
        y = stop site
    '''
    n=[]
    x=[]
    y=[]
    regions = MERS_REGIONS if virus.lower() == "mers" else SARS_REGIONS
    for sr in regions.split('\n'):
        srs = sr.split()
        if len(srs) == 3 and srs[0][0]!='#':
            n.append(srs[0])
            x.append(int(srs[1]))
            y.append(int(srs[2]))
    return n,x,y

## EXAMPLE USAGE:
#        f = max(e)/2
#        dy = -0.04*f
#        dyl = -0.12*f
#                
#        for srx,sry,srn in zip(srlistx,srlisty,srlistn):
#            plt.plot([srx,sry],[dy,dy],label=srn,linewidth=4)
#            delta = 8 if sry-srx < 15 else 3 if sry-srx < 25 else 0
#            plt.text(srx-delta,dyl,srn)
#        plt.legend(loc="upper left",title="Regions",bbox_to_anchor=(1,1))




