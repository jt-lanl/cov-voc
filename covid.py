'''covid-specific utilities and hardcoded values'''

import sys
import re
import datetime
import itertools
from pathlib import Path
from collections import Counter,namedtuple
import pickle
import warnings

import readseq
import sequtil
import intlist
import mutant

DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20210728.ipkl.gz"

def corona_args(ap):
    ''' call this in the getargs() function, and these options will be added in '''
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        default=Path(DEFAULTFASTA),
        help="input fasta file with aligned sequences (first is master)")
    paa("--filterbyname","-f",nargs='+',
        help="Only use sequences whose name matches this pattern")
    paa("--xfilterbyname","-x",nargs='+',
        help="Do not use sequences whose name matches this pattern")
    paa("--keepdashcols",action="store_true",
        help="Do not strip columns with dash in ref sequence")
    paa("--keeplastchar",action="store_true",
        help="Do not strip final stop codon from end of sequences")
    paa("--dates","-d",nargs=2,
        help="range of dates (two dates, yyyy-mm-dd format)")
    paa("--days",type=int,default=0,
        help="Consider date range of DAYS days ending on the last sampled date")
    paa("--fixsiteseventy",action="store_true",
        help="Sites 68-70 should be I--, not -I- or --I")

    return

def get_isl(fullname):
    epi_patt = re.compile(r"EPI_ISL_\d+")
    g = epi_patt.search(fullname)
    return g[0] if g else "X"


def UKvariant():
    UK = {  69:   '-',
            70:   '-',
            144:  '-',
            501:  'Y',
            570:  'D',
            681:  'H',
            716:  'I',
            982:  'A',
            1118: 'H',
            }
    return UK

def ZAvariant():
    ZA = {   80:  'A',
            215:  'G',
            242: '-',
            243: '-',
            244: '-',
            417: 'N',
            484: 'K',
            501: 'Y',
            701: 'V',
    }
    return ZA

def BRvariant():
    BR = {  20: 'N',
            26: 'S',
            138: 'Y',
            417: 'T',
            484: 'K',
            501: 'Y',
            655: 'Y',
            1027: 'I',
            1176: 'F',
            }
    return BR
            
def USvariant():
    US = { 13: 'I',
           152: 'C',
           452: 'R',
           }
    return US

def UGvariant():
    UG = { 157: "L",
           367: "F",
           613: "H",
           614: "D", # G
           681: "R"
           }
    return UG

the_variant = {
    'UK': UKvariant,
    'ZA': ZAvariant,
    'BR': BRvariant,
    'US': USvariant,
    'UG': UGvariant,
    }

def fixsiteseventy_gen(seqs,args):
    for s in seqs:
        if s.seq[67:70] == "--I" or s.seq[67:70] == "-I-":
            s.seq = s.seq[:67] + "I--" + s.seq[70:]
        yield s

def fixsiteseventy(seqs,args):
    count=0
    for s in seqs:
        if s.seq[67:70] == "--I" or s.seq[67:70] == "-I-":
            count += 1
            s.seq = s.seq[:67] + "I--" + s.seq[70:]
            if args.verbose > 1:
                print("fix70: ",s.name,file=sys.stderr)
    return count

site_specifications = {
    "RBD"       : "330-521",
    "NTD"       : "13-20,140-158,242-264",
    "NTD-18"    : "13-17,19,20,140-158,242-264",
    "NTD-13-18" : "14-17,19,20,140-158,242-264",
    "NTD-13"    : "14-20,140-158,242-264",
    "NTD-TRUNC" : "140-158,242-264",
    "NTD+RBD"   : "13-20,140-158,242-264,330-521",
    "NTD-13-18+RBD": "14-17,19,20,140-158,242-264,330-521",
    "NTD-18+RBD": "13-17,19,20,140-158,242-264,330-521",
    "NTD+6970-18+RBD" : "13-17,19,20,69,70,140-158,242-264,330-521",
    "NTD+6970+RBD"    : "13-20,69,70,140-158,242-264,330-521",
    "NTDISH"    : "13,18,20,69,70,141,142,143,144,152,153,157,242,243,244,253,254,255,256,257,262",
    "RBDISH"    : "367,417,439,440,452,453,477,478,484,490,494,501,520,614",
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
    

def get_first_item(seqs):
    '''get first item in iterable, and and put it back'''
    if isinstance(seqs,list):
        first = seqs[0]
    else:
        first = next(seqs)
        seqs = itertools.chain([first],seqs)
    return first,seqs
        
def checkseqlengths(seqs):
    first,seqs = get_first_item(seqs)
    seqlen = len(first.seq)
    for s in seqs:
        if len(s.seq) != seqlen:
            raise RuntimeError(f"Sequence {s.name} has inconsistent length: {len(s.seq)} vs {seqlen}")
        yield s

def read_seqfile(args,**kwargs):
    seqs = readseq.read_seqfile(args.input,badchar='X',**kwargs)
    return seqs

def fix_seqs(seqs,args):

    firstseq,seqs = get_first_item(seqs)

    if "-" in firstseq.seq and not args.keepdashcols:
        warnings.warn("Stripping sites with dashes in first sequence")
        sequtil.stripdashcols(firstseq.seq,seqs)

    if args.fixsiteseventy:
        seqs = fixsiteseventy_gen(seqs,args)

    if not args.keeplastchar and firstseq.seq and firstseq.seq[-1] in "$X":
        seqs = striplastchar(seqs)

    return seqs

def striplastchar(seqs):
    for s in seqs:
        s.seq = s.seq[:-1]
        yield s

def filter_seqs(seqs,args):
    seqs = filter_seqs_by_date(seqs,args)
    seqs = filter_seqs_by_pattern(seqs,args)
    seqs = fix_seqs(seqs,args)
    return seqs

def filter_seqs_by_date(seqs,args):

    if args.days and args.dates:
        raise RuntimeError("Cannot specify both --days AND --dates")
    if not args.days and not args.dates:
        return seqs

    if args.days:
        seqs = list(seqs)
        lastdate = sequtil.range_of_dates(seqs)[1]
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


def init_lineages(filename,firstseq):
    NamedPattern = namedtuple('NamedPattern',['name','pattern','color'])
    lineages = []
    if not filename:
        return lineages
    with open(filename) as f:
        for line in f:
            line = re.sub("#.*","",line).strip()
            if not line:
                #ignore empty and commented-out lines
                continue
            ## Match: Color [Mutation]! Name, with "!" optional and Name optional
            m = re.match("(\S+)\s+(\[.*\])(!?)\s*(\S*).*",line)
            if not m:
                warnings.warn(f"No match: {line}")
                continue
            mpattern = mutant.Mutation(m[2]).regex_pattern(firstseq,exact=bool(m[3]))
            lineages.append( NamedPattern(name=m[4],pattern=mpattern,color=m[1]) )
    return lineages

def match_lineages(lineages,fullseq):
    lineage_name=""
    for lineage in lineages:
        if re.match(lineage.pattern,fullseq):
            lineage_name = lineage.name
            break ## match first available
    return lineage_name #,lineage_color

def match_lineage_name_color(lineages,fullseq):
    lineage_name="other"
    lineage_color="Gray"
    for lineage in lineages:
        if re.match(lineage.pattern,fullseq):
            lineage_name = lineage.name
            lineage_color = lineage.color
            break ## match first available
    return lineage_name,lineage_color



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




