'''covid-specific utilities and hardcoded values'''

import sys
import re
import datetime
from pathlib import Path
from collections import Counter
import warnings

import readseq
import sequtil
import intlist
from hamming import hamming

#DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20210225.fasta"
#DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20201201.20210227.fasta.gz"
#DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20201201-20210302.fasta.gz"
#DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.20201001.20210309.Human.fasta"
#DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20210403.fasta.gz"
DEFAULTFASTA="Data/RX-REG_COMP.SPIKE.protein.Human.20210101-20210426.fasta.gz"

def corona_args(ap):
    ''' call this in the getargs() function, and these options will be added in '''
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        default=Path(DEFAULTFASTA),
        help="input fasta file with aligned sequences (first is master)")
    paa("--filterbyname","-f",
        help="Only use sequences whose name matches this pattern")
    paa("--xfilterbyname","-x",
        help="Do not use sequences whose name matches this pattern")
    paa("--keepdashcols",action="store_true",
        help="Do not strip columns with dash in ref sequence")
    #paa("--pattern","-p",
    #    help="If specified requires sequence name to match pattern")
    #paa("--xp",
    #    help="If specified, excludes sequence if name matches this pattern")

    #paa("--xvariant",
    #    help="Remove sequences that exhibit specified variant")
    #paa("--variant",
    #    help="Keep only sequences that exhibit specified variant")
    #paa("--noukvariant",action="store_true",
    #    help="Remove sequences that exhibit UK variant")
    #paa("--ukvariant",action="store_true",
    #    help="Keep only sequences that exhibit UK variant")
    #paa("--fixsiteseventy",action="store_true",
    #    help="Replace '--I' with 'I--' at sites 68-70")
    paa("--dates","-d",nargs=2,
        help="range of dates (two dates, yyyy-mm-dd format)")
    paa("--daysago",type=int,default=0,
        help="Consider date range from DAYSAGO days ago until today")

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

def fixsiteseventy(seqs):
    count=0
    for s in seqs:
        if s.seq[67:70] == "--I" or s.seq[67:70] == "-I-":
            count += 1
            s.seq = s.seq[:67] + "I--" + s.seq[70:]
            print("fix70: ",s.name,file=sys.stderr)
    return count

site_specifications = {
    "RBD"       : "330-521",
    "NTD"       : "14-20,140-158,242-264",
    "NTD-18"    : "14-17,19,20,140-158,242-264",
    "NTD+13-18" : "13-17,19,20,140-158,242-264",
    "NTD+13"    : "13-20,140-158,242-264",
    "NTD-TRUNC" : "140-158,242-264",
    "NTD+RBD"   : "14-20,140-158,242-264,330-521",
    "NTD-18+RBD": "14-17,19,20,140-158,242-264,330-521",
    "NTD+136970-18+RBD" : "13-17,19,20,69,70,140-158,242-264,330-521",
    "NTD+136970+RBD"    : "13-20,69,70,140-158,242-264,330-521",
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
    title = args.filterbyname or "Global"
    if title==".":
        title = "Global"
    if args.xfilterbyname:
        title = title + " w/o " + args.xp
    if 0: #args.noukvariant:
        title = title + " w/o UK-variant"
    return title

def read_seqfile(args,**kwargs):

    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    
    seqs = readseq.read_seqfile(args.input,badchar='X',**kwargs)
    vprint(len(seqs),"sequences read")

    clen = Counter([len(s.seq) for s in seqs])
    for l in clen:
        vprint(clen[l],"sequences of length",l)
    if len(clen)>1:
        warnings.warn("Not all sequences are the same length")

    return seqs

def filter_seqs(seqs,args):
    seqs = filter_seqs_by_date(seqs,args)
    seqs = filter_seqs_by_pattern(seqs,args)
    
    firstseq = seqs[0].seq
    if "-" in firstseq and not args.keepdashcols:
        warnings.warn("Stripping sites with dashes in first sequence")
        sequtil.stripdashcols(firstseq,seqs)
    return seqs

def filter_seqs_by_date(seqs,args):

    if args.daysago and args.dates:
        raise RuntimeError("Cannot specify both --daysago AND --dates")
    if not args.daysago and not args.dates:
        return seqs

    firstseq = seqs[:1]
    seqs = seqs[1:]

    if args.daysago:
        t = datetime.date.today()
        f = t - datetime.timedelta(days=args.daysago)
        args.dates = f,t

    if args.dates:
        seqs = sequtil.filter_by_date(seqs,args.dates[0],args.dates[1],keepfirst=True)

    return firstseq + seqs

def filter_seqs_by_pattern(seqs,args):

    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)

    if args.filterbyname and args.filterbyname != "Global":
        patt,wo,xpatt = args.filterbyname.partition("-w/o-")
        seqs = sequtil.filter_by_pattern(seqs,patt,keepfirst=True)
        vprint(len(seqs),"sequences fit pattern:",patt)
        if xpatt:
            seqs = sequtil.filter_by_pattern_exclude(seqs,xpatt,keepfirst=True)
            vprint(len(seqs),"sequences after removing x-pattern:",xpatt)

    if args.xfilterbyname:
        seqs = sequtil.filter_by_pattern_exclude(seqs,args.xfilterbyname,keepfirst=True)
        vprint(len(seqs),"sequences after removing x-pattern:",args.xfilterbyname)

    ## if there are any dashes in first sequence, then strip those
    ## columns from all the sequences
    if "-" in seqs[0].seq:
        warnings.warn("Sites with dashes in first sequence")

    if 0: #args.fixsiteseventy:
        fixes = fixsiteseventy(seqs)
        print("Fixed sites 68-70 for",fixes,"sequences")
        if fixes>0:
            readseq.write_fasta("RX.fasta",seqs)                                
            return

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




