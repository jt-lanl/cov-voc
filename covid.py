'''covid-specific utilities and hardcoded values'''

import sys
import os
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

DEFAULTSEQFILE="Latest.ipkl"

def default_seqfile(seqfilename=DEFAULTSEQFILE):
    '''
    return the default file for input sequences;
    hunt around in various directories until you find it
    '''
    for d in [os.getenv('DATA'),
              '.',
              'data',
              '..',
              '../data',
    ]:
        if not d:
            continue
        seqfile = Path(d) / seqfilename
        if seqfile.exists():
            return seqfile
        seqfile = Path(d) / (seqfilename + ".gz")
        if seqfile.exists():
            return seqfile
        
    return None

def datestring(yyyymmdd):
    ''' check format of dates; should be yyyy-mm-dd '''
    ## (but '.' is also allowed to indicated a default)
    if (not re.match(r'\d\d\d\d-\d\d-\d\d',yyyymmdd)
        and yyyymmdd != '.'):
        raise ValueError(f'Invalid date {yyyymmdd}; '
                         f'should be in yyyy-mm-dd format')
    return yyyymmdd

def corona_args(ap):
    ''' 
    call this in the getargs() function, 
    and these options will be added in 
    '''
    paa = ap.add_argument
    paa("--input","-i",type=Path,
        default=default_seqfile(),
        help="input file with aligned sequences (first is reference)")
    paa("--title",
        help="use this TITLE in plots")
    paa("--filterbyname","-f",nargs='+',
        help="Only use sequences whose name matches this pattern")
    paa("--xfilterbyname","-x",nargs='+',
        help="Do not use sequences whose name matches this pattern")
    paa("--dates","-d",nargs=2,type=datestring,
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

#### Routines for parsing sequence names


xpand_names = {
    ## dict to convert WHO lineage names to pango pattern
    'Alpha':  r'(B\.1\.1\.7)|(Q\.[1-9].*)',
    'Beta':   r'B\.1\.351',
    'Gamma':  r'P\.1.*',
    'Delta':  r'(B\.1\.617\.2)|(AY\.[1-9].*)',
    'Lambda': r'C\.37',
    'Mu':     r'B\.1\.621(\.1)?',
    'Omicron': r'(B\.1\.1\.529)|(BA\.[1-9].*)',
}

def get_isl(fullname):
    '''return EPI_ISL number from the sequence name'''
    epi_patt = re.compile(r"EPI_ISL_\d+")
    g = epi_patt.search(fullname)
    return g[0] if g else "X"

def get_lineage_from_name(name):
    '''get pango lineage by parsing the sequence name'''
    #re.sub is maybe more robust...but slower
    linpatt = re.sub(r".*EPI_ISL_\d+\.","",name)
    #try:
    #    tokens = name.split('.',6)
    #    lintok = tokens[6]
    #except IndexError:
    #    lintok = None
    #if linpatt != lintok:
    #    print("name=",name,"patt:",linpatt,"token:",lintok)
    return linpatt

def date_fromiso(s):
    if type(s) == datetime.date:
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError:
        if s == ".":
            return None
        return None #raise RuntimeError(f"Invalid Date {s}")

def date_from_seqname(s):
    #try:
    #    tokens = s.name.split('.')
    #    datestring = tokens[4]
    #except IndexError:
    #    datestring = s.name
    ## the following statement is more robust ... but slower!
    datestring = re.sub(".*(\d\d\d\d-\d\d-\d\d).*",r"\1",s.name)
    return date_fromiso(datestring)
    
def count_bad_dates(seqlist):
    return(sum(date_from_seqname(s) is None for s in seqlist))

def range_of_dates(seqlist):
    dates = [date_from_seqname(s) for s in seqlist]
    dates = [d for d in dates if d is not None]
    return min(dates).isoformat(),max(dates).isoformat()

def filter_by_date(seqs,fromdate,todate,keepfirst=False):
    '''input seqs is iterable (list or iterator); output is generator'''

    f_date = date_fromiso(fromdate)
    t_date = date_fromiso(todate)

    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        d = date_from_seqname(s)
        if not d:
            continue
        if f_date and f_date > d:
            continue
        if t_date and t_date < d:
            continue
        yield s


def mstring_brackets(mstring):
    '''make sure mstring has brackets around it'''
    mstring = re.sub(r'\s','',mstring) ## remove extra spaces between ssms
    mstring = re.sub(r'\s*ancestral\s*','',mstring) ## take out 'ancestral'
    mstring = re.sub(r'[\[\]]','',mstring) ## remove if already there
    mstring = f'[{mstring}]'
    return mstring

def mstring_fix(mstring):
    '''
    1. fix G-- vs --G 
    2. remove 'ancestral' string
    3. make a new mstring that treats G142x as G142.
    '''
    ### could we do D215A,+215AGY  => +214AAG,D215Y
    ## fix G-- vs --G
    if re.search('E156G,F157-,R158-',mstring):
        warnings.warn(f'fixing mstring with G-- pattern: {mstring}')
        mstring = re.sub('E156G,F157-,R158-','E156-,F157-,R158G',mstring)

    if re.search('ancestral',mstring):
        warnings.warn(f'removing "ancestral" from "{mstring}"')
        mstring = re.sub(r'\s*ancestral\s*','',mstring)

    ## fix G142G or G142_ or G142D -> G142.
    if re.search('G142',mstring):
        mstring = re.sub(r'G142[GD_]','G142.',mstring)

    ## fix N.10: I210V,N211-,L212I -> I210V,N211I,L212-
    #if re.search('I210V,N211-,L212I',mstring):
    #    warnings.warn(f'fixing N.10, V-I -> VI-')
    #    mstring = re.sub('N211-,L212I','N211I,L212-',mstring)

    return mstring

    #newmut = []
    #for ssm in mutant.Mutation(mstring):
    #    if ssm.site == 142 and ssm.mut in "GD_":
    #        newmut.append(mutant.SingleSiteMutation.from_ref_site_mut('G',142,'.'))
    #    else:
    #        newmut.append(ssm)
    #return str(mutant.Mutation(newmut))        

site_specifications = {
    "RBD"         : "330-521",
    "NTD"         : "14-292",
    "NTDss"       : "13-20,140-158,242-264",
    "xNTDss"      : "13-17,19,20,69,70,138-141,143-158,242-264",
    "NTDss_trunc" : "140-158,242-264",
    "NTDish"      : "13,18,20,69,70,141-144,152,153,157,242-244,253-257,262",
    "RBDish"      : "367,417,439,440,452,453,477,478,484,490,494,501,520,614",
    "RBDplus"     : "330-521,655,675,679,681",  ## obsolete, now use RBD+Furin
    "Furin"       : "655,675,677,679,681,950",
}

def spike_sites(sitespec):
    '''return list of integers, site numbers, corresponding to sitespec string'''
    ## sitespec may be of the form:
    ##   '13-20,140-158' -- integer list
    ##   'RBD' -- receptor binding domain 
    ##   'NTD', 'NTDss' -- N-terminal domain (supersite)
    ##   Combinations of the above using "+" and "-"
    ##   'RBD+NTD' -- include sites in either domain
    ##   'NTD-18' -- all sites in NTD except 18
    ##   'RBD+18' -- all sites in RBD plus site 18
    ## Problematic (don't do this, if you can avoid it):
    ##   Because '-' can indicate a 'remove site' or can indicate a numeric range
    ##   there are some unfortunate ambiguities
    ##   'RBD+1-10' -- Ambiguous:
    ##                 RBD plus 1 minus 10 /or/ RBD plus 1 through 10 ??
    ##                 You should get RBD plus 1 thru 10, but be careful
    ##    NTDss+69+70-18+RBD -- worse than ambiguous, will not include 70, will include 18
    ##                          because "70-18" -> "" since it's an invalid range
    ##    NTDss-18+69+70+RBD -- will give that right answer
    ##    69-70+NTDss-18+RBD -- may give right answer, but awkward formulation
    ##
    site_specs_internal = {key.upper(): value for key,value in site_specifications.items()}
    def get_intlist(s):
        return intlist.string_to_intlist( site_specs_internal.get(s.upper(),s) )
    xsite = set() ## set of sites to be excised
    sites = set() ## set of sites to be included
    for spec in sitespec.split('+'):
        if re.match(r'\d.*',spec):
            sites.update(get_intlist(spec))
            continue
        specx = spec.split('-')
        sites.update(get_intlist(specx[0]))
        for sx in specx[1:]:
            xsite.update( get_intlist(sx) )
    sites -= xsite
    return sorted(sites)
        

CONTINENTS = ["United-Kingdom",
              "Europe-minus-United-Kingdom",
              "North-America",
              "Asia",
              "Africa",
              "South-America",
              "Oceania",
]
ABBREV_CONTINENTS = {"United-Kingdom"           : "UK",
                     "Europe-minus-United-Kingdom": "Eu-UK",
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
        if "-minus-" in cx:
            c,x = cx.split("-minus-")
        else:
            c,x = cx,None
        ConExclude.append((cx,c,x))
    return ConExclude

def filename_prepend(pre,file):
    '''prepend a string to a file name; eg
       "pre","file" -> "prefile", but also
       "pre","dir/file" -> "dir/prefile"
    '''
    ## alt: dir,base = os.path.split(file)
    ##      return os.path.join(dir,pre+base)
    if not file:
        return file
    return re.sub(r"(.*/)?([^/]+)",r"\1"+pre+r"\2",file)

def get_title(args):
    ## Get title for plots and tables
    if args.title:
        return args.title
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


def get_first_item(items,keepfirst=True):
    '''
    get first item in iterable, and and put it back;
    works when the iterable is a list or an iterator
    '''
    warnings.warn("use sequtil.get not covid.get")
    return sequtil.get_first_item(items,keepfirst=keepfirst)

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

    first,seqs = sequtil.get_first_item(seqs)

    if "-" in first.seq and args.stripdashcols:
        seqs = sequtil.stripdashcols(first.seq,seqs)

    if not args.keeplastchar and "$" in first.seq:
        first.seq = re.sub('\$[^\$]*$','',first.seq)
        seqs = striplastchars(seqs,len(first.seq))

    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)

    return seqs

def striplastchars(seqs,seqlen):
    for s in seqs:
        s.seq = s.seq[:seqlen]
        s.seq = re.sub('\$','X',s.seq)
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
    first,seqs = sequtil.get_first_item(seqs)
    ref = first.seq
    for s in seqs:
        if X in s.seq:
            ss = list(s.seq)
            for n in range(len(ss)):
                if ss[n] == X:
                    ss[n] = ref[n]
            s.seq = "".join(ss)
        yield s

def lastdate_byfile(file,seqs=None):
    ## to get last date
    ## 1/ get modification date of input file
    ## 2/ if that doesn't work (eg, if file not found) get today's date
    ## 3/ if that doesn't work, get last date in dataset

    lastdate=None ## in case nothing works!
    try:
        mtime = os.path.getmtime(file)
        lastdate = datetime.date.fromtimestamp(mtime).isoformat()
    except FileNotFoundError:
        try:
            lastdate = datetime.date.today().isoformat()
        except:
            if seqs:
                seqs = list(seqs)
                _,lastdate = range_of_dates(seqs)
    return lastdate
        
def filter_seqs_by_date(seqs,args):

    if not args.days and not args.dates:
        return seqs
    if args.days and args.dates:
        raise RuntimeError("Cannot specify both --days AND --dates")

    if args.days:
        lastdate = lastdate_byfile(args.input,seqs)
        t = date_fromiso(lastdate) 
        f = t - datetime.timedelta(days=args.days)
        seqs = filter_by_date(seqs,f,t,keepfirst=True)
        
    if args.dates:
        seqs = filter_by_date(seqs,args.dates[0],args.dates[1],keepfirst=True)
        
    return seqs

def filter_seqs_by_pattern(seqs,args):
    '''input is iterable seqs; output is generator seqs'''

    keepers = []
    xcludes = []
    if args.filterbyname:
        for name in args.filterbyname:
            patt,wo,xpat = name.partition("-minus-")
            if patt in xpand_names:
                patt = xpand_names[patt]
            if patt != "Global":
                keepers.append(patt)
            if xpat:
                xcludes.append(xpat)
    if args.xfilterbyname:
        for name in args.xfilterbyname:
            if name in xpand_names:
                name = xpand_names[name]
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




