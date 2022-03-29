'''covid-specific utilities and hardcoded values'''

import os
import re
import argparse
import datetime
import itertools as it
from pathlib import Path
import warnings

import wrapgen
import readseq
import sequtil
import intlist
from verbose import verbose as v

MAX_TITLE_LENGTH=60 ## truncate long title names
ISO_DATE_REGEX = re.compile(r'\d\d\d\d-\d\d(-\d\d)?')
ISO_DATE_REGEX_DOTS = re.compile(r'\.(\d\d\d\d-\d\d-\d\d)\.')
EPI_ISL_REGEX = re.compile(r'EPI_ISL_\d+')
LINEAGE_REGEX = re.compile(r'EPI_ISL_\d+\.(.*)')

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


def corona_args(ap):
    '''
    call this in the getargs() function,
    and these options will be added in
    '''

    ### Check --dates option; is format okay? are they in order?
    def datestring(yyyymmdd):
        ''' check format of dates; should be yyyy-mm-dd '''
        ## (but '.' is also allowed to indicated a default)
        if (not re.match(r'\d\d\d\d-\d\d-\d\d',yyyymmdd)
            and yyyymmdd != '.'):
            raise ValueError(f'Invalid date {yyyymmdd}; '
                             f'should be in yyyy-mm-dd format')
        return yyyymmdd

    class CheckDates(argparse.Action):
        ''' check that the dates are good and that they are in order '''
        def __call__(self,parser,namespace,dates,option_string=None):
            dates = [datestring(date) for date in dates]
            if not any(bool(date == ".") for date in dates):
                if dates[0] > dates[1]:
                    parser.error(f'Dates out of order: {dates}')
            setattr(namespace,self.dest,dates)

    paa = ap.add_argument
    #faa = ap.add_argument_group('File input options').add_argument
    paa("--input","-i",type=Path,
        default=default_seqfile(),
        help="input file with aligned sequences (first is reference)")
    paa("--nseq",type=int,default=0,
        help="read at most NSEQ sequences")
    paa("--filterbyname","-f",nargs='+',
        help="Only use sequences whose name matches this pattern")
    paa("--xfilterbyname","-x",nargs='+',
        help="Do not use sequences whose name matches this pattern")
    paa("--dates","-d",nargs=2,type=datestring,action=CheckDates,
        help="Only use seqs in range of dates (two dates, yyyy-mm-dd format)")
    paa("--days",type=int,default=0,
        help="Consider date range of DAYS days ending on the last sampled date")
    paa("--stripdashcols",action="store_true",
        help="Strip columns with dash in reference sequence")
    paa("--keeplastchar",action="store_true",
        help="Do not strip final stop codon from end of sequences")
    paa("--keepx",action="store_true",
        help="Keep sequences that include bad characters, denoted X")
    paa("--skipx",action="store_false",dest='keepx',
        help="Skip sequences that include bad characters, denoted X")
    ap.set_defaults(keepx=False)
    paa("--title",
        help="use this TITLE in plots")

#### Routines for parsing sequence names

xpand_WHO_Pangolin = {
    ## dict to convert WHO lineage names to pango pattern
    'Alpha':  r'(B\.1\.1\.7)|(Q\.[1-9].*)',
    'Beta':   r'B\.1\.351',
    'Gamma':  r'P\.1.*',
    'Delta':  r'(B\.1\.617\.2)|(AY\.[1-9].*)',
    'Lambda': r'C\.37',
    'Mu':     r'B\.1\.621(\.1)?',
    'Omicron': r'(B\.1\.1\.529)|(BA\.[1-9].*)',
}
def expand_who_name_to_pangolin_pattern(patt):
    '''
    if patt is one of the WHO names, then
    replace it with its associated pangolin
    pattern
    '''
    return xpand_WHO_Pangolin.get(patt,patt)

def get_isl(fullname):
    '''return EPI_ISL number from the sequence name'''
    epi_match = EPI_ISL_REGEX.search(fullname)
    return epi_match[0] if epi_match else "X"

def get_lineage_from_name(name):
    '''get pango lineage by parsing the sequence name'''
    m = LINEAGE_REGEX.search(name)
    linpatt = m[1] if m else None
    #try:
    #    tokens = name.split('.',6)
    #    lintok = tokens[6]
    #except IndexError:
    #    lintok = None
    #if linpatt != lintok:
    #    print("name=",name,"patt:",linpatt,"token:",lintok)
    return linpatt

def date_fromiso(s):
    '''return datetime.date object from date string in yyyy-mm-dd format'''
    ## if "." or invalid (quite different cases!), return None
    ## problem with raising error is that many badly formatted dates out there
    if not s:
        return None
    if isinstance(s, datetime.date):
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError:
        if s == ".":
            return None
        return None #raise RuntimeError(f"Invalid Date {s}")

def date_from_seqname(sname):
    '''extract date string from sequence name'''
    try:
        tokens = sname.split('.')
        datestr = tokens[4]
    except IndexError:
        v.vprint_only(5,"Invalid name:",sname)
        datestr = sname
    if not ISO_DATE_REGEX.match(datestr):
        v.vprint_only(5,"Invalid date:",datestr,sname)
        m = ISO_DATE_REGEX_DOTS.search(sname)
        datestr = m[1] if m else None
    return date_fromiso(datestr)

def count_bad_dates(seqlist):
    return sum(date_from_seqname(s.name) is None for s in seqlist)

def range_of_dates(seqlist):
    '''return tuple of iso-formatted dates'''
    assert isinstance(seqlist,list)
    dates = [date_from_seqname(s.name) for s in seqlist]
    dates = [d for d in dates if d is not None]
    v.vprint("Range of dates bsed on",len(dates),"sequences")
    return (min(dates).isoformat(),
            max(dates).isoformat())

def filter_by_date(seqs,fromdate,todate,keepfirst=False):
    '''input seqs is iterable (list or iterator); output is generator'''

    f_date = date_fromiso(fromdate)
    t_date = date_fromiso(todate)

    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        d = date_from_seqname(s.name)
        if not d:
            continue
        if f_date and f_date > d:
            continue
        if t_date and t_date < d:
            continue
        yield s

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
ABBREV_CONTINENTS = {"United-Kingdom"             : "UK",
                     "Europe-minus-United-Kingdom": "Eu-UK",
                     "North-America"              : "NAmer",
                     "Asia"                       : "Asia",
                     "Africa"                     : "Africa",
                     "South-America"              : "SAmer",
                     "Oceania"                    : "Ocean",
}

def parse_continents(withglobal=False):
    '''
    returns a list of three-element tuples (cx,c,x)
    cx: full name of region, possibly including '-minus-'
    c: included part of region, before the '-minus-'
    x: excluded part of region, after the '-minus-'
    [in most cases, cx=c and x=None]
    '''
    cx_c_x=[]
    if withglobal:
        cx_c_x.append(("Global","Global",None))
    for cx in CONTINENTS:
        if "-minus-" in cx:
            c,x = cx.split("-minus-")
        else:
            c,x = cx,None
        cx_c_x.append((cx,c,x))
    return cx_c_x

def filename_prepend(pre,file):
    '''prepend a string to a file name; eg
       "pre","file" -> "prefile", but also
       "pre","directory/file" -> "directory/prefile"
    '''
    ## alt: re.sub(r"(.*/)?([^/]+)",r"\1"+pre+r"\2",file)
    if not file:
        return file
    directory,base = os.path.split(file)
    return os.path.join(directory,pre+base)

def get_title(args):
    '''produce default title for plots and tables'''
    if args.title:
        return "Global" if args.title=='.' else args.title
    if args.filterbyname:
        newfilterbynames = [re.sub('-minus-',' w/o ',name)
                            for name in args.filterbyname]
        title = "+".join(newfilterbynames)
    else:
        title = "Global"
    if args.xfilterbyname:
        title = title + " w/o " + "+".join(args.xfilterbyname)
    if len(title) > MAX_TITLE_LENGTH:
        title = title[:MAX_TITLE_LENGTH-3]+"..."
    return title

def get_first_item(items,keepfirst=True):
    '''
    get first item in iterable, and and put it back;
    works when the iterable is a list or an iterator
    if keepfirst==False, then don't put it back
    '''
    warnings.warn("use sequtil.get not covid.get")
    return sequtil.get_first_item(items,keepfirst=keepfirst)

def read_filter_seqfile(args,**kwargs):
    '''
    read sequence file from args.input,
    and filter according to args,
    return a generator of sequences
    '''
    seqs = read_seqfile(args,**kwargs)
    seqs = filter_seqs(seqs,args)

    return seqs

def read_seqfile(args,**kwargs):
    '''
    read sequences file from args.input
    return a generator of sequences
    '''
    seqs = readseq.read_seqfile(args.input,badchar='X',**kwargs)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences read:")
    return seqs


def fix_seqs(seqs,args):
    '''
    seqs = fix_seqs(seqs,args)
    will return sequences with stripdashcols (obsolete)
    and with last character (if it's $) stripped
    and
    '''
    ## ...a bit heavy-handed, and assumes first seq is still there

    ## we peek at first sequence, but do not remove it from seqs
    first,seqs = sequtil.get_first_item(seqs)

    if "-" in first.seq and args.stripdashcols:
        seqs = sequtil.stripdashcols(first.seq,seqs)

    if not args.keeplastchar and "$" in first.seq:
        first.seq = re.sub(r'\$[^\$]*$','',first.seq)
        seqs = striplastchars(seqs,len(first.seq))

    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)

    return seqs

def striplastchars(seqs,seqlen):
    '''
    truncate all sequences to length seqlen;
    replace any remaining '$' with 'X'
    '''
    for s in seqs:
        s.seq = s.seq[:seqlen]
        s.seq = re.sub(r'\$','X',s.seq)
        yield s

def filter_seqs(seqs,args):
    '''filter sequences according to args: by date, by pattern, by nseq'''
    ## by date first so multiple runs will have the same date range
    ## with --days option
    seqs = filter_seqs_by_date(seqs,args)
    seqs = filter_seqs_by_pattern(seqs,args)
    seqs = fix_seqs(seqs,args)
    if args.nseq:
        seqs = it.islice(seqs,args.nseq+1)
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
    '''return the last date in the range of dates'''
    ## to get last date
    ## 1/ get modification date of input file
    ## 2/ if that doesn't work (eg, if file not found), get today's date
    ## 3/ if that doesn't work (but why wouldn't it?), get last date in dataset
    ## 4/ if that doesn't work, raise RuntimeError
    ## Note, if seqs is None, then skip step 3

    lastdate=None ## in case nothing works!
    try:
        mtime = os.path.getmtime(file)
        lastdate = datetime.date.fromtimestamp(mtime).isoformat()
    except FileNotFoundError:
        try:
            lastdate = datetime.date.today().isoformat()
        except: ## don't think this will ever happen
            if seqs:
                seqs = list(seqs)
                _,lastdate = range_of_dates(seqs)

    if lastdate is None:
        raise RuntimeError('Cannot find last date for date range')

    return lastdate

def date_range_from_args(args):
    '''return list of two iso-formatted dates (start and stop);
    obtain what /should/ be in args.dates, but if it is not set,
    then infer what it should be from args.days'''
    if args.dates:
        return [(None if date=='.' else date)
                for date in args.dates]
    if args.days:
        lastdate = lastdate_byfile(args.input)
        t = date_fromiso(lastdate)
        f = t - datetime.timedelta(days=args.days)
        return [f.isoformat(),t.isoformat()]
    return [None,None]

def expand_date_range(daterange,daysperweek=7):
    '''pad the date range by a few (or many) days out front
    so that weekly and/or cumulative counts are mainteined
    correctly'''
    start_date,stop_date = daterange
    if start_date is None or start_date == '.':
        return daterange
    if daysperweek == 0:
        start_date = None
    elif daysperweek > 1:
        start_date = date_fromiso(start_date)
        start_date = start_date - datetime.timedelta(days=daysperweek-1)
        start_date = start_date.isoformat()
    return [start_date,stop_date]


def filter_seqs_by_date(seqs,args,keepfirst=True):
    '''passes through seq's whose date is in range specified by args;
    also, ensures that args.dates is set to range of dates (eg, if range
    specified by --days, then set args.dates to be consistent with that)
    '''

    if not args.days and not args.dates:
        return seqs
    if args.days and args.dates:
        raise RuntimeError("Cannot specify both --days AND --dates")

    if args.days:
        args.dates = date_range_from_args(args)

    if args.dates:
        seqs = filter_by_date(seqs,args.dates[0],args.dates[1],keepfirst=keepfirst)

    return seqs

def filter_seqs_by_pattern(seqs,args,keepfirst=True):
    '''input is iterable seqs; output is generator seqs'''

    keepers = []
    xcludes = []
    if args.filterbyname:
        for name in args.filterbyname:
            patt,_,xpat = name.partition("-minus-")
            patt = expand_who_name_to_pangolin_pattern(patt)
            if patt != "Global":
                keepers.append(patt)
            if xpat:
                xcludes.append(xpat)
    if args.xfilterbyname:
        for name in args.xfilterbyname:
            name = expand_who_name_to_pangolin_pattern(name)
            xcludes.append(name)

    ## Use r"\."+ to ensure that names have to be preceded by a dot
    keepers = [r"\."+patt for patt in keepers]
    xcludes = [r"\."+xpat for xpat in xcludes]

    if keepers:
        seqs = sequtil.filter_by_patternlist(seqs,keepers,
                                             keepfirst=keepfirst)
    if xcludes:
        seqs = sequtil.filter_by_patternlist_exclude(seqs,xcludes,
                                                     keepfirst=keepfirst)

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
