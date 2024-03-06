'''covid-specific utilities and hardcoded values'''

import os
import re
import argparse
import datetime
import itertools as it
from pathlib import Path
import warnings

import verbose as v
import intlist
import wrapgen
import readseq
import sequtil
import mstringfix
import mutant

MAX_TITLE_LENGTH=60 ## truncate long title names

DEFAULTSEQFILE_KEEPX="Latest-keepx.fasta"
DEFAULTSEQFILE_SKIPX="Latest-skipx.fasta"
REF_SEQUENCE_NAME="NC_045512" #first part of name for all proteins


def find_seqfile(seqfilename,skipx=False):
    '''
    return the default file for input sequences;
    hunt around in various directories until you find it
    '''
    if not seqfilename or str(seqfilename) == ".":
        seqfilename = (DEFAULTSEQFILE_SKIPX if skipx
                       else DEFAULTSEQFILE_KEEPX)

    for special in ['-', '/dev/stdin']:
        if str(seqfilename) == special:
            return seqfilename

    for directory in ['.',os.getenv('DATA'),'data','..','../data',]:
        if not directory:
            continue
        for ext in ("",".zst", ".gz",".xz"):
            seqfile = Path(directory) / Path(str(seqfilename) + ext)
            if seqfile.exists():
                return seqfile

    v.print(f'Input sequence file "{seqfilename}" not found')
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

    paa = ap.add_argument_group('Input Options').add_argument
    paa("--input","-i",type=Path,
        default=None,
        help="input file with aligned sequences (first is reference)")
    paa("--fastread",action="store_true",
        help="assume input file is already fixed/cleaned-up")
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
    paa("--keeplastchar",action="store_true",
        help="Do not strip final stop codon from end of sequences")
    paa("--keepx",action="store_true",default=False,
        help="Keep sequences that include bad characters, denoted X")
    paa("--skipx",action="store_false",dest='keepx',
        help="Skip sequences that include bad characters, denoted X")
    ap.set_defaults(keepx=False)

def corona_fixargs(args):
    '''after argparser.parse_args() has been called,
    then call this for cleanup'''
    if not args.input or str(args.input) == '.':
        args.fastread = True
    args.input = find_seqfile(args.input,
                              skipx=(not args.keepx))
    return args


## Routine for post-arg manipulation
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


#### Routines for dealing with lineages, expand name, baselin seqs
#### (could some of this be merged with lineagenotes.py ??)

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

BASELINE_MSTRINGS = {
    'Wuhan' : "",
    'BA.2' : "T19I,L24-,P25-,P26-,A27S,G142D,V213G,G339D,S371F,S373P,S375F,T376A,D405N,R408S,K417N,N440K,S477N,T478K,E484A,Q493R,Q498R,N501Y,Y505H,D614G,H655Y,N679K,P681H,N764K,D796Y,Q954H,N969K",
    'BA.5' : "T19I,L24-,P25-,P26-,A27S,H69-,V70-,G142D,V213G,G339D,S371F,S373P,S375F,T376A,D405N,R408S,K417N,N440K,L452R,S477N,T478K,E484A,F486V,Q498R,N501Y,Y505H,D614G,H655Y,N679K,P681H,N764K,D796Y,Q954H,N969K",
    'BA.2.75': "T19I,L24-,P25-,P26-,A27S,G142D,K147E,W152R,F157L,I210V,V213G,G257S,G339H,S371F,S373P,S375F,T376A,D405N,R408S,K417N,N440K,G446S,N460K,S477N,T478K,E484A,Q498R,N501Y,Y505H,D614G,H655Y,N679K,P681H,N764K,D796Y,Q954H,N969K",
    'XBB.1.5': "T19I,L24-,P25-,P26-,A27S,V83A,G142D,Y144-,H146Q,Q183E,V213E,G252V,G339H,R346T,L368I,S371F,S373P,S375F,T376A,D405N,R408S,K417N,N440K,V445P,G446S,N460K,S477N,T478K,E484A,F486P,F490S,Q498R,N501Y,Y505H,D614G,H655Y,N679K,P681H,N764K,D796Y,Q954H,N969K",
    'BA.2.86': "+16MPLF,T19I,R21T,L24-,P25-,P26-,A27S,S50L,H69-,V70-,V127F,G142D,Y144-,F157S,R158G,N211-,L212I,V213G,L216F,H245N,A264D,I332V,G339H,K356T,S371F,S373P,S375F,T376A,R403K,D405N,R408S,K417N,N440K,V445H,G446S,N450D,L452W,N460K,S477N,T478K,N481K,V483-,E484K,F486P,Q498R,N501Y,Y505H,E554K,A570V,D614G,P621S,H655Y,N679K,P681R,N764K,D796Y,S939F,Q954H,N969K,P1143L",
    'JN.1': "+16MPLF,T19I,R21T,L24-,P25-,P26-,A27S,S50L,H69-,V70-,V127F,G142D,Y144-,F157S,R158G,N211-,L212I,V213G,L216F,H245N,A264D,I332V,G339H,K356T,S371F,S373P,S375F,T376A,R403K,D405N,R408S,K417N,N440K,V445H,G446S,N450D,L452W,L455S,N460K,S477N,T478K,N481K,V483-,E484K,F486P,Q498R,N501Y,Y505H,E554K,A570V,D614G,P621S,H655Y,N679K,P681R,N764K,D796Y,S939F,Q954H,N969K,P1143L",
}

def get_baseline_mstring(lineage='Wuhan'):
    '''return mstring for the given lineage'''
    mstring = BASELINE_MSTRINGS.get(lineage,"")
    mstring = mstringfix.mstring_brackets(mstring)
    return mstring

def reset_baseline(firstseq,lineage):
    '''use lineage instead of Wuhan as baseline "first" sequence'''
    mut_mgr = mutant.MutationManager(firstseq)
    base_mstring = get_baseline_mstring(lineage)
    if not base_mstring:
        return firstseq
    base_mutant = mutant.Mutation.from_mstring(base_mstring)
    base_seq = mut_mgr.seq_from_mutation(base_mutant)
    return base_seq

def ndxlist_from_sites(m_mgr,sites,compact=False,contig=False):
    '''return list of indices associated with sites'''
    ndxlist = []
    for site in sites:
        if compact:
            ndxlist.append(m_mgr.index_from_site(site))
        else:
            ndxlist.extend(m_mgr.indices_from_site(site))
    ndxlist = sorted(ndxlist)
    if contig:
        ndxlist = list(range(min(ndxlist),max(ndxlist)+1))
    return ndxlist

def keepsites(mut_mgr,seqs,sitestring):
    '''replace sequence with shorter sequence including only sites in list'''

    sitelist = intlist.string_to_intlist(sitestring)
    keepndx = ndxlist_from_sites(mut_mgr,sitelist)
    keepranges = intlist.intlist_to_rangelist(keepndx)
    v.vprint('keepranges:',keepranges)
    for s in seqs:
        s.seq = "".join(s.seq[lo:hi] for lo,hi in keepranges)
        yield s

## Routines for parsing sequence names, which for the LANL database look like:
## HDF-IPP45030.Europe.France.Hauts-de-France_Denain.2021-11-22.EPI_ISL_9089887.AY.42
## corresponding to dot-separated components:
## idstring.Location.Location.Location.Date.EPI-number.Lineage
## Note that the locations are (generally): Continent.Country.Region(State,Province,etc)
## And that the Lineage may have further dots (thus, only the first six dots are deliminters).

ISO_DATE_REGEX = re.compile(r'\d\d\d\d-\d\d(-\d\d)?')
ISO_DATE_REGEX_DOTS = re.compile(r'\.(\d\d\d\d-\d\d-\d\d)\.')
EPI_ISL_REGEX = re.compile(r'EPI_ISL_\d+')
LINEAGE_REGEX = re.compile(r'EPI_ISL_\d+\.(.*)')

def parse_seqname(myparser):
    '''decorator function that allows argument to be s or s.name
    where s is a SequenceSample'''
    def wrapper(seq_or_name,*args,**kwargs):
        if isinstance(seq_or_name, sequtil.SequenceSample):
            name = seq_or_name.name
        elif isinstance(seq_or_name, str):
            name = seq_or_name
        else:
            raise TypeError(f'argument seq_or_seqname in function {myparser.__name__} '
                            f'is of type {type(seq_or_name)}, but it needs to be '
	                    'either a string or a SequenceSample')
        return myparser(name,*args,*kwargs)
    return wrapper

def deprecated(usefcn=None):
    '''decrator that indentifies functions as deprecated'''
    def decorator(myfcn):
        def wrapper(*args,**kwargs):
            warning = f'Warning: "{myfcn.__name__}" is deprecated; '
            if usefcn:
                warning += f'use function "{usefcn}" instead.'
            v.print_only(1,warning)
            return myfcn(*args,**kwargs)
        return wrapper
    return decorator


@parse_seqname
def get_isl(seqname,nomatch="X"):
    '''return EPI_ISL number from the sequence name'''
    epi_match = EPI_ISL_REGEX.search(seqname)
    return epi_match[0] if epi_match else nomatch

@parse_seqname
def get_lineage(seqname):
    '''get pango lineage by parsing the sequence name'''
    mat = LINEAGE_REGEX.search(seqname)
    linpatt = mat[1] if mat else None
    #try:
    #    tokens = seqname.split('.',6)
    #    lintok = tokens[6]
    #except IndexError:
    #    lintok = None
    #if linpatt != lintok:
    #    print("name=",seqname,"patt:",linpatt,"token:",lintok)
    return linpatt

@deprecated(usefcn='get_lineage')
def get_lineage_from_name(seqname):
    '''same as get_lineage(), this is deprecated longer name'''
    return get_lineage(seqname)

def date_fromiso(yyyymmdd):
    '''return datetime.date object from date string in yyyy-mm-dd format'''
    ## if "." or invalid (quite different cases!), return None
    ## problem with raising error is that many badly formatted dates out there
    ## this routine is necessary since date.fromisoformat() function is
    ## not available until python version 3.7
    if isinstance(yyyymmdd, datetime.date):
        return yyyymmdd
    try:
        yyyy,mm,dd = yyyymmdd.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except (ValueError,AttributeError,TypeError):
        ## various things can go wrong
        return None

@parse_seqname
def get_date(sname,as_string=False):
    '''extract date string from sequence name'''
    try:
        tokens = sname.split('.')
        datestr = tokens[4]
    except IndexError:
        v.vvprint_only(5,"Invalid name:",sname)
        datestr = sname
    if not ISO_DATE_REGEX.match(datestr):
        v.vvprint_only(5,"Invalid date:",datestr,sname)
        mat = ISO_DATE_REGEX_DOTS.search(sname)
        datestr = mat[1] if mat else None
    if as_string:
        return datestr
    return date_fromiso(datestr)

@deprecated(usefcn='get_date')
def date_from_seqname(seqname,**kwargs):
    '''deprecated: use get_date'''
    return get_date(seqname,**kwargs)

## Parse info in all the seqnames

def count_bad_dates(seqlist):
    '''return a count of how many bad dates in the seqlist'''
    if not hasattr(seqlist, '__len__'):
        raise TypeError(f'argument is of type {type(seqlist)} and should be a list')
    return sum(get_date(s.name) is None for s in seqlist)

def range_of_dates(seqlist):
    '''return tuple of iso-formatted dates'''
    assert isinstance(seqlist,list)
    dates = map(get_date,seqlist)
    dates = [d for d in dates if d is not None]
    v.vprint("Range of dates based on",len(dates),"sequences")
    return (min(dates).isoformat(),
            max(dates).isoformat())

### Routines for filtering sequences (usually based on parsing sequence name)

def test_isref(first):
    '''is this sequence actually the reference sequence?'''
    assert isinstance(first,sequtil.SequenceSample)
    if REF_SEQUENCE_NAME not in first.name:
        return False
    return True

def get_first_item(seqs,keepfirst=None):
    '''wraps the sequtil version, checking if the first really is a ref'''
    kwargs = {'keepfirst': keepfirst} if keepfirst is not None else {}
    first,seqs = sequtil.get_first_item(seqs,**kwargs)
    if keepfirst is False and test_isref(first) is False:
        warnings.warn('first seq is not reference sequence')
    return first,seqs

def get_first_item_ifref(seqs):
    '''return seqs[0],seqs[1:] if seqs[0] is ref sequence;
       return None, seqs       if seqs[0] is not ref sequence'''
    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    if not test_isref(first):
        first = None
    else:
        _,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    return first,seqs

def divert_if_firstisref(myfilter):
    '''
    decorator function adds 'firstisref' to the keywordlist of
    myfilter, then checks for 'firstisref=True' keyword in the call to
    myfilter and if that it the case, then it "diverts" that first
    sequence from the filtering done myb my filter, pulls the first
    sequence out, doing thefiltering, and then putting that unfiltered
    first sequence back in
    '''

    def wrapper(seqs,*args,**kwargs):
        firstisref = kwargs.pop('firstisref',None)
        keepfirst = kwargs.pop('keepfirst',None)
        if keepfirst is not None:
            warnings.warn(f'The keepfirst={keepfirst} keyword is deprecated; '
                          f'use refisfirst={~bool(keepfirst)} instead.')
        if firstisref is False:
            seqs = myfilter(seqs,*args,**kwargs)
        elif firstisref is True:
            first,seqs = get_first_item(seqs,keepfirst=False)
            seqs = myfilter(seqs,*args,**kwargs)
            seqs = it.chain([first],seqs)
        elif firstisref is None:
            first,seqs = get_first_item_ifref(seqs)
            seqs = myfilter(seqs,*args,**kwargs)
            if first:
                seqs = it.chain([first],seqs)
        else:
            raise ValueError(f'keyword firstisref={firstisref} '
                             'should be True, False, or None')
        return seqs
    return wrapper

@divert_if_firstisref
def filter_by_date(seqs,fromdate,todate):
    '''input seqs is iterable (list or iterator); output is generator'''

    f_date = date_fromiso(fromdate) or datetime.date.min
    t_date = date_fromiso(todate) or datetime.date.max

    for s in seqs:
        the_date = get_date(s)
        if not the_date:
            continue
        if f_date <= the_date <= t_date:
            yield s

@divert_if_firstisref
def filter_seqs_by_date(seqs,args):
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
        seqs = filter_by_date(seqs,args.dates[0],args.dates[1])

    return seqs

@divert_if_firstisref
def filter_seqs_by_pattern(seqs,args):
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

    ## Use word boundary to ensure that whole names are matched;
    ## thus patt="India" should not match "Indiana"
    ## Note that std \b is replaced by (?=\b|_)
    ## so that _ can be used as a word boundary
    word_boundary=r"(?=\b|_)"
    keepers = [word_boundary+patt+word_boundary for patt in keepers]
    xcludes = [word_boundary+xpat+word_boundary for xpat in xcludes]

    if keepers:
        seqs = sequtil.filter_by_patternlist(seqs,keepers)
    if xcludes:
        seqs = sequtil.filter_by_patternlist_exclude(seqs,xcludes)

    return seqs

@divert_if_firstisref
def filter_seqs(seqs,args):
    '''filter sequences according to args: by date, by pattern, by nseq'''
    ## by date first so multiple runs will have the same date range
    ## with --days option
    seqs = filter_seqs_by_date(seqs,args)
    seqs = filter_seqs_by_pattern(seqs,args)
    if args.nseq:
        seqs = it.islice(seqs,args.nseq)
    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences filtered:")

    return seqs

def truncate_seqs(seqs,seqlen):
    '''
    truncate all sequences to length seqlen;
    '''
    for s in seqs:
        s.seq = s.seq[:seqlen]
        yield s

def fix_seqs(first,seqs,args):
    '''
    Assumes that seqs still contains the first sequence
    '''

    if not(hasattr(first,'seq') and first.seq):
        ## since .nm file doesn't even have .seq attribute
        return seqs

    if not args.keeplastchar and first.seq[-1] in "$*X":
        seqs = truncate_seqs(seqs,len(first.seq[:-1]))

    return seqs

def read_seqfile(args,**kwargs):
    '''
    read sequences file from args.input
    return a generator of sequences
    '''
    if args.fastread:
        ## bypass all the "fixing" and low-level filtering
        ## can still filter by location, date, etc
        ## assumes input file has already been fixed/cleaned-up
        v.vprint('In --fastread mode, sequences are not filtered or "fixed"')
        if kwargs:
            warnings.warn(f'Ignoring keyword args: {kwargs}')
        seqs = readseq.read_seqfile(args.input,nofilter=True)
        return seqs

    seqs = readseq.read_seqfile(args.input,badchar='X',**kwargs)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    seqs = fix_seqs(first,seqs,args)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences read:")
    return seqs

def read_filter_seqfile(args,firstisref=True,**kwargs):
    '''
    read sequence file from args.input,
    and filter according to args,
    return a generator of sequences
    '''
    seqs = read_seqfile(args,**kwargs)
    seqs = filter_seqs(seqs,args,firstisref=firstisref)

    return seqs

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
        except Exception as err: ## don't think this will ever happen
            v.print(f'Exception should never happen: {err}')
            if seqs:
                seqs = list(seqs)
                _,lastdate = range_of_dates(seqs)

    if lastdate is None:
        raise RuntimeError('Cannot find last date for date range')

    return lastdate

## Filtering seqs, helper function
## Parsing args (used for filter-by-date, coming later
def date_range_from_args(args):
    '''return list of two iso-formatted dates (start and stop);
    obtain what /should/ be in args.dates, but if it is not set,
    then infer what it should be from args.days'''
    if args.dates:
        return [(None if date=='.' else date)
                for date in args.dates]
    if args.days:
        lastdate = lastdate_byfile(args.input)
        tdate = date_fromiso(lastdate)
        fdate = tdate - datetime.timedelta(days=args.days)
        return [fdate.isoformat(),tdate.isoformat()]
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
    returns a list of three-element tuples (cx,c_include,c_exclude)
    cx: full name of region, possibly including '-minus-'
    c_include: included part of region, before the '-minus-'
    c_exclude: excluded part of region, after the '-minus-'
    [in most cases, cx=c_include and c_exclude=None]
    '''
    cx_c_x=[]
    if withglobal:
        cx_c_x.append(("Global","Global",None))
    for cx in CONTINENTS:
        if "-minus-" in cx:
            c_include,c_exclude = cx.split("-minus-")
        else:
            c_include,c_exclude = cx,None
        cx_c_x.append((cx,c_include,c_exclude))
    return cx_c_x


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
    ''' provides a tuple of three lists: nom,beg,end
        nom = name of region
        beg = start site
        end = stop site
    '''
    nom=[]
    beg=[]
    end=[]
    regions = MERS_REGIONS if virus.lower() == "mers" else SARS_REGIONS
    for sr in regions.split('\n'):
        srs = sr.split()
        if len(srs) == 3 and srs[0][0]!='#':
            nom.append(srs[0])
            beg.append(int(srs[1]))
            end.append(int(srs[2]))
    return nom,beg,end

### generic utility (only used by xspike)
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
