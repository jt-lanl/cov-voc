'''
Make table/spreadsheet with various counts
of pango forms and associated mutant strings

Input is two columns: pango lineage name, and m-string
(it's ok if there's spaces in the m-string)
'''

## parallel usage: (48 variants per process; ie, 48 lines of variantes-in.tsv)
## parallel -k -N48 --pipe \
##  pm mktable --jobno {#} -M - < variants-in.tsv > variants-out.tsv

import sys
import re
from collections import Counter
import datetime
import argparse
import warnings

import verbose as v
import intlist
import breakpipe
import sequtil
import covid
import mutant
import pseq
import mstringfix

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(argparser)
    paa = argparser.add_argument
    paa("--pmfile","-M",
        help="file with pango lineages and mstrings")
    paa("--rows","-R",
        help="(0-based int)List of rows to compute")
    paa("--mutant","-m",
        help="mutant mstring; eg, [Q414K,D614G,T716I]")
    paa("--pango","-p",
        help="name of Pango lineage; eg, B.1.1.7")
    paa("--nopango",action="store_true",
        help="do no pango lineage computation")
    paa("--tweakfile",
        help="Use tweakfile to fix mstrings")
    paa("--jobno",type=int,default=1,
        help="Use --jobno={#} in parallel mode")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)

    if args.days:
        warnings.warn("Probably you do not want to use --days; "
                      "both full and days=60 are computed")

    return args

## ColumnNames
Total = 'Total'
Pango = 'Pango'
Pattern = 'Pattern'
Count = 'Count'
Fraction = 'Fraction'
Recent = 'SixtyDays'
TotalPangoCount = 'TotalPangoCount'
TotalPangoFraction = 'TotalPangoFraction'
TotalPatternFullCount = 'TotalPatternFullCount'
SixtyDaysPatternFullCount = 'SixtyDaysPatternFullCount'
TotalPatternFullFraction = 'TotalPatternFullFraction'
TotalPatternInclusiveCount = 'TotalPatternInclusiveCount'
SixtyDaysPatternInclusiveCount = 'SixtyDaysPatternInclusiveCount'
TotalPatternInclusiveFraction = 'TotalPatternInclusiveFraction'
PangoPatternFullCount = 'PangoPatternFullCount'
PangoPatternFullFraction = 'PangoPatternFullFraction'
PangoPatternInclusiveCount = 'PangoPatternInclusiveCount'
PangoPatternInclusiveFraction = 'PangoPatternInclusiveFraction'
CommonPango1 = 'CommonPango1'
CommonPango1Count = 'CommonPango1Count'
CommonPango2 = 'CommonPango2'
CommonPango2Count = 'CommonPango2Count'
CommonPango3 = 'CommonPango3'
CommonPango3Count = 'CommonPango3Count'
ExampleISL = 'ExampleISL'

ColumnHeaders = {
    Pango: "Pango lineage",
    Pattern: "Spike backbone for variant...",
    TotalPatternFullCount: "Number of sequences that exactly match this pattern",
    TotalPatternInclusiveCount: "Number of sequences that contain this pattern",
    SixtyDaysPatternFullCount: "Number of sequences that exactly match this pattern, last 60 days",
    SixtyDaysPatternInclusiveCount: "Number of sequences that contain this pattern, last 60 days",
    'TotalLineagesMatchPatternFull': "Lineages with sequences that " \
                                "exactly match this pattern in Spike (and count)",
    'TotalLineagesMatchPatternInclusive': "Lineages with sequences that " \
                                "contain this pattern in Spike (and count)",
    'SixtyDaysLineagesMatchPatternFull': 'Lineages with sequences that exactly match pattern, last sixty days',
    'SixtyDaysLineagesMatchPatternInclusive': 'Lineages with sequenecs that contain pattern, last sixty days',
    TotalPangoCount: "Total number of sequences with this Pango lineage designation",
    PangoPatternFullCount: "Number of sequences in the Pango lineage that " \
                                "exactly match this pattern",
    PangoPatternFullFraction: "Fraction of sequences in the Pango lineage that " \
                                "exactly match this pattern",
    PangoPatternInclusiveCount: "Number of sequences in the Pango lineage that " \
                                "contain this pattern",
    PangoPatternInclusiveFraction: "Fraction of sequences in the Pango lineage that " \
                                "contain this pattern",
    'CountriesTotalPatternFull': "Countries that exactly match this pattern",
    'CountriesTotalPatternInclusive': "Countries that contain this pattern",
    'CountriesSixtyDaysPatternFull': "Countries that exactly match this pattern, "\
    "counts based on last 60 days",
    'CountriesSixtyDaysPatternInclusive': "Countries that contain this pattern, "\
    "counts based on last 60 days",
    ExampleISL: 'Example ISL number for this pattern',
}

def column_header(column_name):
    '''fuller column header given shorter column name'''
    ## if not in ColumnHeaders dict, then just return the name
    return ColumnHeaders.get(column_name,column_name)

def get_lineage_from_name(name):
    '''get pango lineage by parsing the sequence name'''
    return re.sub(r".*EPI_ISL_\d+\.","",name)

def get_country_from_name(name):
    '''get name of country from sequence name'''
    return get_region_from_name(2,name)

def get_region_from_name(level,name):
    '''return the name of the region (based on level) in the full sequence name'''
    m = name.split(".")
    return m[level] if m else None

def pango_seqs(seqs,pango):
    '''return an iterator of seqs whose names indicate the pango type'''
    if not pango:
        yield from []
    else:
        for s in seqs:
            if pango == get_lineage_from_name(s.name):
                yield s

def pango_pseqs(pseqs,pango):
    '''iterator of pseqs that are consistent with pango'''
    if not pango:
        yield from []
    else:
        pango = covid.expand_who_name_to_pangolin_pattern(pango)
        for ps in pseqs:
            try:
                if re.match(pango,ps.lineage):
                    yield ps
            except:
                v.print("pango=",pango)
                v.print("ps.lineage=",ps.lineage)
                raise RuntimeError(f"Invalid pango name: {pango}/{ps.lineage}")

def sixtydays_seqs(seqs,days=60,file=None):
    '''return an iterator of seqs whose dates are in the last 60 days'''
    lastdate = covid.lastdate_byfile(file,seqs)
    t = covid.date_fromiso(lastdate)
    f = t - datetime.timedelta(days=days)
    v.vprint("Sixty days:",f,t)
    return covid.filter_by_date(seqs,f,t,firstisref=False)

def sixtydays_pseqs(pseqs,days=60,file=None):
    '''return an iterator of seqs whose dates are in the last 60 days'''
    lastdate = covid.lastdate_byfile(file,pseqs)
    t = covid.date_fromiso(lastdate)
    f = t - datetime.timedelta(days=days)
    v.vprint("Sixty days:",f,t)
    for ps in pseqs:
        if ps.date and ps.date > f:
            yield ps

def mstring_seqs(seqs,m_mgr,mstring,exact=False):
    '''return an iterator of seqs that match the mstring pattern'''
    mpatt = mutant.Mutation(mstring)
    v.vvprint("mpatt:",mpatt)
    return m_mgr.filter_seqs_by_pattern(mpatt,seqs,exact=exact)

def mstring_pseqs(pseqs,m_mgr,mstring,exact=False):
    '''return an iterator of seqs that match the mstring pattern'''
    mpatt = mutant.Mutation(mstring)
    v.vvprint("mpatt:",mpatt)
    ## if exact==True, then assume we just did inclusive (back when exact was False)
    ## this helps a little, but not that much
    return pseq.filter_pseqs_by_pattern(m_mgr,mpatt,pseqs,
                                        exact=exact,assume_inclusive=exact)

def xopen(filename):
    if filename == '-':
        return sys.stdin
    return open(filename)

def read_input_file(filename):
    '''return list of (pango,mstring) pairs'''
    with xopen(filename) as f_input:
        for line in f_input:
            line = re.sub(r'#.*','',line)
            line = line.strip()
            if not line:
                continue
            try:
                pango,mstring = line.split('\t',1)
                pango = pango.strip()
                pango = re.sub(r' ','',pango)
                mstring = mstringfix.mstring_brackets(mstring)
            except ValueError:
                warnings.warn(f'Problem reading line: {line}')
                continue
            yield (pango,mstring)

def truncate_pango_name(pangofull):
    '''
    return the pango lineage associated ith a full pango name
    it is the SECOND field, as delineated by underscores; eg
    Iota_B.1.526 -> B.1.526
    Omicron_BA.1.1__BA.1_add[R346K] -> BA.1.1

    Note: currently, we treat 'xxx_Omicron_BA.1' as an invalid string,
    because 'Omicron' is not a strict pango lineage name

    Note also: the pango string is assumed to be of the form
        [A-Z]+(\.\d+)+
    ie, alphabetical.numerical.numerical...numrical
    with an indeterminate number (possibly even zero)
    of '.numerical' strings
    '''
    match = re.match(r'[^_]+_+([A-Z]+(\.\d+)*)',pangofull)
    if not match:
        raise RuntimeError(f'Invalid designation: {pangofull}')
    pango = match[1].strip()
    return pango

def get_row(seqs,seqs_sixtydays,m_mgr,pangofull,mstring):
    '''produce a single row of the spreadsheet as a dict()'''

    row = dict()

    row[Pango] = pangofull ## give it the full name
    pango = truncate_pango_name(pangofull)

    row[Pattern] = re.sub(r'[\[\]]','',mstring) ## brackets off

    seqdict=dict()
    seqdict[Total] = seqs
    total_count = len(seqs)
    v.vvprint("Read",len(seqs),"sequences")
    seqdict[Pango] = list(pango_pseqs(seqs,pango))
    v.vvprint("Of which,",len(seqdict[Pango]),
             "sequences matched pango form",pangofull)
    if len(seqdict[Pango])==0:
        v.print(f'No matches to pango form {pango}, full name {pangofull}')

    row[TotalPangoCount] = len(seqdict[Pango])
    row[TotalPangoFraction] = row[TotalPangoCount]/total_count
    seqdict[Recent] = seqs_sixtydays

    for seqtype in [Total,Pango,Recent]:
        seqs = seqdict[seqtype]

        matched_seqs = seqs
        for exact in [False,True]:  ## do False first, then use True to further filter

            column_name = "PatternFull" if exact else "PatternInclusive"
            column_name = seqtype + column_name
            matched_seqs = list(mstring_pseqs(matched_seqs,m_mgr,mstring,exact=exact))
            v.vvprint(pango,seqtype,f'exact={exact}',len(matched_seqs),len(seqs))

            row[column_name + Count] = len(matched_seqs)
            row[column_name + Fraction] = len(matched_seqs)/len(seqs) if len(seqs) else 0

            if seqtype == Total and exact is True:
                if len(matched_seqs):
                    #exampleseq = matched_seqs[-1]
                    #row[ExampleISL] = exampleseq.ISL
                    row[ExampleISL] = pseq.lowest_isl(matched_seqs)
                else:
                    row[ExampleISL] = "NA"
                    warnings.warn(f"For pango={pangofull}, no sequences "
                                  f"exactly match: {mstring}")

            if seqtype in [Total,Recent]:
                lineages = [ps.lineage for ps in matched_seqs]
                cnt_lineages = Counter(lineages)
                sorted_lineages = sorted(cnt_lineages,key=cnt_lineages.get,reverse=True)
                str_lineages = ",".join("%s(%d)" % (lin,cnt_lineages[lin])
                                  for lin in sorted_lineages)
                column_name = "PatternFull" if exact else "PatternInclusive"
                row[seqtype + "LineagesMatch" + column_name] = str_lineages

            if seqtype in [Total,Recent]:
                countries = [get_country_from_name(s.name) for s in matched_seqs]
                cnt_countries = Counter(countries)
                sorted_countries = sorted(cnt_countries,key=cnt_countries.get,reverse=True)
                str_countries = ",".join("%s(%d)" % (c,cnt_countries[c])
                                         for c in sorted_countries[:3])
                row["Countries" + seqtype + column_name] = str_countries

    v.vvprint(list(row))
    return row

def format_row(row=None,sepchar='\t'):
    '''
    convert row into a string that is a tab-separated list;
    if row is None, then print the header
    '''
    ##comma-separated doesn't work for m-string patterns
    column_name_list = list(ColumnHeaders)
    if row is None:
        return sepchar.join(map(column_header,column_name_list))
    return sepchar.join(str(row[column_name])
                        for column_name in column_name_list)


@breakpipe.no_broken_pipe
def _main(args):
    '''mktable main'''

    plist=[]
    mlist=[]
    if args.pango and args.mutant:
        plist.append(args.pango)
        mlist.append(args.mutant)

    if args.pmfile:
        for pango,mstring in read_input_file(args.pmfile):
            if args.nopango:
                pango=""
            v.vprint(pango,'\t',mstring)
            plist.append(pango)
            mlist.append(mstring)

    ## If at this point mlist is empty, then something is wrong!
    if not mlist:
        raise RuntimeError('No mstrings are input: check -M or (-m and -p)')

    ## Do this now just to make sure we don't get an error
    pango_truncated = [truncate_pango_name(pango) for pango in plist]

    fixer = mstringfix.MStringFixer(args.tweakfile)
    fixer.append(r'\s*ancestral\s*','')
    fixer.append(r'G142[GD_]','G142.')
    v.vprint(f"FIXER:\n{fixer}")
    mlist = [fixer.fix(mstring) for mstring in mlist]

    if args.rows:
        ndxlist = intlist.string_to_intlist(args.rows)
        plist = [plist[ndx] for ndx in ndxlist]
        mlist = [mlist[ndx] for ndx in ndxlist]

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    m_mgr = mutant.MutationManager(first.seq)

    seqs = [pseq.ProcessedSequence(m_mgr,s) for s in seqs]

    v.vprint("Sequences processed:",len(seqs))

    seqs_sixtydays = sixtydays_seqs(seqs,file=args.input)
    seqs_sixtydays = list(seqs_sixtydays)
    v.vprint("Read",len(seqs),"sequences")
    v.vprint("Of which",len(seqs_sixtydays),"are from the last 60 days")

    if args.jobno == 1:
        print(format_row(None)) ## print the header
    for pango,mstring in zip(plist,mlist):
        row = get_row(seqs,seqs_sixtydays,m_mgr,pango,mstring)
        print(format_row(row),flush=True)
        #print(end="",flush=True) ## what does this do? just flush?

    if args.jobno == 1:
        v.print("Having made table, you may want to run 'mutisl' to get fasta file with DNA sequences")

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
