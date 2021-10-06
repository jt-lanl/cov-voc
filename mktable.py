'''
Make table/spreadsheet with various counts
of pango forms and associated mutant strings

Input is two columns: pango lineage name, and m-string
(it's ok if there's spaces in the m-string)
'''
import sys
import re
import itertools as it
from collections import Counter
import argparse
import warnings

import sequtil
import covid
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(argparser)
    paa = argparser.add_argument
    paa("--pmfile","--pm",
        help="file with pango lineages and mstrings")
    paa("--mutant","-m",
        help="mutant mstring; eg, [Q414K,D614G,T716I]")
    paa("--pango","-p",
        help="name of Pango lineage; eg, B.1.1.7")
    paa("--nopango",action="store_true",
        help="do no pango lineage computation")
    paa("-N",type=int,default=0,
        help="for debugging purposes, only look at first N sequences in datafile")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

## ColumnNames
Total = 'Total'
Pango = 'Pango'
Pattern = 'Pattern'
Count = 'Count'
Fraction = 'Fraction'
TotalPangoCount = 'TotalPangoCount'
TotalPangoFraction = 'TotalPangoFraction'
TotalPatternFullCount = 'TotalPatternFullCount'
TotalPatternFullFraction = 'TotalPatternFullFraction'
TotalPatternInclusiveCount = 'TotalPatternInclusiveCount'
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

ColumnHeaders = {
    Pango: "Pango lineage",
    Pattern: "Spike backbone for variant...",
    TotalPatternFullCount: "Number of sequences thta exactly match this pattern",
    TotalPatternInclusiveCount: "Number of sequences that contain this pattern",
    'LineagesMatchPatternFull': "Lineages with sequences that " \
                                "exactly match this pattern in Spike (and count)",
    'LineagesMatchPatternInclusive': "Lineages with sequences that " \
                                "contain this pattern in Spike (and count)",
    TotalPangoCount: "Total number of sequences with this Pango lineage designation",
    PangoPatternFullCount: "Number of sequences in the Pango lineage that " \
                                "exactly match this pattern",
    PangoPatternFullFraction: "Fraction of sequences in the Pango lineage that " \
                                "exactly match this pattern",
    PangoPatternInclusiveCount: "Number of sequences in the Pango lineage that " \
                                "contain this pattern",
    PangoPatternInclusiveFraction: "Fraction of sequences in the Pango lineage that " \
                                "contain this pattern",
    'CountriesPatternFull': "Countries that exactly match this pattern",
    'CountriesPatternInclusive': "Countries that contain this pattern",    
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

def mstring_seqs(seqs,MM,mstring,exact=False):
    '''return an iterator of seqs that match the mstring pattern'''
    mpatt = mutant.Mutation(mstring)
    vvprint("mpatt:",mpatt)
    return MM.filter_seqs_by_pattern(mpatt,seqs,exact=exact)

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
    2. make a new mstring that treats G142x as G142.
    '''
    ### could we do D215A,+215AGY  => +214AAG,D215Y
    ## fix G-- vs --G
    if re.search('E156G,F157-,R158-',mstring):
        warnings.warn(f'fixing mstring with G-- pattern: {mstring}')
        mstring = re.sub('E156G,F157-,R158-','E156-,F157-,R158G',mstring)

    ## fix G142G or G142_ or G142D -> G142.
    newmut = []
    for ssm in mutant.Mutation(mstring):
        if ssm.site == 142 and ssm.mut in "GD_":
            newmut.append(mutant.SingleSiteMutation.from_ref_site_mut('G',142,'.'))
        else:
            newmut.append(ssm)
    return str(mutant.Mutation(newmut))

def read_input_file(filename):
    items=[]
    with open(filename) as f_input:
        for line in f_input:
            line = re.sub(r'#.*','',line)
            line = line.strip()
            if not line:
                continue
            try:
                pango,mstring = line.split(None,1)
                mstring = mstring_brackets(mstring)
            except ValueError:
                warnings.warn("Problem reading line:",line)
                continue
            items.append((pango,mstring))
    return items


def get_row(seqs,MM,pango,mstring):
    '''produce a single row of the spreadsheet as a dict()'''
    
    row = dict()

    row[Pango] = pango
    row[Pattern] = re.sub(r'[\[\]]','',mstring) ## brackets off

    seqdict=dict()
    seqdict[Total] = seqs
    TotalCount = len(seqs)
    vprint("Read",len(seqs),"sequences")
    seqdict[Pango] = list(pango_seqs(seqs,pango))
    vprint("Of which,",len(seqdict[Pango]),"sequences matched pango form",pango)
    row[TotalPangoCount] = len(seqdict[Pango])
    row[TotalPangoFraction] = row[TotalPangoCount]/TotalCount

    for seqtype in [Total,Pango]:
        seqs = seqdict[seqtype]
        for exact in [False,True]:
            column_name = "PatternFull" if exact else "PatternInclusive"
            column_name = seqtype + column_name
            mstring_adj = mstring_fix(mstring)
            if mstring_adj != mstring:
                print(f"M-String adjusted: {mstring} -> {mstring_adj}",file=sys.stderr)
            matched_seqs = list(mstring_seqs(seqs,MM,mstring_adj,exact=exact))
            vprint(seqtype,exact,len(matched_seqs),len(seqs))

            row[column_name + Count] = len(matched_seqs)
            row[column_name + Fraction] = len(matched_seqs)/len(seqs) if len(seqs) else 0

            if seqtype == Total:
                lineages = [get_lineage_from_name(s.name) for s in matched_seqs]
                cnt_lineages = Counter(lineages)
                sorted_lineages = sorted(cnt_lineages,key=cnt_lineages.get,reverse=True)
                str_lineages = ",".join("%s(%d)" % (lin,cnt_lineages[lin])
                                  for lin in sorted_lineages)
                column_name = "PatternFull" if exact else "PatternInclusive"
                row["LineagesMatch" + column_name] = str_lineages

                countries = [get_country_from_name(s.name) for s in matched_seqs]
                cnt_countries = Counter(countries)
                sorted_countries = sorted(cnt_countries,key=cnt_countries.get,reverse=True)
                str_countries = ",".join("%s(%d)" % (c,cnt_countries[c])
                                         for c in sorted_countries[:3])
                row["Countries" + column_name] = str_countries

    return row

def format_row(row,header=False,sepchar='\t'):
    '''convert row into a string that is a tab-separated list'''
    ##comma-separated doesn't work for m-string patterns
    column_name_list = list(ColumnHeaders)
    if header:
        return sepchar.join(map(column_header,column_name_list))
    return sepchar.join(str(row[column_name])
                        for column_name in column_name_list)


def _main(args):
    '''mktable main'''

    seqs = covid.read_seqfile(args)
    seqs = covid.filter_seqs(seqs,args)
    seqs = sequtil.checkseqlengths(seqs)

    if args.N:
        ## just grab the first N (for debugging w/ shorter runs)
        seqs = it.islice(seqs,args.N+1)

    first,seqs = sequtil.get_first_item(seqs)
    firstseq = first.seq

    MM = mutant.MutationManager(firstseq)

    seqs = list(seqs)

    if args.pango and args.mutant:
        row = get_row(seqs,MM,args.pango,args.mutant)
        print(format_row(row,header=True))
        print(format_row(row))

    if args.pmfile:
        header_yet=False
        for pango,mstring in read_input_file(args.pmfile):
            if args.nopango:
                pango=""
            row = get_row(seqs,MM,pango,mstring)
            if not header_yet:
                print(format_row(row,header=True))
                header_yet=True
            print(format_row(row))
            print(end="",flush=True)

if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)

    _main(_args)
