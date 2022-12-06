'''routines for dealing with OWID (Our World In Data) database of case counts'''

## The OWID data is available directly at
## https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv
## You can actually use that URL as the "--input ..." argument
## The human-readable site is at
## https://github.com/owid/covid-19-data/tree/master/public/data

import re
import argparse
import datetime
import numpy as np
import pandas as pd

import verbose as v
import embersutil as emu

## Routines for harmonizing raw OWID file

def _getargs():
    '''get command line arguments'''
    parser = argparse.ArgumentParser(description=__doc__)
    paa = parser.add_argument
    paa("--input","-i",required=True,
        help="input csv file, use 'OWID' for direct input of OWID website")
    paa("--output","-o",
        help="output GISAID-harmoinzed case counts csv file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = parser.parse_args()
    return args

STD_URL=\
    "https://raw.githubusercontent.com/owid/covid-19-data/" \
    "master/public/data/owid-covid-data.csv"

ATTRIBUTES=('continent','location','date','total_cases')

OWID_TO_GISAID_DICT = {
    'Bahamas': 'The Bahamas',
    'Congo': 'Republic of the Congo',
    'Democratic Republic of Congo': 'Democratic Republic of the Congo',
    'United States Virgin Islands': 'U.S. Virgin Islands',
    'United States': 'USA',
    'Timor': 'Timor-Leste',
    'Sint Maarten (Dutch part)': 'Sint Maarten',
    'Czechia': 'Czech Republic',
    'Micronesia (country)': 'Micronesia',
    'Cape Verde': 'Cabo Verde',
    'Faeroe Islands': 'Faroe Islands',
    'Samoa': 'American Samoa',
    'Wallis and Futuna': 'Wallis and Futuna Islands',
}

def owid_to_gisaid(name):
    '''convert OWID name to GISAID name'''
    name = OWID_TO_GISAID_DICT.get(name,name)
    if name and not pd.isna(name):
        ## remove dots
        name = re.sub(r'\.','',name)
        ## change spaces and apostrophes into hyphens
        name = re.sub(r"'",'-',name)
        name = re.sub(r' ','-',name)
    return name

def harmonize_dataframe(df_owid):
    '''from input OWID dataframe, return harmonized GISAID dataframe'''
    df = df_owid.copy()
    ## harmonize names
    for locate in ['location','continent']:
        df.loc[:,locate] = list(map(owid_to_gisaid,df[locate]))
    ## set NA's to zero, and make counts integers
    df.loc[:,'total_cases'] = \
        list(map(lambda x: 0 if np.isnan(x) else int(x),
                 df['total_cases']))
    ## only keep the columns of interest
    df = df.loc[:,ATTRIBUTES]
    return df

def read_dataframe(filename):
    '''return a pandas dataframe, with case counts'''
    if filename:
        df = pd.read_csv(filename,sep=',',low_memory=False)
        for att in ATTRIBUTES:
            if att not in df.columns:
                v.print(f'Invalid OWID file {filename} '
                        f'does not contain attribute {att}')
                df = None
                break
    else:
        df = None
    return df

### Routines for filtering/counting cases from harmonized csv file

def filter_cases(df_cases,filterbyname=None,xfilterbyname=None):
    '''
    filter input dataframe with case counts so that
    it only includes the locationa and/or continents desired
    '''
    if df_cases is None:
        return None
    keepers = []
    xcludes = []
    for name in (filterbyname or []):

        ## special case for "Europe-minus-United-Kingdom"
        patt,_,xpat = name.partition("-minus-")
        if patt != "Global":
            keepers.append(patt)
        if xpat:
            xcludes.append(xpat)
    xcludes.extend(xfilterbyname or [])

    df_keepers = []
    for name in keepers:
        if "." in name:
            ## special case for "Asia.India"
            ## as a way to avoid "India" matching "Indiana"
            cont,_,loca = name.partition(".")
            dfk = df_cases[(df_cases.location == loca)&
                           (df_cases.continent == cont)]
        else:
            dfk = df_cases[(df_cases.location == name)|
                           (df_cases.continent == name)]
        df_keepers.append(dfk)

    if df_keepers:
        df_cases = pd.concat(df_keepers)
    for name in xcludes:
        df_cases = df_cases[(df_cases.location != name)&
                            (df_cases.continent != name)]

    return df_cases

def case_counts(df_cases,ord_range,daysperweek=7):
    '''return array of total weekly (or daily or cumulative) case counts'''
    if df_cases is None:
        return None
    ord_min,ord_max = ord_range
    cum_cases = [0]*(ord_max+1-ord_min+daysperweek)
    for day,count in zip(df_cases.loc[:,'date'],
                         df_cases.loc[:,'total_cases']):
        date = emu.date_fromiso(day)
        ord_day = date.toordinal()
        if ord_min-daysperweek <= ord_day <= ord_max:
            cum_cases[ord_day-ord_min+daysperweek] += count
    total_cases = max(cum_cases) - cum_cases[daysperweek]
    v.vprint('total cases:',total_cases,max(cum_cases),
             cum_cases[ord_max-ord_min+daysperweek],
             cum_cases[daysperweek])
    v.vprint('dates:',
             datetime.date.fromordinal(ord_min-daysperweek),
             datetime.date.fromordinal(ord_min),
             datetime.date.fromordinal(ord_max))

    if daysperweek==0:
        return cum_cases,total_cases
    num_cases = [0] * (ord_max-ord_min+1)
    for n in range(ord_max-ord_min+1):
        num_cases[n] = cum_cases[n+daysperweek] - cum_cases[n]
    return num_cases,total_cases

def _main(args):
    '''read raw OWID file and write GISAID harmonized file'''
    if args.input == "OWID":
        args.input = STD_URL
    df_owid = read_dataframe(args.input)
    df_gisaid = harmonize_dataframe(df_owid)
    if args.output:
        v.print(f"Writing harmonized file to: {args.output}")
        df_gisaid.to_csv(args.output,index=False)
    else:
        v.print("Harmonized data NOT written to file")

if __name__ == "__main__":
    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
