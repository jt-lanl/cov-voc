'''Expand the table/summary from newmuts by adding columns from Bloom's database'''
import re
import argparse
import pandas as pd
import pathlib

import verbose as v
import numu
import wildtypes

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="name of input newmuts summary file")
    paa("--bloom","-b",
        help="name of Bloom's database")
    paa("--clade",
        help="name of clade used as reference in Bloom's database")
    paa("--output","-o",
        help="name of output file that has Bloom data merged in")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    df = pd.read_table(args.input)
    df.insert(0,"tmstring", None)
    for nr,row in df.iterrows():
        tmstring = re.sub(r'^[A-Z]','',row["mutation"])
        df.loc[nr,"tmstring"] = tmstring

    if args.verbose:
        df.to_csv('df.tsv',sep='\t',index=False)

    init_columns = [numu.match_column_name(df.columns,name)
                    for name in ["parent","child","lineage_seq","total_seq"]]

    if args.bloom and pathlib.Path(args.bloom).is_file():
        ## merge based on truncated m-string (site+mutation, no wild-type)
        ## will need to add tmstring columns for both dataframes

        jdf = pd.read_table(args.bloom)
        jdf.insert(0, "tmstring", None)
        for nr,row in jdf.iterrows():
            tmstring = str(row.site) + row.mutant
            jdf.loc[nr,"tmstring"] = tmstring

        if args.verbose:
            jdf.to_csv('jdf.tsv',sep='\t',index=False)
        df = df.merge(jdf,how='outer',on="tmstring")
        if args.verbose:
            df.to_csv('dfm.tsv',sep='\t',index=False)

    ## replace NA with 0 for init_columns
    for col in init_columns:
        for nr,row in df.iterrows():
            if pd.isna(df.loc[nr,col]):
                df.loc[nr,col] = 0

    if args.clade:
        wild = numu.match_column_name(df.columns,'wild')
        if wild:
            df = df.rename(columns={wild:"clade-wildtype"})
        else:
            df.insert(0,"clade-wildtype",None)
        df.insert(df.columns.get_loc("clade-wildtype"),
                  "wuhan-wildtype",None)
            
    ## maybe not complete -- what about clade wildtype for
    ## items in df not in jdf ??  will need to make a table
    ## of wildtype by site (either via covid's table, or
    ## learning from jdf file

    if "lineage_transitions" in df.columns:
        ## put lineage transitions as the last column
        col = df.pop("lineage_transitions")
        df.insert(len(df.columns),"lineage_transitions",col)

    if "clade-wildtype" not in df.columns: ## really, why?
        df.insert(1,"clade-wildtype",None)
    if "mutant" not in df.columns:
        df.insert(1,"mutant",None)
    if "mutation" not in df.columns:
        df.insert(1,"mutation",None)
    if "site" not in df.columns:
        df.insert(0,"site",None)
    for nr,row in df.iterrows():
        site = int(re.sub('\D*(\d+)\D*',r'\1',row["tmstring"]))
        df.loc[nr,"site"] = site
        ww = wildtypes.wuhan_wildtype(site)
        df.loc[nr,"wuhan-wildtype"] = ww
        mutant = re.sub('\D*\d+(\D*)',r'\1',row["tmstring"])
        df.loc[nr,"mutant"] = mutant
        if not row["tmstring"].startswith("+"):
            df.loc[nr,"mutation"] = ww + row["tmstring"]
        if args.clade:
            theclade = re.sub(r'\-.*','',args.clade) ## strip off anything after dash
            df.loc[nr,"clade-wildtype"] = \
                wildtypes.clade_wildtype(theclade,site)
        
    df = df.drop(labels=['site_x','site_y','tmstring'],
                 axis='columns',errors='ignore')

    df = df.sort_values(by="site")

    if args.output:
        df.to_csv(args.output,sep="\t",index=False)
        

if __name__ == "__main__":
    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
