'''Expand the table/summary from newmuts by adding columns from Bloom's database'''
import re
import argparse
import pathlib
import pandas as pd

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
    paa("--backbone",
        help="name of lineage used for comparison")
    paa("--output","-o",
        help="name of output file that has Bloom data merged in")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

COLUMN_NAMES = "site mutation wuhan-wildtype clade-wildtype mutant parent child transitory lineage_seq total_seq denominator lineage_trans"


def match_column_names(df: pd.DataFrame, names: list[str]) -> list[str]:
    return [fullname for name in names if (fullname := numu.match_column_name(df.columns, name))]


def re_order_columns(df: pd.DataFrame, col_names: str | list[str] = None) -> pd.DataFrame:
    """re-order the columns on the dataframe so they are all in a standard order"""
    col_names = col_names or COLUMN_NAMES
    if isinstance(col_names,str):
        col_names = col_names.split()

    ## column names at the front of the list
    column_order = match_column_names(df, col_names)

    ## add any further column names not yet included
    for column_name in df.columns:
        if column_name not in column_order:
            column_order.append(column_name)

    return df.reindex(columns=column_order)


def _main(args):
    '''main'''
    v.vprint(args)

    df = pd.read_table(args.input)

    ## Create site_mut column used for merging
    df.insert(0,"merge_on_site_mut", None)
    for nr,row in df.iterrows():
        df.loc[nr,"merge_on_site_mut"] = re.sub(r'^[A-Z]','',row["mutation"])

    if args.bloom and pathlib.Path(args.bloom).is_file():
        ## merge based on truncated m-string (site+mutation, no wild-type)
        ## will need to add merge_on_site_mut columns for both dataframes

        jdf = pd.read_table(args.bloom)

        ## Create site_mut column used for merging
        jdf.insert(0, "merge_on_site_mut", None)
        for nr,row in jdf.iterrows():
            jdf.loc[nr,"merge_on_site_mut"] = str(row.site) + row.mutant

        if args.verbose:
            jdf.to_csv('jdf.tsv',sep='\t',index=False)
        df = df.merge(jdf,how='outer',on="merge_on_site_mut")
        if args.verbose:
            df.to_csv('dfm.tsv',sep='\t',index=False)

    ## replace NA with 0 for init_columns
    if args.verbose:
        df.to_csv('df.tsv',sep='\t',index=False)

    zero_empty_columns = match_column_names(df, ["parent","child","lineage_seq","transitory","total_seq"])
    for col in zero_empty_columns:
        for nr,row in df.iterrows():
            if pd.isna(df.loc[nr,col]):
                df.loc[nr,col] = 0

    ## Backbone: 1 or 0 if ssm is in the backbone (lineage) 
    ## Disabled for now
    if False and wildtypes.clade_defined(args.backbone):
        bbname = f'backbone_{args.backbone}'
        df.insert(1,bbname,None)
        backbone = [str(ssm)
                    for ssm in wildtypes.clade_ssmlist(args.backbone)]
        for nr,row in df.iterrows():
            df.loc[nr,bbname] = int(bool(row["mutation"] in backbone))

    theclade = re.sub(r'\-.*','',args.clade) if args.clade else "ALL" ## strip off anything after dash
    if args.clade:
        wild = numu.match_column_name(df.columns,'wild')
        if wild:
            df = df.rename(columns={wild:"clade-wildtype"})
        else:
            df.insert(0,"clade-wildtype",None)
        df.insert(df.columns.get_loc("clade-wildtype"),
                  "wuhan-wildtype",None)

    if "lineage_transitions" in df.columns:
        ## put lineage transitions as the last column
        col = df.pop("lineage_transitions")
        df.insert(len(df.columns),"lineage_transitions",col)

    if "clade-wildtype" not in df.columns:  ## really, why?
        df.insert(1, "clade-wildtype", None)
    if "mutant" not in df.columns:
        df.insert(1, "mutant", None)
    if "mutation" not in df.columns:
        df.insert(1, "mutation", None)
    if "site" not in df.columns:
        df.insert(0, "site", None)
    for nr, row in df.iterrows():
        site = int(re.sub("\D*(\d+)\D*", r"\1", row["merge_on_site_mut"]))
        ## if "+" in row["mutation"] then site is really site+ and ww is empty (sort of)
        df.loc[nr, "site"] = site
        ww = wildtypes.wuhan_wildtype(site)
        df.loc[nr, "wuhan-wildtype"] = ww
        df.loc[nr, "clade-wildtype"] = ww
        mutant = re.sub("\D*\d+(\D*)", r"\1", row["merge_on_site_mut"])
        df.loc[nr, "mutant"] = mutant
        if not row["merge_on_site_mut"].startswith("+"):
            df.loc[nr, "mutation"] = ww + row["merge_on_site_mut"]
        if args.clade:
            df.loc[nr, "clade-wildtype"] = wildtypes.clade_wildtype(theclade, site)

    df = df.drop(labels=['site_x','site_y','merge_on_site_mut'],
                 axis='columns',errors='ignore')

    df = re_order_columns(df,COLUMN_NAMES)
    df.insert(0,'clade',theclade)
    df = df.sort_values(by="site")

    if args.output:
        df.to_csv(args.output,sep="\t",index=False)


if __name__ == "__main__":
    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
