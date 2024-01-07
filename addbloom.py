'''Expand the table/summary from newmuts by adding columns from Bloom's database'''
import re
import argparse
import pandas as pd

import verbose as v
import bloom

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="name of input newmuts summary file")
    paa("--bloom","-b",
        help="name of Bloom's database")
    paa("--output","-o",
        help="name of output file that has Bloom data merged in")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

    jdf = pd.read_table(args.bloom)
    bloom_summary = bloom.BloomSummary(args.bloom)
    jdf = bloom_summary.df
    wild = bloom_summary.match_column_name('wild')
    jdf.insert(0, "mstring", None)
    for nr,row in jdf.iterrows():
        jdf.loc[nr,"mstring"] = row[wild]+str(row.site)+row.mutant

    df = pd.read_table(args.input)
    BLOOM_HEADINGS = []
    for h in ("escape","entry","binding"):
        try:
            hx = bloom_summary.match_column_name(h)
            BLOOM_HEADINGS.append(hx)
        except RuntimeError:
            v.print(f'No "{h}" column in {args.bloom}')


    if "site" not in df.columns:
        df.insert(0,"site",None)
    if "clade-wildtype" not in df.columns:
        df.insert(1,"clade-wildtype",None)
    for h in BLOOM_HEADINGS:
        df.insert(len(df.columns),h,None)
        
    for nr,row in df.iterrows():
        ## strip first letter since Bloom uses a different
        ## definition of wildtype (he calls wildtype-BA2)
        jmstring = [s[1:] for s in jdf["mstring"]]
        v.vvprint(f'{row["mutation"]=}')
        ddmstring = row["mutation"].strip()
        v.vvprint(f'{ddmstring=}')
        dmstring = ddmstring[1:]
        v.vvprint(f'{dmstring=}')
        v.vvprint(f'jmstring: {len(jmstring)} {jmstring[:5]}')
        v.vvprint(f'dmstring: {len(dmstring)} {dmstring[:5]}')
        v.vvvprint(f'True or False:',[jms == dmstring for jms in jmstring])
        jwild = [s[0] for s in jdf["mstring"]]
        the_wild=None
        for wild,jms,mstring in zip(jwild,jmstring,jdf["mstring"]):
            if jms == dmstring:
                the_wild = wild
            if jms == dmstring and wild != ddmstring[0]:
                v.vprint('wild:',wild,jms,ddmstring,dmstring,mstring)
                if mstring != ddmstring:
                    v.vvprint('WILD:',wild,jms,ddmstring,dmstring,mstring)
            
        jrow = jdf[[jms == dmstring for jms in jmstring]]
        
        #row["mutation"].strip()[1:]]
        v.vvprint_only(5,'jrow',f'{nr=},{row["mutation"]=}')
        
        site = int(re.sub('\D*(\d+)\D*',r'\1',df.loc[nr,"mutation"]))
        df.loc[nr,"site"] = site
        if jrow.empty:
            continue
        df.loc[nr,"clade-wildtype"] = the_wild
        for h in BLOOM_HEADINGS:
            df.loc[nr,h] = jrow.iloc[0][h]

    if args.output:
        df.to_csv(args.output,sep="\t",index=False)
        

if __name__ == "__main__":
    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
