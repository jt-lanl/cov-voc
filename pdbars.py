'''Make a bar plot of some column in a tsv file as fcn of site'''

import argparse
import pandas as pd
import matplotlib.pyplot as plt

import verbose as v
import numu

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input tsv file")
    paa("--column","-c",
        help="name of the column")
    paa("--sortcolumn",
        help="use this column to define top values")
    paa("--bysite",action="store_true",default=True,
        help="Plot by site number")
    paa("--bymstring",action="store_false",dest='bysite',
        help="Plot by individual mstring")
    paa("--top","-t",type=int,default=50,
        help="Only plot the highest-valued points")
    paa("--rbd",action="store_true",
        help="Plot only RBD region (328-528) to plotfile")
    paa("--title",
        help="Title to put on plot")
    paa("--ylog",action="store_true",
        help="plot y-axis with log scale")
    paa("--output","-o",
        help="write plot to file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args


def _main(args):
    '''main'''
    v.vprint(args)

    df = pd.read_table(args.input)
    
    column_name = numu.match_column_name(df.columns,args.column)
    if column_name is None:
        v.print(f'Output file {args.output} wil not be created')
        return

    plt.figure(figsize=(18,5))
    bysite = None
    if args.bysite:
        bysite = numu.BySiteSummary()
        for ssm,val in zip(df["mutation"],df[column_name]):
            bysite.add_mutation(ssm,val)
        sitelist = sorted(bysite.val.keys())
        if args.sortcolumn and not args.rbd:
            sort_column = numu.match_column_name(df.columns,args.sortcolumn)
            xbysite = numu.BySiteSummary()
            for ssm,val in zip(df["mutation"],df[sort_column]):
                xbysite.add_mutation(ssm,val)
        else:
            xbysite = bysite
            
        if args.top:
            sitelist = xbysite.topsites(n=args.top,thresh=1)
        if args.rbd:
            sitelist = list(range(328,529))
        v.vprint(f'{len(sitelist)=}')
        if not sitelist:
            v.print('Empty sitelist, no sites satisfy threshold')
            v.print(f'File {args.output} will not be created')
            return
        bysite.plotsites(sitelist)
        plt.xlabel('site')
    else:
        sorted_df = df.sort_values(by=column_name,
                                   ascending=False,
                                   na_position='last')
        sorted_df = sorted_df[:args.top]
        sorted_df = sorted_df.sort_values(by="site",ascending=True)
        
        xvals = []
        yvals = []
        for _,row in sorted_df.iterrows():
            xvals.append(row["mutation"])
            yvals.append(row[column_name])

        plt.bar(range(len(xvals)),yvals)
        plt.xticks(range(len(xvals)),xvals,
                   rotation='vertical')
        plt.xlim([-1,len(xvals)])
        plt.xlabel('mutation')
        
        plt.gca().spines[['right', 'top']].set_visible(False)



    plt.ylabel(column_name)
    if args.ylog:
        plt.yscale('log')
    if args.title:
        plt.title(args.title)
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
        
if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
