'''Make a scatterplot of two columns of a pandas dataframe'''

import argparse
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import verbose as v
import numu

_colormap = mpl.colormaps.get_cmap('gist_rainbow')
_colors = _colormap(np.linspace(0,1,256))

def set_colormap(name,ntotal,noff=0):
    global _colormap,_colors
    _colormap = plt.get_cmap(name,ntotal+noff)
    _colors = _colormap(np.linspace(0,1,ntotal+noff))
    if noff > 0:
        _colors = _colors[noff:]
        _colormap = mpl.colors.ListedColormap(_colors)

    
def nth_color(n,ntotal):
    n = n % ntotal
    return _colors[n]
    

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input tsv file")
    paa("--columns","-c",nargs=2,
        help="name of the two columns")
    paa("--ncolors",type=int,default=1,
        help="how many colors")
    paa("--colorscolumn","--cc",
        help="column name of quantity used for assigning colors")
    paa("--cmap",default="gist_rainbow",
        help="color map for ncolors>1")
    paa("--xlog",action="store_true",
        help="plot x-axis with log scale")
    paa("--ylog",action="store_true",
        help="plot y-axis with log scale")
    paa("--output","-o",
        help="output plot to file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args


def segment(N,M):
    '''Given an array [0,...,N] with N+1 elements, produce a shorter array
    with m segments (so array with have m+1 elements), 
    that are roughly equally spaced, 
    begining with 0 and ending with N'''
    nstart=0
    nsegs=M
    yield nstart
    while nstart < N:
        delta_n = (N-nstart)//nsegs
        nstart = nstart + delta_n
        nsegs -= 1
        v.vprint('segment:',delta_n,nstart)
        yield nstart
        
    
        
def _main(args):
    '''main'''
    v.vprint(args)

    try:
        df = pd.read_table(args.input)
    except FileNotFoundError:
        v.print(f'File {args.input} not found!')
        return
    
    column_names = [numu.match_column_name(df.columns,name)
                    for name in args.columns]

    xname,yname = column_names

    if xname is None or yname is None:
        v.print(f'Will not output {args.output}')
        return

    xvals=[]
    yvals=[]
    zvals=[]
    for nr,row in df.iterrows():
        xvals.append( df.iloc[nr][xname] )
        yvals.append( df.iloc[nr][yname] )

    zcolors = ['Black'] * len(xvals)

    if args.ncolors and not args.colorscolumn:
        args.colorscolumn = yname

    if args.colorscolumn:
        colorscolumn = numu.match_column_name(df.columns,args.colorscolumn)
        zvals=[]
        for nr,row in df.iterrows():
            zvals.append( float(df.iloc[nr][colorscolumn]) )
        ncolors = args.ncolors
        noff = max([1,int(ncolors/10)])
        set_colormap(args.cmap,ncolors,noff=noff)
        pct = [100*i/ncolors for i in range(ncolors+1)]
        zcut = np.nanpercentile(zvals,pct,
                                method='closest_observation')

        for n in range(len(zcut)-1):
            ## hack, to avoid duplications
            if zcut[n+1] <= zcut[n]:
                zcut[n+1] = zcut[n]+1
                
        v.vprint('z cuts:',zcut)
        for nc in range(ncolors):
            for nr,zval in enumerate(zvals):
                if zcut[nc] <= zval <= zcut[nc+1]:
                    zcolors[nr] = nth_color(nc,ncolors+1)

    plt.figure()
    plt.scatter(xvals,yvals,c=zcolors)
    if args.xlog:
        plt.xscale('log')
    if args.ylog:
        plt.yscale('log')
    plt.xlabel(xname)
    plt.ylabel(yname)

    if args.colorscolumn:
        if ncolors > 10:
            ndx = list(segment(ncolors,7))
            ticks = [n/ncolors for n in ndx]
            zcut = [zcut[n] for n in ndx]
        else:
            ticks = np.linspace(0,1,len(zcut))
            
        fmt_zcut = ["%d" % z for z in zcut]
        plt.colorbar(mpl.cm.ScalarMappable(cmap=_colormap),
                     ax=plt.gca(),
                     ticks=ticks,
                     format=mpl.ticker.FixedFormatter(fmt_zcut),
                     label=colorscolumn)

        
    plt.tight_layout()
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
    

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
