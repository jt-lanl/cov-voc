'''
Requires python 3.6+

Typical use:
  python twobars.py variant.S.614.country.txt --plot
Full list of command-line options:
  python twobars.py -h
'''

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
import argparse
import datetime
from collections import Counter
import scipy.stats as stats

def getargs():
    ap = argparse.ArgumentParser()
    paa = ap.add_argument
    paa("infile",
        help="Input table")
    paa("--delay","-d",type=int,default=14,
        help="Number of days to wait to observe possible growth")
    paa("--range","-r",nargs=2,help="range of dates (two dates, ISO format)")
    paa("--mutants","-m",nargs=2,default=["D614","G614"],
        help="list of (exactly 2) mutant variants to be compared")
    paa("--mincount",type=int,default=15,
        help="minimum count before and after onset")
    paa("--verbose","-v",action="count",help="verbose")
    paa("--plot",action="store_true",help="make plot as well as table")
    paa("--nolegend",action="store_true",help="do not include legend in plots")
    paa("--output","-o",help="write plot to file")
    paa("--title","-t",help="Title")
    paa("--pvalue",type=float,default=0.05,
        help="Threshold for significance (use pvalue=1.1 for no threshold)")
    paa("--level",type=int,default=0,
        help="level=1, X; level=2, X_Y; level=3, X_Y_Z")
    args = ap.parse_args()
    return args

def date_fromiso(s):
    yyyy,mm,dd = s.split("-")
    dt = datetime.date(int(yyyy),int(mm),int(dd))
    return dt

def date_fromshef(s):
    m,d,yy = s.split("/")
    dt = datetime.date(2000+int(yy),int(m),int(d))
    return dt

def date_friendly(dt):
    #return dt.strftime("%b %-d")
    mmm = dt.strftime("%b")
    dd = dt.strftime("%-d")
    return "%3s %2d" % (mmm, int(dd))


def ord2date(ord):
    return date_friendly(datetime.date.fromordinal(ord))

def uniq_patterns(plist):
    '''if two or more items in the list match, then only the last one is kept;
    note that "match" is case insensitive'''
    puniq = []
    for i,p in enumerate(plist):
        for q in plist[i+1:]:
            if q.lower() == p.lower():
                vprint("Match:",p,q)
                break
        else:
            puniq.append(p)
    return puniq

def pattshort(p):
    Abbrevs = [
        ## note: applied sequentially so order matters
        ("Democratic-Republic-of-the-Congo","DRC"),
        ("New-South-Wales","NSW"),
        ("Washington","WA"),
        ("New-York","NY"),
        ("California","CA"),
        ("San-Francisco","SF"),
        ("United-Kingdom","UK"),
        ("Comunitat-Valenciana","Com-Valenciana"),
        ("Wisconsin","WI"),
        ("Hauts-de-France","HdF"),
        ("England_LIVE","England_Liverpool"),
        ("England_NOTT","England_Nottingham"),
        ("England_SHEF","England_Sheffield"),
        ("England_NORW","England_Norwich"),
        ("England_NORT","England_Northampton"),
        ("England_EXET","England_Exeter"),
        ("England_BRIS","England_Bristol"),
        ("England_CAMB","England_Cambridge"),
        ("England_LOND","England_London"),
        ("England_OXON","England_Oxford"),
        ("England_PORT","England_Portsmouth"),
    ]
    for longpatt,shortpatt in Abbrevs:
        p = re.sub(longpatt,shortpatt,p)
    
    return p

def XY_pattern(patt):
    ''' input: X_Y_Z_W; return: X_Y '''
    tk = patt.split("_")
    if len(tk) > 2:
        return "_".join(tk[:2])
    else:
        return patt

def X_pattern(patt):
    ''' input: X_Y_Z_W; return: X '''
    tk = patt.split("_")
    if len(tk) > 1:
        return tk[0]
    else:
        return patt

def filter_patterns(vList):
    ## This finds extra patterns, like USA_CA, which don't appear in v.info strings
    ## It takes longer, and introduces new patterns, eg: Australia_NSW even
    ## though virtually all of Australia_NSW is Australia_NSW_Sydney
    ## which leads to a kind of duplication, but some of these duplicates
    ## are then removed, using the logic below

    ## Begin with all patterns that appear in info strings
    pattset = set(v.info for v in vList)
    vprint("Number of patterns:",len(pattset))



    ## Mostly checks for multiple naems for the same place;
    ## Currently only based on lower case, so example duplicate
    ## is "de-la-Loire" vs "De-La-Loire"
    puniq = uniq_patterns(list(pattset))
    vprint("Uniq patterns:",len(puniq))

    ## subpatt: X_Y_Z -> X_Y
    XYpatt = set( XY_pattern(p) for p in pattset )
    newpatt = XYpatt - pattset
    XYZpatt = pattset - XYpatt

    vprint("  X_Y_Z patterns:",len(XYZpatt))
    vprint("All X_Y patterns:",len(XYpatt))
    vprint("New X_Y patterns:",len(newpatt))
    vvprint("New X_Y:",sorted(newpatt))

    pattcount = Counter()
    for p in pattset | newpatt:
        for v in vList:
            if re.match(p,v.info):
                pattcount[p] += 1

    ## Find duplicate (and near duplicate) patterns; eg
    ## if all of X_Y is counted as X_Y_Z,
    ##    then prefer X_Y_Z, consider X_Y a dup
    ## if X_Y is mostly X_Y_Z, but there are other X_Y_*'s,
    ##    then prefer X_Y, consider X_Y_Z a dup
    ## we'll keep those other X_Y_{non-Z} patterns, though
    ## most of them will probably not have enough counts to
    ## make it to the final list
    dup_patt = set()
    for XYZ in XYZpatt:
        XY = XY_pattern(XYZ)
        if pattcount[XYZ] == pattcount[XY]:
            vprint("Identical:",XY,XYZ,pattcount[XY])
            dup_patt.add(XY)
        elif pattcount[XYZ] > 0.8*pattcount[XY]: ## Hardcoded 0.8
            dup_patt.add(XYZ)

    if 0: ## Hmm, I guess we are not going to do this
        ## Find cases of X_Y where Y is exclusive; ie no other X_*'s
        ## eg Norway_Oslo; consider X_Y a duplicate since it'll be in
        ## the X item
        Xpatt = set(X_pattern(p) for p in XYpatt)
        for X in Xpatt:
            Xchildren = [p for p in XYpatt if X_pattern(p)==X]
            if len(Xchildren)==1:
                print("Single region in country:",X,Xchildren)
                dup_patt.add(Xchildren[0])

    ## Ok, eliminate all those dup's
    vprint("Dup:",sorted(dup_patt))
    keeppatt = pattset | newpatt
    keeppatt = keeppatt - dup_patt

    return keeppatt,pattcount

def truncate_pattern(patt,level):
    '''return the subpattern that is at the specified level; eg
       level=1; X_Y_Z -> X
       level=2; X_Y_Z -> X_Y, X_Y -> X_Y, X -> None
    '''
    tk = patt.split("_")
    if len(tk) >= level:
        return "_".join(tk[:level])
    else:
        return None

def filter_patterns_bylevel(vList,level):
    if level==0:
        keep,count = filter_patterns(vList)
        return keep
    if 1<=level<=4:
        lpat = [truncate_pattern(v.info,level) for v in vList]
        lpat = set(p for p in lpat if p is not None)
        #print("levelpatt:",lpat)
        return lpat
    raise RuntimeError("level=%d, should be 0<=level<=4" % level)

class virion:
    def __init__(self,dg,info,n,date):
        self.dg = dg      ## D614 or G614
        self.info = info  ## Locatin information
        self.n = n        ## EPI_ISL number
        self.date = date

    def __repr__(self):
        return "%4s %25s %s %s" % (self.dg, self.info, self.n, self.date.isoformat() )

def read_virions(filename):
    vList = []
    with open(filename,"r") as fin:
        for line in fin:
            tokens=line.strip().split()
            if len(tokens) < 4:
                vprint("Incomplete line:",line.strip())
                continue
            dg,info,n,d = tokens
            try:
                dt = date_fromiso(d)
            except:
                try:
                    dt = date_fromshef(d)
                except:
                    vvprint("Incomplete date:",d,n)
                    continue

            vList.append( virion(dg,info,n,dt ) )
    return vList
            
def main(args):

    D,G = args.mutants if args.mutants else ("D614","G614")
    #DG = [D,G]

    vList = read_virions(args.infile)
    vprint("Number of entries:",len(vList))

    if 1:
        keeppatt = filter_patterns_bylevel(vList,args.level)
    else:
        keeppatt = set(pattcount.keys())

    ## Print Header
    print(" "*16,"          Onset    |--Before-|   Delay   |--After--|   Last    Delta   Fisher")
    #hy 2020-08-17 print(" "*16,"Location  Date     G/(G+D) G+D   Date    G/(G+D) G+D   Sample  G/(G+D)  p-val")
    d = D[0:1]
    g = G[0:1]
    print(" "*16,"Location  Date     "+g+"/("+g+"+"+d+") "+g+"+"+d+"   Date    "+g+"/("+g+"+"+d+") "+g+"+"+d+"   Sample  "+g+"/("+g+"+"+d+")  p-val")
    
    ## Initialize lists (we'll be plotting these)
    good_patt = []
    f_befor = []
    f_after = []
    pvalues = []

    for patt in sorted( keeppatt ):

        pattcount = sum(re.match(patt,v.info) is not None for v in vList)

        if pattcount < 2*args.mincount:
            vvprint(patt,"Too few counts:",pattcount)
            continue

        DG_datelist = {D: [], G: []}
        for v in vList:
            if not re.match(patt,v.info, re.IGNORECASE | re.ASCII):
                continue
            if v.dg in DG_datelist:
                DG_datelist[v.dg].append(v.date)

        count = sum( len(datelist) for datelist in DG_datelist.values() )
        vvprint(f"{patt}: Counted {count}/{pattcount} good entries;",end=" ")
        for m in D,G:
            vvprint(m+":",len(DG_datelist[m]),end=" ")
        vvprint("")

        if count==0:
            continue

        alldates = set( dt for m in (D,G) for dt in DG_datelist[m] )
        vvprint("Range dates:",min(alldates),max(alldates))
        ordmin = min(alldates).toordinal()
        ordmax = max(alldates).toordinal()

        if args.range:
            ordplotmin = date_fromiso(args.range[0]).toordinal()
            ordplotmax = date_fromiso(args.range[1]).toordinal()
            ordmin = min([ordmin,ordplotmin])
        else:
            ordplotmin = ordmin
            ordplotmax = ordmax

        ordmax = max([ordmax,ordplotmax])

        Ndays = ordmax+1-ordmin

        ## Get cumulative counts
        DG_cum = {D: [], G: []}  
        for ord in range(ordmin,ordmax+1):
            day = datetime.date.fromordinal(ord)
            for m in D,G:
                DG_cum[m].append( len([dt for dt in DG_datelist[m] if dt <= day]) )

        ## Only keep data that is within the specified date range
        if args.range:
            for m in D,G:
                ztmp = []
                for i,cnt in enumerate(DG_cum[m]):
                    if ordplotmin <= ordmin+i <= ordplotmax:
                        ztmp.append(cnt)
                DG_cum[m] = ztmp
            ordmin = ordplotmin
            Ndays = ordplotmax+1-ordmin

        ## Now that DG_cum is computed, we can define GD_ratio
        def GD_ratio(day,start_day=None):
            ''' return ratio of G/(G+D), as well as G,D counts;
            cumulative total between start_day and day
            '''
            count_G = DG_cum[G][day]
            count_D = DG_cum[D][day]
            if start_day:
                count_G -= DG_cum[G][start_day]
                count_D -= DG_cum[D][start_day]
            f = count_G / (count_G + count_D) if count_G>0 else 0
            return f,count_G,count_D

        ## From now on, refer to shortened names
        patt = pattshort(patt)

        ## Find onset: first day with G,D both >= ONSET_DG, and G+D > mincount
        ONSET_DG = 3
        first_day = None
        for day in range(Ndays):
            f,g,d = GD_ratio(day)
            if d+g >= args.mincount and (g>=ONSET_DG and d>=ONSET_DG):
                first_day = day
                break
        if first_day is None:
            vprint(patt,"Onset criteria not met")
            continue
        
        last_day = Ndays-1
        later_day = first_day + args.delay
        if not first_day < later_day < last_day:
            vprint(patt,"not enough time")
            continue

        ## Counts that make up the second bar
        ff,gg,dd = GD_ratio(last_day,later_day)
        if gg+dd < args.mincount:
            vprint(patt,"not enough counts:",gg+dd)
            continue

        ## Get p-value
        table = [[d, g],[dd, gg]]
        _,pval = stats.fisher_exact(table)
        if pval >= args.pvalue:
            vprint(patt,"not significant: p=",pval)
            continue

        print( "%25s %7s   %.3f %5d  %7s  %.3f %5d  %7s  %5.2f %.6f" %
               (patt,
                ord2date(first_day + ordmin), f, d+g,
                ord2date(later_day + ordmin), ff, dd+gg,
                ord2date(last_day + ordmin), ff-f, pval))

        ## Update all our lists
        good_patt.append(patt)
        f_befor.append(f)
        f_after.append(ff)
        pvalues.append(pval)

    if len(good_patt)==0:
        print("No regions satisfied constraints")
        return
        
    if not args.plot:
        return

    ## Convert lists into numpy arrays
    x = np.arange(len(good_patt))
    f_befor = np.array(f_befor)
    f_after = np.array(f_after)
    pvalues = np.array(pvalues)

    ## Re-sort based on fraction in first bar
    ndx = np.argsort(f_befor)
    good_patt = [good_patt[n] for n in ndx]
    f_befor = f_befor[ndx]
    f_after = f_after[ndx]
    pvalues = pvalues[ndx]

    max_patt_length = max([len(p) for p in good_patt])
    vprint("max pattern length:",max_patt_length)

    figheight = 2.6 + 0.07*max_patt_length
    figwidth = 2.5+0.35*len(ndx)
    if args.nolegend:
        figwidth -= 1.0

    fig,ax=plt.subplots(figsize=(figwidth,figheight))
    
    blue = '#1f77b4'
    orange = '#ff7f0e'
    if 1:
        ## double bar plots with orange and blue
        barwid=0.33
        plt.bar(x-0.2,1-f_befor,width=barwid,color=orange,label=D)
        plt.bar(x-0.2,f_befor,bottom=1-f_befor,width=barwid,color=blue,label=G)

        plt.bar(x+0.2,1-f_after,width=barwid,color=orange)
        plt.bar(x+0.2,f_after,bottom=1-f_after,width=barwid,color=blue)

    else:
        ## single bar with only blue
        plt.bar(x-0.2,f_befor,width=0.38,color=blue)
        plt.bar(x+0.2,f_after,width=0.38,color=blue)
        plt.xticks(x,good_patt,rotation=45,ha='right',position=(0,0.02))

        ndxrev = f_after < f_befor
        plt.bar(x[ndxrev]-0.2,f_befor[ndxrev],width=0.38,color='violet')
        plt.bar(x[ndxrev]+0.2,f_after[ndxrev],width=0.38,color='violet')
        
    ## write p-value diagonally across top of bar
    for i in range(len(x)):
        if pvalues[i] < 0.00001:
            pvalstr = " p<0.00001"
        else:
            pvalstr = " p=%.5f"%pvalues[i]
        plt.text(x[i]-0.4,1,pvalstr,
                 ha='left', va='bottom',
                 transform=ax.transData,rotation=45)

    ## write names of regions diagonally across bottom of bars
    plt.xticks(x,good_patt,rotation=45,ha='right',position=(0,0.02))

    ## adjust xlimit to make room for legend, if legend desired
    if args.nolegend:
        plt.xlim([-0.7,len(x)+0.75])
    else:
        plt.xlim([-3,len(x)+0.75])

    plt.ylim([0,1.5])
    plt.yticks(np.arange(0,1.01,0.2))
    plt.ylabel("fraction %s vs %s" % (D,G))
    if args.title:
        title = args.title
    else:
        title = "onset-threshold = %d/%d, delay = %d days" %\
            (ONSET_DG, args.mincount, args.delay)
        if args.pvalue <= 1.0:
            title = title + ", p<%.2g"%(args.pvalue,)
    #plt.title(title)
    if not args.nolegend:
        plt.legend(loc="upper left")
    plt.tight_layout()    
    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()
    

if __name__ == "__main__":
    args = getargs()

    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,**kw)
        
    def vvprint(*p,**kw):
        if args.verbose and args.verbose > 1:
            print(*p,file=sys.stderr,*kw)
            
    main(args)

