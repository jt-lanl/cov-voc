import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import datetime

def date_friendly(dt):
    #return dt.strftime("%b %-d")
    mmm = dt.strftime("%b")
    dd = dt.strftime("%-d")
    return "%3s %2d" % (mmm, int(dd))

def half_labels(labels,n=2):
    '''replace list of labels ["one","two","three","four","five"]
    with ["one", "", "three", "", "five"]
    '''
    return [label if i%n==0 else "" 
            for i,label in enumerate(labels)]

def embersplot(counts,
               fullnames,
               mcolors,
               ordmin,
               ordplotrange=None,
               title=None,legendtitle=None,legend=0,onsets=False,
               fraction=False,lineplot=False,show=False):
    '''
    embersplot is an embers-style (stacked bar or line) plot of counts vs date for multiple variants
    counts: dictionary of arrays; each array is counts vs date; keys of dictionary correspond to different variants
    fullnames: dictionary of names associated with dictionary keys (=None if keys are names)
    mcolors: dictionary of colors
    ordmin: first day in array, as ordinal number
    ordplotrange: tuple (ordplotmin, ordplotmax) indicates range of days to plot
    master: if specified, plot a dummy with this as name in legend (can be used as a lable)
    '''

    patterns = list(counts)
    fullnames = fullnames if fullnames else {m:m for m in patterns}
    maxnamelen = max(len(n) for n in fullnames.values())
    namefmt = "%%-%ds" % maxnamelen
    mnames = {m : namefmt % fullnames[m] for m in patterns}

    Ndays = len(counts[patterns[0]]) 
    for m in patterns: ## should make sure len is same for all patterns
        assert Ndays == len(counts[m])

    if legend == 0:
        plt.figure(figsize=(6,3))
    elif legend == 1:
        F = 4.7 if lineplot else 5.0
        plt.figure(figsize=(12,1+len(patterns)/F))
    elif legend == 2:
        F = 4.7 if lineplot else 4.5
        plt.figure(figsize=(12,1+len(patterns)/F))
        #plt.figure(figsize=(12,0.5+0.5*len(patterns)/F))
    else:
        raise RuntimeError(f"legend should be 0, 1, or 2: not {legend}")

    if title:
        #title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')

    counts_bottom = dict()
    counts_total = [0] * len(counts[patterns[0]])
    for m in counts:
        counts_bottom[m] = counts_total
        counts_total = [x+y for x,y in zip(counts_total,counts[m])]

    ## at this point counts_total is array (vs date) of total counts over all patterns

    if legendtitle and legend>1 and not lineplot:
        plt.bar(range(Ndays),[0]*Ndays,width=1,
                label=legendtitle,color="white")

    name_color_sofar = set()
    for m in  patterns[::-1]:

        name = mnames[m] ## mnames
        #if legend>1:   ### take this outside
        #    name += " " + mrelpatt[m] 
        name = " " + name ## hack! leading underscore doesn't make it to legend??

        if lineplot:
            kwargs = dict(color=mcolors[m],label=name)
        else:
            ## For bar plot, repeated name,color tuples only appear once in legend
            name_color = (name,mcolors[m])
            kwargs = dict(color=mcolors[m])
            if name_color not in name_color_sofar:
                name_color_sofar.add(name_color)
                kwargs['label']=name            

        if fraction:
            fm = [a/(b+0.001) for a,b in zip(counts[m],counts_total)]
            bm = [a/(b+0.001) for a,b in zip(counts_bottom[m],counts_total)]
        else:
            fm = counts[m]
            bm = counts_bottom[m]

        if lineplot:
            dy = np.array(range(Ndays))
            fm = np.array(fm)
            dy = dy[fm>0]
            fm = fm[fm>0]
            plt.semilogy(dy,fm,lw=2, **kwargs)
        else:
            plt.bar(range(Ndays),fm,bottom=bm,width=1,**kwargs)

    if fraction and not lineplot:
        plt.ylim([0,1.05])

    #if fraction:
    #    plt.yticks([0.0001,0.001,0.01,0.1,1])
    #    plt.yticks([0.0001,0.01,1])

    if legend:
        plt.legend(bbox_to_anchor=(1.02, 1),
                   #handlelength=3,
                   #markerfirst=False,
                   frameon=False,
                   handletextpad=0,
                   labelspacing=0.45,
                   loc='upper left', borderaxespad=0.,
                   prop={'family' : 'monospace'})

    if not ordplotrange:
        ordplotmin = ordmin+1
        ordplotmax = ordmin+Ndays
    else:
        ordplotmin, ordplotmax = ordplotrange
        
    plt.xlim(ordplotmin-ordmin-1,ordplotmax-ordmin+1)
    xticks = list(range(ordplotmin-ordmin,ordplotmax-ordmin+1,7)) ## was n+6
    xlabels = [datetime.date.fromordinal(int(ord+ordmin)) for ord in xticks]
    xlabels = [date_friendly(dt) for dt in xlabels]
    nhalf = 1 + len(xlabels)//16
    if nhalf>1:
        xlabels = half_labels(xlabels,nhalf)
    plt.xticks(xticks,xlabels,fontsize='medium', ## was small
               rotation=45,ha='right',position=(0,0.01))

    if onsets:
        ylo,yhi = plt.gca().get_ylim()
        for m in patterns:
            if m not in onsets:
                continue
            x = onsets[m].toordinal() - ordmin
            kwargs=dict(lw=1,color=mcolors[m])
            if not fraction:
                kwargs['zorder']=0
            plt.plot([x,x],[ylo,yhi],**kwargs)

    if lineplot:
        if legend==0:
            plt.subplots_adjust(bottom=0.15,right=0.95) ## hardcoded hack!
        else:
            #plt.subplots_adjust(bottom=0.15,right=0.8) ## still hardcoded!
            #plt.subplots_adjust(bottom=0.15,right=0.7) ## still hardcoded!
            pass

        def ytickformat(x,pos):
            if x>1:
                return "" if fraction else "%g" % x
            else:
                return "%.1g" % x

        plt.gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(ytickformat))

    plt.yticks(fontsize=12)

    plt.ylabel("Fraction" if fraction else "Counts")
    plt.tight_layout()
    if show:
        plt.show()
