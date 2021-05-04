import numpy as np
import itertools as it
import re
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib import cm
import mandalacolor as mc

def circleplot(values,names,nodevalues=None,**kwargs):
    ne = len(names)
    assert values.shape == (ne,ne)
    ## make a list of non-diagonal values, used for scaling
    vnd = [values[i,j] for i,j in it.combinations(range(ne),2)]
    vlo = np.median(vnd)  ## deliberately cutting off half of the low ones
    vlo = np.percentile(vnd, max([50, 100*(1 - 10/ne)]))
    vhi = np.max(vnd)
    
    theta = 2*np.pi / ne
    x = 0.8*np.sin( theta*np.arange(ne) )
    y = 0.8*np.cos( theta*np.arange(ne) )    
    
    ## set up to plot edges
    ## (but plot them in order of least value first)
    xyc = []
    for i,j in it.combinations(range(ne),2):
        if values[i,j] > vlo:
            ## nb: x[[i,j]] = [x[i],x[j]]
            xyc.append( (x[[i,j]],y[[i,j]],values[i,j]) )
    xyc = sorted(xyc, key=lambda w: w[2])

    plt.figure(figsize=(10,10))

    lines=[]
    colors=[]
    for xx,yy,c in xyc:
        #plt.plot(xx,yy,c=cm.Blues(c),**kwargs,lw=0.5)
        lines.append( ((xx[0],yy[0]),(xx[1],yy[1])) )
        colors.append(c)

    segments = LineCollection(lines,cmap=cm.Blues,lw=3)
    segments.set_array(np.array(colors))
    plt.gca().add_collection(segments)

    if nodevalues is None:  
        #default node values are diagonal elements of input table
        nodevalues = [values[i,i] for i in range(ne)]

    ## But nodevalues are on a different color scale, generally, 
    ## than the lines, so go ahead and rescale them to [0,1]...
    ## make that [0.2,0.8] so they are all gray
    ndx = np.argsort(nodevalues)
    nodevalues = np.linspace(0.2,0.8,len(nodevalues))[ndx]

    ## pull the integer part out of the name; eg "L452R" becomes 452.
    sites = [int(re.sub("\D*(\d+)\D*",r"\1",name)) for name in names]
    
    ## plot nodes
    plt.plot(x,y,'o',ms=9,c=cm.Blues(1.0))
    for i in range(ne):
        plt.plot([x[i]],[y[i]],'o',ms=7,c = cm.Greys(nodevalues[i]))
        plt.text(1.2*x[i],1.15*y[i]-0.01,names[i],
                 va="center",ha="center",
                 color=mc.getcolor(sites[i]),
                 fontsize='small',
                 fontweight=mc.getweight(sites[i]))

    plt.plot([-1,-1,1,1,-1],[-1,1,1,-1,-1],'-',lw=0)
    plt.xticks([])
    plt.yticks([])
    plt.axis('scaled')
    for spine in plt.gca().spines.values():
        spine.set_visible(False)

    plt.colorbar(segments,shrink=0.8,orientation='horizontal',pad=-0.05)
        
def heatmap(values,names,nodevalues=None,**kwargs):
    ne = len(names)
    assert values.shape == (ne,ne)

    vnd = [values[i,j] for i,j in it.combinations(range(ne),2)]
    vrange = [np.min(vnd), np.median(vnd), np.max(vnd)]


    plt.figure(figsize=(10,10))
    plt.imshow(values,vmin=vrange[0],vmax=vrange[2],**kwargs)
    plt.xticks(np.arange(ne),names,rotation=90)
    plt.yticks(np.arange(ne),names)
    
    plt.colorbar(shrink=0.8)


