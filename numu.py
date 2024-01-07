'''Utilities for tracking new mutations'''

import re
from collections import defaultdict,Counter
import numpy as np
import matplotlib.pyplot as plt

import verbose as v

## Default color cycle: blue, orange, etc
_prop_cycle = plt.rcParams['axes.prop_cycle']
_colors = _prop_cycle.by_key()['color']
def nth_color(n):
    return _colors[n % len(_colors)]

def ssmlabel(ssmlist,maxlen=20):
    '''
    from a list of the form [A5B, A5C, A5X, +5DEF, ...]
    produce a label fo the form "A5: B,C,X,+DEF,..."
    '''
    xssmlist = [re.sub(r'[A-Z]*\d+','',ssm) for ssm in ssmlist]
    noins = [ssm for ssm in ssmlist if ssm[0]!='+']
    if noins:
        site = re.sub('[A-Z]*(\d+).*',r'\1',noins[0])
        initletter = noins[0][0] + site + ": "
    else:
        initletter = ''
    label = initletter+",".join(xssmlist)
    if len(label) > maxlen:
        label = label[:maxlen]
        label = re.sub(r'[^,*]$','',label) + "..."
    return label

def match_column_name(true_names,input_name):
        '''find the true column name that input name aspires to'''
        the_matches = [] ## list of true names that match
        if input_name in true_names:
            return input_name
        for true_name in true_names:
            if input_name.strip().lower() in true_name.strip().lower():
                the_matches.append(true_name)
        if len(the_matches) == 1:
            return the_matches[0]
        if len(the_matches) == 0:
            v.print(f'No match found for name={input_name}')
            return None
        if len(the_matches) > 1:
            v.print(f'Too many matches for name={input_name}:\n'
                    f'{the_matches}')
            return None



class BySiteSummary:
    '''organizes mutations by site (as dict's using site as key):
    ssmval_bysite is a dict of dicts, with ssmval_bysite[9]['V9F']
    being some value (eg, 'escape') that characterizes the mutation
    ssm_bysite is equivalent to list(ssmval_bysite.keys())
    val_bysite is some of values of all mutation at that site
    '''
    def __init__(self):
        '''
        empty initialization -- actual initializing 
        done with add_mutation()
        '''
        self.val = Counter()
        self.ssms = defaultdict(list)
        self.ssmvals = defaultdict(dict)

    def add_mutation(self,ssm,val):
        site = int(re.sub('\D*(\d+)\D*',r'\1',ssm))
        self.val[site] += max([0,val])
        self.ssms[site].append(ssm)
        if ssm in self.ssmvals[site]:
            raise RuntimeError(f'Already added {ssm=}')
        self.ssmvals[site][ssm] = val

    def topsites(self,n=10,thresh=0):
        sites = list(self.val.keys())
        sorted_sites = sorted(sites,key=self.val.get,
                              reverse=True)
        if thresh:
            sorted_sites = [site for site in sorted_sites
                            if self.val[site]>=thresh]
            
        return sorted(sorted_sites[:n])

    def plotsites(self,sitelist):
        '''main plotting, but should be sandwiched between
        and introductory plt.figure()
        and a final plt.show()
        '''
    
        maxval = max(self.val[site] for site in sitelist)
        textspace = maxval*0.02
            
        for ndx,site in enumerate(sitelist):
            maxbot = bot = 0
            ssms = sorted(self.ssms[site],
                          key=self.ssmvals[site].get,
                          reverse=True)
            ssmstr = ssmlabel(ssms)
            for n,ssm in enumerate(ssms):
                val = self.ssmvals[site][ssm]
                if val<0 and bot>0:
                    v.vprint(f'{site=}, {ssm=}, {val=}, {bot=}')
                    bot=0
                plt.bar(ndx,val,bottom=bot,color=nth_color(n))
                bot = bot+val
                maxbot = max([maxbot,bot])
            plt.text(ndx,maxbot+textspace,ssmstr,
                 horizontalalignment='center',
                 rotation='vertical',
                 fontsize='x-small')

        plt.xlabel('site')
        nsites = len(sitelist)
        nmodulo = 1 + nsites // 120
        tickndx = [nmodulo*i for i in range(nsites)
                   if nmodulo*i < nsites]
        plt.xticks(tickndx,
                   [str(sitelist[ndx]) for ndx in tickndx],
                   rotation = 'vertical')
        plt.xlim([-1,len(sitelist)])
        plt.gca().spines[['right', 'top']].set_visible(False)
