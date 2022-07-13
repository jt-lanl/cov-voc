'''Plot locations of mutations colored by APOBEC pattern

Each sequence gets a horizontal line, corrsponding to position
Forward G->A mutations are indicated with an above-the-line tick mark
Reverse-complement C->T mutations use below-the-line tick marks
All other mutations are indicated with smaller gray tick marks
Blue/Cyan indicates Primary/Secondary style of apobec mutations
Red indicates non-apobec mutations (by loose criterion)
Magenta indicates mutations that would be considered apobec by the
loose rules, but non-apobec by the strict rules.

...pairs.py: reads in pairs of sequence names and then makes
plots according to the rseq/seq pairs

'''

import argparse
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

import verbose as v
import sequtil
import apobec as apo

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input sequence file")
    paa("--nseq",type=int,default=0,
        help="if specified, then only plot this many sequences")
    paa("--reversecomplement","-r",action="store_true",
        help="reverse complement the strings as they are read in")
    paa("--strict",action="store_true",
        help="Use stricter APOBEC pattern definition")
    paa("--loose",action="store_false",dest='strict',
        help="Use looser APOBEC pattern definition")
    paa("--merge",action="store_true",
        help="Add an extra merge sequence")
    paa("--pairs",
        help="file with pairs of sequence names")
    paa("--output","-o",
        help="write plot to output file")
    paa("--showgaps",action="store_true",
        help="show where the reference sequence gaps are in the alignment")
    paa("--nosnp",
        help="write sequences names with no SNPs into file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args.loose = not args.strict
    return args

def get_mutation_ref(n,rseq,seq):
    '''
    for a given mutation at site n, specify the three-character ref string
    that it is part of
    '''
    if (rseq[n],seq[n]) == ("G","A"):
        return rseq[n:n+3]
    if (rseq[n],seq[n]) == ("C","T"):
        return rseq[n:n-3:-1]
    return "xxx"

class LineStyle:
    def __init__(self,color,style='solid',ylo=-0.3,yhi=0.3):
        self.color = color
        self.style = style
        self.ylo = ylo
        self.yhi = yhi
    def copy(self):
        return LineStyle(self.color,self.style,self.ylo,self.yhi)

    def __str__(self):
        return f'({self.color},{self.style}) [{self.ylo},{self.yhi}]'

_colorstyle_of_mutation = {
    ## two styles of forward apobec mutaions
    'GA' : LineStyle('Blue',ylo=0),
    'GG' : LineStyle('Cyan',ylo=0),
    'GC' : LineStyle('#BB0000',ylo=0),
    'GT' : LineStyle('#BB0000',ylo=0),
    ## two styles of backward apobec mutations
    'CT' : LineStyle('Blue',yhi=0),
    'CC' : LineStyle('Cyan',yhi=0),
    'CG' : LineStyle('#BB0000',yhi=0),
    'CA' : LineStyle('#BB0000',yhi=0),
    ##
    'other' : LineStyle('Gray',ylo=-0.1,yhi=0.1)
    }

def is_mut_ctx_apobec(mut_ctx):
    if mut_ctx[0] == 'C':
        mut_ctx = mut_ctx[::-1]
    return apo.is_apobec(mut_ctx,fwd='B')

def type_to_color(mut_ctx,loose=True):
    '''
    specify color (and other linestyle properties) for given mutation type
    '''
    mut_type = mut_ctx[:2]
    if mut_type not in _colorstyle_of_mutation:
        return _colorstyle_of_mutation['other'].copy()
    mut_linestyle =  _colorstyle_of_mutation[mut_type].copy()
    ## if not apobec, then re-color it red
    if not is_mut_ctx_apobec(mut_ctx) and mut_linestyle.color != '#BB0000':
        mut_linestyle.color='Magenta'
    return mut_linestyle

def get_figure_size(nitems):
    figsize = (12,1+0.75*nitems)
    if figsize[1] > 9:
        figsize = (12,1+0.4*nitems)
    return figsize

def main_pairs(args,seqs):
    '''if arg.pairs specified'''

    seqs_byname = {s.name:s for s in seqs}

    pairs = []
    with open(args.pairs,'r') as fpairs:
        for line in fpairs:
            line=line.strip()
            v.vprint('line=',line)
            try:
                firstname,seqname = line.strip().split()
                pairs.append( (seqs_byname[firstname],
                               seqs_byname[seqname]) )
            except ValueError:
                continue
    v.vprint('pairs:',[(a.name,b.name) for a,b in pairs])

    plt.figure(figsize=get_figure_size(len(pairs)))

    pairs = pairs[::-1] ## reverse order so first is on top (ie, plotted last)
    for ns,(first,s) in enumerate(pairs):
        v.vprint('ns=',ns,first.name,s.name)
        plt.plot([0,len(s.seq)],[ns,ns],'k-')
        mut_sites = apo.get_all_mutsites(first.seq,s.seq)
        for n in mut_sites:
            site_ref = get_mutation_ref(n,first.seq,s.seq)
            ls = type_to_color(site_ref)
            plt.plot([n,n],[ns+ls.ylo,ns+ls.yhi],color=ls.color,linestyle=ls.style,lw=1)

    plt.yticks(range(len(pairs)),
               labels=[first.name + "/" + s.name for (first,s) in pairs],
               fontsize='x-small')

def re_order_sequences(first,seqs,nosnpfile=None):
    '''re-order sequences so that identical sequences are listed
    last (ie, plotted at the bottom); remaining sequences are
    reversed'''
    seqs_reordered = []
    seqs_identical = []
    for s in seqs:
        mut_sites = apo.get_all_mutsites(first.seq,s.seq)
        if len(mut_sites) > 0:
            seqs_reordered.append(s)
        else:
            seqs_identical.append(s)

    if nosnpfile:
        with open(nosnpfile,'w') as fnosnp:
            for s in seqs_identical:
                print(s.name,file=fnosnp)
        n_identical = len(seqs_identical)
        seqs_identical = seqs_identical[:1]
        seqs_identical[0].name = f'{n_identical} sequences with no SNPs'

    return seqs_reordered + seqs_identical


def main_nopairs(args,seqs):
    '''if args.pairs is not called, then run first against the rest'''

    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    seqs = re_order_sequences(first,seqs,args.nosnp)

    if args.nseq:
        seqs = list(seqs)[:args.nseq]

    ## get the locations of all the -'s in the first seqence
    r_dash_indices = set(sequtil.str_indexes(first.seq,'-'))

    count_types = Counter()
    plt.figure(figsize=get_figure_size(len(seqs)))
    seqs = seqs[::-1]
    all_mut_sites = set()
    color_type = dict()
    for ns,s in enumerate(seqs):
        v.vprint('seq: ',ns,s.name)
        if args.showgaps:
            s_dash_indices = set(sequtil.str_indexes(s.seq,'-'))
            s_n_indices = set(sequtil.str_indexes(s.seq,'N'))
            s_minus_r = s_dash_indices - r_dash_indices
            r_minus_s = r_dash_indices - s_dash_indices
            for din in s_minus_r:
                ## Green deletions
                plt.plot([din,din],[ns-0.4,ns+0.4],'-',color='#99FF99')
            for din in r_minus_s - s_n_indices:
                ## Blue insertions
                plt.plot([din,din],[ns-0.4,ns+0.4],'-',color='#BBBBFF')
            for din in s_n_indices:
                ## Magenta N's
                plt.plot([din,din],[ns-0.2,ns+0.2],'-',color='#FF88FF')
                              
        plt.plot([0,len(s.seq)],[ns,ns],'k-')
        mut_sites = apo.get_all_mutsites(first.seq,s.seq)
        all_mut_sites.update(mut_sites)
        v.vprint(s.name,end=':')
        for n in mut_sites:
            site_ref = get_mutation_ref(n,first.seq,s.seq)
            count_types[site_ref] += 1
            ls = type_to_color(site_ref)
            color_type[n] = ls
            v.vprint("<",n,site_ref,ls.color,">",end=' ')
            plt.plot([n,n],[ns+ls.ylo,ns+ls.yhi],
                     color=ls.color,linestyle=ls.style,lw=2)

        v.vprint()
    if args.merge:
        ns += 1
        seqs.append(seqs[-1].copy())
        seqs[-1].name="Merged SNPs"
        plt.plot([0,len(seqs[-1].seq)],[ns,ns],'k-')
        for n in all_mut_sites:
            ls = color_type[n]
            plt.plot([n,n],[ns+ls.ylo,ns+ls.yhi],
                            color=ls.color,linestyle=ls.style,lw=2)

    for t,cnt in count_types.items():
        v.vprint(f'{cnt:>6d} {t} ({type_to_color(t)})')


    plt.yticks(range(len(seqs)),
               labels=[s.name for s in seqs])
    plt.ylim([-1,len(seqs)])


def _main(args):
    '''main'''
    v.vprint(args)

    if args.loose:
        apo.loosen_apobec_rules()

    v.vprint("APOBEC_FD:",apo.APOBEC_FD)

    seqs = sequtil.read_seqfile(args.input)
    seqs = list(seqs)

    if args.reversecomplement:
        for s in seqs:
            s.seq = apo.reverse_complement(s.seq)

    if args.pairs:
        main_pairs(args,seqs)
    else:
        main_nopairs(args,seqs)

    plt.xlabel('Position in sequence')
    plt.tight_layout()

    if args.output:
        plt.savefig(args.output)
    else:
        plt.show()

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
