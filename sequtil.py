'''utilities for seqs;
   where seqs is a list (or generator) of SequenceSample's
'''

import re
import itertools
from collections import Counter
import numpy as np
from scipy import stats

import verbose as v
from seqsample import SequenceSample
from readseq import read_seqfile,write_seqfile
import intlist

def copy_seqlist(seqlist):
    '''return a 'deep' copy of the sequence list'''
    ## nb, only copies .name and .seq attributes
    return [SequenceSample(s.name,s.seq) for s in seqlist]

def gen_columns_seqlist(seqlist):
    '''generator that produces columns,
       where each column is a tuple of aa's at a given site.
       note, this requires that seqlist be an alignment.  equivalent to:
       ( [s.seq[n] for s in seqlist] for n in range(len(seqlist[0].seq)) )'''
    return zip(*[s.seq for s in seqlist])

def getcolumn(seqs,n,keepx=False):
    '''
    return list of aa's at given column index
    Equiv: multicolumn(seqs,[n],keepx=keepx)
    '''
    ## note, site number is one + zero-based index; so n+1 is site number
    if keepx:
        return [s.seq[n] for s in seqs]
    return [s.seq[n] for s in seqs if "X" not in s.seq[n]]

def multicolumn(seqs,nlist,keepx=False):
    ''' return a list of aa strings, corresponding to the columns indicated by nlist '''
    aaslist = [ "".join(s.seq[n] for n in nlist) for s in seqs ]
    if not keepx:
        aaslist = [aas for aas in aaslist if "X" not in aas]
    return aaslist


def numpy_from_seqlist(seqlist):
    ''' Experimental!'''
    x = np.array([list(s.seq) for s in seqlist])
    ## next line 1byte/character, use chr() to get characters back
    #x = np.asarray(np.vectorize(ord)(x),dtype=np.byte)
    return x

def xentropy(clist,keepx=False):
    '''entropy from a list of characters'''
    cnt = Counter(clist)
    if not keepx:
        cnt.pop('X',None)
    return stats.entropy(list(cnt.values()))

def zentropy(clist,wlist=None,keepx=False):
    '''entropy from a list of characters and a list of weights'''
    cnt = Counter()
    wlist = wlist or [1]*len(clist)
    for c,w in zip(clist,wlist):
        cnt[c] += w
    if not keepx:
        cnt.pop('X',None)
    return stats.entropy(list(cnt.values()))

def chunked_entropy(seqs,chunk=500,keepx=False):
    '''compute entropy one chunk at a time'''
    ## a chunk is /all/ the sequences, and /some/ of the sites
    ## use Counter() the reduce list of subseq's to
    ## a shorter list of unique subseq's
    E = []
    L = len(seqs[0].seq)
    for k in range(0,L,chunk):
        subseqs = [s.seq[k:k+chunk] for s in seqs]
        cnt = Counter(subseqs)
        substrs,wts = zip(*cnt.items())
        E.extend(zentropy(chars,wts,keepx=keepx)
                 for chars in zip(*substrs))
    return E

def mostcommonchar(clist):
    '''return the most common item in the list'''
    [(c,_)] = Counter(clist).most_common(1)
    return c

def consensus(seqlist):
    '''create a consesnsus sequence from the sequence list'''
    return "".join(mostcommonchar(clist)
                   for clist in gen_columns_seqlist(seqlist))


def relativename(master,mutant,matchchar="."):
    '''Mutant with blanks(matchchar's) for matches against a master;
    eg, master="ABC", mutant="ABD", then output="..D"
    '''
    s=""
    for a,b in zip(master,mutant):
        s += matchchar if a==b else b
    return s

def get_first_item(items,keepfirst=True):
    '''
    get first item in an iterable, and return (item,iterable) as a tuple
    if keepfirst==True, then first item remains in iterable
    note that the iterable can be a list or an iterator
    returned iterable will be list or iterator, depending on input
    '''
    if isinstance(items,list):
        first = items[0]
        if not keepfirst:
            items = items[1:]
    else:
        first = next(items)
        if keepfirst:
            items = itertools.chain([first],items)
    return first,items

def checkseqlengths(seqs):
    '''ensure that all sequences are of the same length'''
    first,seqs = get_first_item(seqs)
    seqlen = len(first.seq)
    for s in seqs:
        if len(s.seq) != seqlen:
            raise RuntimeError(f"Sequence {s.name} has inconsistent length: "
                               f"{len(s.seq)} vs {seqlen}")
        yield s

def str_indexes(s,c):
    '''
    similar to s.index(c) but returns list of /all/ indexes n such that s[n]==c;
    '''
    ## equiv one-liner: [m.span(0)[0] for m in re.finditer(c,s)]
    ## and yes i know the plural of index is indices
    ndx = []
    n=0
    while True:
        try:
            n = s.index(c,n)
            ndx.append(n)
            n += 1
        except ValueError:
            break
    return ndx

## stripdashcols: strips columns where dash appears in the master sequence
## remove_gap_columns: strips columns where dash appears in all sequences

def stripdashcols(master,seqs,dashchar="-"):
    '''strips positions from each sequence in seqs array,
    based on dashes in master sequence'''
    if dashchar not in master:
        yield from seqs
    ndx = str_indexes(master,dashchar)
    keep = [n for n in range(len(master)) if n not in ndx]
    for s in seqs:
        s.seq = "".join(s.seq[n] for n in keep)
        yield s

def get_gap_columns(first,seqs,dashchar='-'):
    '''find indexes associated with gap-only columns'''
    ## ie, indices such that s.seq[n]==dashchar for all s in seqs
    dashindexes = str_indexes(first.seq,dashchar)
    distinct_sequences = set(s.seq for s in seqs)
    for seq in distinct_sequences:
        dashindexes = [n for n in dashindexes if seq[n]==dashchar]
        if not dashindexes:
            break
    return dashindexes

def remove_gap_columns(first,seqs,dashchar='-'):
    '''remove the gap-only columns'''
    seqs = list(seqs) ## to ensure get_gap_column() doesn't consume seqs
    dashindexes = get_gap_columns(first,seqs,dashchar=dashchar)
    keepndx = [n for n in range(len(first.seq)) if n not in dashindexes]
    keepranges = intlist.intlist_to_rangelist(keepndx)
    v.vprint("ndx:",len(keepndx))
    v.vprint("rng:",len(keepranges))
    first.seq = "".join(first.seq[lo:hi] for lo,hi in keepranges)
    for s in seqs:
        #s.seq = "".join(s.seq[n] for n in keepndx)
        s.seq = "".join(s.seq[lo:hi] for lo,hi in keepranges)
    return first,seqs



######################## Filter based on pattens in the name of seq

def filter_by_patternlist(seqs,patternlist,
                          exclude=False,
                          ignorecase=True):
    '''
    return an iterator of SequenceSample's whose names
    match any of the patterns in the patternlist
    '''
    if not patternlist:
        ## empty list means no filtering
        yield from seqs

    flags = re.I if ignorecase else 0
    for s in seqs:
        ## if exclude is False, then yield if any matches
        ## if exclude is True, then yield if not any matches
        anymatches = any(re.search(pattern,s.name,flags)
                         for pattern in patternlist)
        if bool(exclude) ^ bool(anymatches):
            yield s

def filter_by_patternlist_exclude(seqs,patternlist,**kwargs):
    '''
    return an iterator of SequenceSample's that do NOT match
    any of the patterns in the patternlist
    '''
    return filter_by_patternlist(seqs,patternlist,exclude=True,**kwargs)

def filter_by_pattern(seqs,pattern,**kwargs):
    '''
    return an iterator of SequenceSample's that match the pattern
    '''
    return filter_by_patternlist(seqs,[pattern],**kwargs)

def filter_by_pattern_exclude(seqs,pattern,**kwargs):
    '''
    return an iterator of SequenceSample's that do not match the pattern
    '''
    return filter_by_patternlist_exclude(seqs,[pattern],**kwargs)

if __name__ == "__main__":

    import argparse
    import readseq

    def getargs():
        '''get arguments from command line'''
        ap = argparse.ArgumentParser()
        paa = ap.add_argument
        paa("--input","-i",
            help="input fasta file")
        paa("--verbose","-v",action="count",default=0,
            help="verbose")
        args = ap.parse_args()
        return args

    xargs = getargs()
    v.verbosity(xargs.verbose)

    if xargs.input:
        sequences = readseq.read_seqfile(xargs.input)
        print("sequences:",type(sequences))
        sequences = list(sequences)
        print("sequences:",type(sequences))
        print("Sequences:",len(sequences))
