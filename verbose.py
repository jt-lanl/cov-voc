'''verbose: print statements that depend on a user-set level of verbosity'''

import sys
from collections import Counter

## Global variable (well, global to the verbose module)
VERBOSE=0

def verbosity(vlevel):
    '''user sets the level of verbosity'''
    global VERBOSE
    VERBOSE = vlevel

def vnprint(vlevel,*p,**kw):
    '''
    if verbosity level is n or higher, then print
    to sys.stderr and flush every print
    '''
    if VERBOSE >= vlevel:
        print(*p,file=sys.stderr,flush=True,**kw)

def vprint(*p,**kw):  vnprint(1,*p,**kw)
def vvprint(*p,**kw): vnprint(2,*p,**kw)

## Variant of vprint that only prints the warning a fixed
## maximum number of times; thus warning you that something
## is happening, but not scrolling pages and pages of it by
## your screen.  Later, ie outside the loop you can (optionally)
## call for a summary of how many times the warning was triggered.

## common usage might be:
## for data_sample in data:
##    if is_bad(data_sample):
##         vprint_only(5,'Bad data',data_sample.info)
##    ...
## vprint_only_summary('Bad data')

## Note that the warning string is a key.  That means you
## can have several of these going at the saime time, identifying
## different potential causes for a warning.  But it also means
## that you cannot put changing information into that key string

## Okay:
## for data_sample in data:
##     if is_bad(data_sample):
##         vprint_only(5,'Bad data',data_sample.info)
##     if too_big(data_sample):
##         vprint_only(5,'Too big','Data sample has',data_sample.size,'units'))
##    ...
## vprint_only_summary('Bad data')
## vprint_only_summary('Too big')

## Mistake:
##    vprint_only(5,f'Bad data: {data_sample.info}')
##

def vprint_only_summary(msg):
    '''summarize how many times vprint_only was called with msg as first arguemnt'''
    vprint_only(None,msg)

def vprint_only(maxcount,msg,*p,**kw):
    '''vprint, but only maxcount times, at most'''
    try:
        vprint_only.count[msg] += 1
    except AttributeError:
        vprint_only.count = Counter()
        vprint_only.count[msg] += 1

    if maxcount is None:
        vprint_only.count[msg] -= 1 # to un-count the call with None
        if vprint_only.count[msg] > 0:
            vprint(msg,vprint_only.count[msg],"warnings")
    elif vprint_only.count[msg] < maxcount:
        vprint(msg,*p,*kw)

## should there also be a vvprint_only ???
