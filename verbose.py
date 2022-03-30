'''verbose: print statements that depend on a user-set level of verbosity'''
## an alternative implementation using class instead of just module
import sys
from collections import Counter

class verbose:
    '''
    encapsulates vprint functions
    that write messages based on user-set level of verbosity
    Typical use:
       from verbose import verbose as v
       ...
       v.verbosity(1)
       ...
       v.vprint("read",n_sequences,"seqeunces from file")
       ...
    As an alternative to:
       verbosity_level=1
       ...
       if verbosity_level > 0:
           print("read",n_sequences,"sequences from file")
       ...
    '''

    ## class variables
    verbose_level=0
    count=Counter()

    @classmethod
    def verbosity(cls,vlevel):
        '''user sets the level of verbosity'''
        cls.verbose_level = vlevel

    @classmethod
    def vnprint(cls,vlevel,*p,**kw):
        '''
        if verbosity level is n or higher, then print
        to sys.stderr and flush every print
        '''
        if cls.verbose_level >= vlevel:
            print(*p,file=sys.stderr,flush=True,**kw)

    @classmethod
    def print(cls,*p,**kw):
        '''print (to stderr)'''
        cls.vnprint(0,*p,**kw)

    @classmethod
    def vprint(cls,*p,**kw):
        '''verbose print'''
        cls.vnprint(1,*p,**kw)

    @classmethod
    def vvprint(cls,*p,**kw):
        '''very verbose print'''
        cls.vnprint(2,*p,**kw)

    ## Variant of vprint that only prints the warning a fixed
    ## maximum number of times; thus warning you that something
    ## is happening, but not scrolling pages and pages of it by
    ## your screen.  Later, ie outside the loop you can (optionally)
    ## call for a summary of how many times the warning was triggered.

    ## common usage might be:
    ## from verbose import verbose as v
    ## ...
    ## for data_sample in data:
    ##    if is_bad(data_sample):
    ##         v.vprint_only(5,'Bad data',data_sample.info)
    ##    ...
    ## v.vprint_only_summary('Bad data')

    ## Note that the warning string is a key.  That means you
    ## can have several of these going at the saime time, identifying
    ## different potential causes for a warning.  But it also means
    ## that you cannot put changing information into that key string

    ## Okay:
    ## for data_sample in data:
    ##     if is_bad(data_sample):
    ##         v.vprint_only(5,'Bad data',data_sample.info)
    ##     if too_big(data_sample):
    ##         v.vprint_only(5,'Too big','Data sample has',data_sample.size,'units'))
    ##    ...
    ## v.vprint_only_summary('Bad data')
    ## v.vprint_only_summary('Too big')

    ## Mistake:
    ##    vprint_only(5,f'Bad data: {data_sample.info}')
    ##

    @classmethod
    def vnprint_only_summary(cls,vlevel,msg,*p,**kw):
        '''
        summarize how many times vprint_only was called
        with msg as first arguemnt; keep quiet if count was zero
        '''
        if cls.count[msg]:
            cls.vnprint(vlevel,msg,cls.count[msg],*p,**kw)

    @classmethod
    def vnprint_only(cls,vlevel,maxcount,msg,*p,**kw):
        '''vprint, but only maxcount times, at most'''
        cls.count[msg] += 1
        if cls.count[msg] <= maxcount:
            cls.vnprint(vlevel,msg,*p,*kw)

    @classmethod
    def print_only(cls,maxcount,msg,*p,**kw):
        '''vprint, but only maxcount times, at most'''
        cls.vnprint_only(0,maxcount,msg,*p,**kw)

    @classmethod
    def print_only_summary(cls,msg,*p,**kw):
        '''
        summarize how many times vprint_only was called
        with msg as first arguemnt
        '''
        cls.vnprint_only_summary(0,msg,*p,**kw)

    @classmethod
    def vprint_only(cls,maxcount,msg,*p,**kw):
        '''vprint, but only maxcount times, at most'''
        cls.vnprint_only(1,maxcount,msg,*p,**kw)

    @classmethod
    def vprint_only_summary(cls,msg,*p,**kw):
        '''
        summarize how many times vprint_only was called
        with msg as first arguemnt
        '''
        cls.vnprint_only_summary(1,msg,*p,**kw)

    @classmethod
    def vvprint_only(cls,maxcount,msg,*p,**kw):
        '''vprint, but only maxcount times, at most'''
        cls.vnprint_only(2,maxcount,msg,*p,**kw)

    @classmethod
    def vvprint_only_summary(cls,msg,*p,**kw):
        '''
        summarize how many times vprint_only was called
        with msg as first arguemnt
        '''
        cls.vnprint_only_summary(2,msg,*p,**kw)
