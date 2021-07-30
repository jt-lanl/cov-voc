import sys
import re
import datetime
import numpy as np

import readseq

def copy_seqlist(seqlist):
    ## nb, only copies .name and .seq attributes
    return [readseq.SequenceSample(s.name,s.seq) for s in seqlist]

def gen_columns_seqlist(seqlist):
    '''generator that produces columns,
       where each column is a tuple of aa's at a given site.
       note, this requires that seqlist be an alignment.  equivalent to:
       ( [s.seq[n] for s in seqlist] for n in range(len(seqlist[0].seq)) )'''
    return zip(*[s.seq for s in seqlist])

def getcolumn(seqs,n,keepx=False):
    ''' return list of aa's at given site number (column) Equiv: multicolumn(seqs,[n],keepx=keepx)'''
    ## note, site number is one + zero-based index; so n-1 is index of s.seq
    if keepx:
        return [s.seq[n-1] for s in seqs]
    else:
        return [s.seq[n-1] for s in seqs if "X" not in s.seq[n]]

def multicolumn(seqs,nlist,keepx=False):
    ''' return a list of aa strings, corresponding to the columns indicated by nlist '''
    aaslist = [ "".join(s.seq[n-1] for n in nlist) for s in seqs ]
    if not keepx:
        aaslist = [aas for aas in aaslist if "X" not in aas]
    return aaslist

def numpy_from_seqlist(seqlist):
    ''' Experimental!'''
    x = np.array([list(s.seq) for s in seqlist])
    ## next line 1byte/character, use chr() to get characters back
    #x = np.asarray(np.vectorize(ord)(x),dtype=np.byte) 
    return x

    
def relativename(master,mutant,matchchar="."):
    '''Mutant with blanks(matchchar's) for matches against a master;
    eg, master="ABC", mutant="ABD", then output="..D"
    '''
    s=""
    for a,b in zip(master,mutant):
        s += matchchar if a==b else b
    return s

def str_indexes(s,c):
    '''similar to s.index(c) but finds /all/ indexes n such that s[n]==c; 
    return list of n's '''
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
    
def stripdashcols(master,seqs,dashchar="-"):
    '''strips positions from each sequence in seqs array, 
    based on dashes in master sequence'''
    ndx = str_indexes(master,dashchar)
    keep = [n for n in range(len(master)) if n not in ndx]
    for s in seqs:
        s.seq = "".join(s.seq[n] for n in keep)

## date-based utilities

def date_fromiso(s):
    if type(s) == datetime.date:
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError:
        if s == ".":
            return None
        return None #raise RuntimeError(f"Invalid Date {s}")

def date_from_seqname(s):
    datestring = re.sub(".*(\d\d\d\d-\d\d-\d\d).*",r"\1",s.name)
    return date_fromiso(datestring)
    
def add_date_attribute(seqlist):
    for s in seqlist:
        s.date = date_from_seqname(s.name)

def count_bad_dates(seqlist):
    return(sum(date_from_seqname(s) is None for s in seqlist))

def range_of_dates(seqlist):
    dates = [date_from_seqname(s) for s in seqlist]
    dates = [d for d in dates if d is not None]
    return min(dates).isoformat(),max(dates).isoformat()

def filter_by_date(seqs,fromdate,todate,keepfirst=False):
    '''input is iterable (list or iterator); output is generator'''

    f_date = date_fromiso(fromdate)
    t_date = date_fromiso(todate)

    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        d = date_from_seqname(s)
        if not d:
            continue
        if f_date and f_date > d:
            continue
        if t_date and t_date < d:
            continue
        yield s

def filter_by_pattern(seqs,pattern,keepfirst=False,ignorecase=True):

    if "Global" == pattern: ## should test in calling routine
        yield from seqs
    
    flags = re.I if ignorecase else 0        
    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        if re.search(pattern,s.name,flags):
            yield s


def filter_by_patternlist(seqs,patternlist,
                          keepfirst=False,ignorecase=True):

    if "Global" in patternlist: ## should test in calling routine
        yield from seqs
    
    flags = re.I if ignorecase else 0        
    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        if any(re.search(pattern,s.name) for pattern in patternlist):
            yield s

    
def filter_by_pattern_exclude(seqs,pattern,keepfirst=False):
    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        if not re.search(pattern,s.name):
            yield s

def filter_by_patternlist_exclude(seqs,patternlist,keepfirst=False):
    for n,s in enumerate(seqs):
        if keepfirst and n == 0:
            yield s
            continue
        if not any(re.search(pattern,s.name) for pattern in patternlist):
            yield s
            
def mutantlist(reference,variant,returnstring=False,badchar=None):
    '''return a list of mutations, of the form "AnB", where ref[n-1]=A and var[n-1]=B and n is site number'''
    assert( len(reference) == len(variant) )
    mutants = []
    for i,(ref,var) in enumerate(zip(reference,variant),start=1):
        if ref != var and var != badchar:
            mutants.append( ref + str(i) + var )
    if returnstring:
        mutants = "["+",".join(mutants)+"]"
    return mutants
    

if __name__ == "__main__":

    import argparse
    def getargs():
        ap = argparse.ArgumentParser()
        paa = ap.add_argument
        paa("--input","-i",
            help="input fasta file")
        paa("--range",nargs=2,
            help="range of two dates, in yyyy-mm-dd format")
        paa("--verbose","-v",action="count",default=0,
            help="verbose")
        args = ap.parse_args()
        return args

    args = getargs()
    def vprint(*p,**kw):
        if args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        if args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)


    if args.input:
        seqlist = readseq.read_seqfile(args.input)
        print("seqlist:",type(seqlist))
        seqlist = list(seqlist)
        print("seqlist:",type(seqlist))
        print("Sequences:",len(seqlist))

        print("Bad dates:",count_bad_dates(seqlist))

        if args.range:
            seqlist = filter_by_date(seqlist,args.range[0],args.range[1],
                                     keepfirst=True)
            print("Sequences:",len(seqlist)-1,"in date range:",args.range)

            
        


    

