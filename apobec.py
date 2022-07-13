'''utilities for apobec analysis'''


## taken form biomaniac/kutil.py
_rctable = str.maketrans("ACTG","TGAC")
def reverse_complement(seq):
    '''return reverse complement of sequence (eg, AAACCT -> AGGTTT)'''
    ## note that dash (or any non ACTG char) is translated as is
    return seq.translate(_rctable)[::-1]

## default is strict
APOBEC_FD = set(['G'+R+D for R in "AG" for D in "AGT"])
APOBEC_RC = set([H+Y+'C' for Y in "CT" for H in "ACT"])
APOBEC = APOBEC_FD | APOBEC_RC

def loosen_apobec_rules():
    '''if called, re-defines (loosens) global APOBEC criteria'''
    global APOBEC_FD,APOBEC_RC,APOBEC
    APOBEC_FD = set(['G'+R+D for R in "AG" for D in "ACGT"])
    APOBEC_RC = set([H+Y+'C' for Y in "CT" for H in "ACGT"])
    APOBEC = APOBEC_FD | APOBEC_RC

def is_apobec(triplet,fwd='F'):
    '''F = forward, R = reverse complement, B = both'''
    if fwd == 'F':
        return triplet in APOBEC_FD
    if fwd == 'R':
        return triplet in APOBEC_RC
    if fwd in 'B':
        return triplet in APOBEC
    raise RuntimeError('fwd argument should be F, R, or B')

def sites_match_char(seq,c):
    '''find the indices of seq for which seq[n]==c'''
    n=0
    while True:
        try:
            n = seq.index(c,n)
            yield n
            n = n+1
        except ValueError:
            break

def get_all_mutsites(rseq,seq):
    '''
    return a list of every site that is any kind of mutation
    (well, not ones involving N or dash or other ambiguity codes)
    '''
    return [n
            for n,(rc,c) in enumerate(zip(rseq,seq))
            if (rc != c and rc in "ACGT" and c in "ACGT")]

        
#def get_mutation_ref(n,rseq,seq):
def get_mut_context(n,rseq,seq):
    '''
    for a given mutation at site n, 
    specify the three-character context string
    that it is part of the reference sequence
    '''
    if (rseq[n],seq[n]) == ("G","A"):
        return rseq[n:n+3]
    if (rseq[n],seq[n]) == ("C","T"):
        return rseq[n:n-3:-1] 
    return "xxx"
