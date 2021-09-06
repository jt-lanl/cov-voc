'''needleman wunch alignment, with a window'''

import numpy as np

verbose=False
# the three directions you can go in the traceback:
DIAG = 0 
UP = 1 
LEFT = 2


def needleman_wunsch_matrix(seq1, seq2, win=0):
    """
    fill in the DP matrix according to the Needleman-Wunsch algorithm.
    Returns the matrix of scores and the matrix of pointers
    """
 
    match =  1  # match score
    mismatch = -1  # mismatch penalty
    indel = -1 # indel penalty
 
    n = len(seq1)
    m = len(seq2)
    s = np.zeros( (n+1, m+1) ) - 100 # DP matrix (large negative is default)
    s[0,0]=0
    ptr = np.zeros( (n+1, m+1), dtype=int  ) # matrix of pointers
 
    ##### INITIALIZE SCORING MATRIX (base case) #####
 
    for i in range(1, n+1) :
        s[i,0] = indel * i
    for j in range(1, m+1):
        s[0,j] = indel * j
 
    ########## INITIALIZE TRACEBACK MATRIX ##########
 
    # Tag first row by LEFT, indicating initial "-"s
    ptr[0,1:] = LEFT
 
    # Tag first column by UP, indicating initial "-"s
    ptr[1:,0] = UP
 
    #####################################################


    for i in range(1,n+1):
        jrange = range(max(1,i-win),min(i+win,m+1)) if win else range(1,m+1)
        for j in jrange:
            # match
            if seq1[i-1] == seq2[j-1]:
                s[i,j] = s[i-1,j-1] + match
                ptr[i,j] = DIAG
            # mismatch
            else :
                s[i,j] = s[i-1,j-1] + mismatch
                ptr[i,j] = DIAG
            # indel penalty
            if s[i-1,j] + indel > s[i,j] :
                s[i,j] = s[i-1,j] + indel
                ptr[i,j] = UP
            # indel penalty
            if s[i, j-1] + indel > s[i,j]:
                s[i,j] = s[i, j-1] + indel
                ptr[i,j] = LEFT
 
    return s, ptr
 
def needleman_wunsch_trace(seq1, seq2, s, ptr) :
 
    #### TRACE BEST PATH TO GET ALIGNMENT ####
    align1 = ""
    align2 = ""
    n, m = (len(seq1), len(seq2))
    i = n
    j = m
    curr = ptr[i, j]
    while (i > 0 or j > 0):        
        ptr[i,j] += 3
        if curr == DIAG :            
            align1 = seq1[i-1] + align1
            align2 = seq2[j-1] + align2
            i -= 1
            j -= 1            
        elif curr == LEFT:
            align1 = "-" + align1
            align2 = seq2[j-1] + align2
            j -= 1            
        elif curr == UP:
            align1 = seq1[i-1] + align1
            align2 = "-" + align2
            i -= 1
 
        curr = ptr[i,j]
 
    return align1, align2
 
def needleman_wunsch(seq1, seq2, win=0) :
    """
    computes an optimal global alignment of two sequences using the Needleman-Wunsch
    algorithm
    returns the alignment and its score
    """
    s,ptr = needleman_wunsch_matrix(seq1, seq2, win=win)
    align1,align2 = needleman_wunsch_trace(seq1, seq2, s, ptr)
 
    if verbose:
        show_dp_matrix(s, seq1, seq2)
        show_ptr_matrix(ptr, seq1, seq2)
        print("\n"+"~`"*25)
        
        print("Alignment Score: %f\n" % (s[len(seq1),len(seq2)]))
        print("Alignment:")
        print(alignment[0])
        print(alignment[1])
 
    return align1, align2, s[len(seq1), len(seq2)]

def score_alignment(alig1,alig2):
    match =  1  # match score
    mismatch = -1  # mismatch penalty
    indel = -1 # indel penalty
 
    score = 0
    for a,b in zip(alig1,alig2):
        if a=='-' and b=='-':
            score += 0 
        elif a=='-' or b=='-':
            score += indel
        elif a==b:
            score += match
        else:
            score += mismatch
    return score
