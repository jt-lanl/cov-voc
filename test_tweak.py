import pytest

from mutant import Mutation, MutationManager
import tweak as tku


ref = 'MFVFL-VLLPLVSSQCV-------N-L----T-TRT-'
seqs = '''
MFVFL-VLLPLVSSQCV-------N-L----T-TRT-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCV-------N-L----T-TRT-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------N-L----R-TTT-
MFVFL-VLLPLVSSQCVMPLF----------T-T-T-
MFVFL-VLLPLVSSQCV-------MPLFNFIT-T-T-
MFVFL-VLLPLVSSQCV-------N-FR---T-T-T-
'''
mstringpairs='''
[N17M,+17P,L18L,+18F] [+16MPLF,N17-,L18-]
[+18I,T19T,T20T,R21-] [T19I,T20T,R21T]
[+18R,T19T,T20T,R21-] [T19R,T20T,R21T]
'''

seqs = [seq for seq in seqs.split('\n') if seq]

mstringpairs = [tuple(mstringpair.split())
                for mstringpair in mstringpairs.split('\n')
                if mstringpair]

mut_mgr = MutationManager(ref)

def test_extra_dashes():
    xpander = tku.ExpandSeq.from_mstringpairs(mstringpairs)
    xtra = xpander.needed
    assert xtra == {17: 1, 18: 1, 16: 4}


@pytest.mark.filterwarnings("ignore")    
def test_tweaks_from_identical_mstringpairs():
    for ma,mb in mstringpairs:
        ## mstring -> itself should give tweak==None
        assert None == tku.tweak_from_mstringpair(mut_mgr,ma,ma)
        assert None == tku.tweak_from_mstringpair(mut_mgr,mb,mb)

def test_tweaks_from_mstringpairs():
    ma = '[+18I,T19T,T20T,R21-]'
    mb = '[T19I,T20T,R21T]'
    tweak = tku.tweak_from_mstringpair(mut_mgr,ma,mb)
    print(ma,mb,tweak.ndxlo,tweak.sa,tweak.sb)
    assert tweak.ndxlo == 27
    assert tweak.sa == "I---T-T-"
    assert tweak.sb == "----I-TT"


def test_expand():
    mstringpairs = [("[N17M,+17P,L18L,+18F]","[+16MPLF,N17-,L18-]"),
                    ("[+18I,T19T,T20T,R21-]","[T19I,T20T,R21T]")]
    ref = 'MFVFL-VLLPLVSSQCV--N-LT-TRT-'
    
    xpand = tku.ExpandSeq.from_mstringpairs(mstringpairs,ref)
    print('needs=',xpand.needed)
    print('xpand=',xpand.expand)
    assert xpand.needed == {17: 1, 18: 1, 16: 4}
    assert xpand.expand == {21: 1, 18: 2}

    ## add two dashes after "V--", one dash after "L"
    #ref = 'MFVFL-VLLPLVSSQCV--N-LT-TRT-'
    xref = 'MFVFL-VLLPLVSSQCV----N-L-T-TRT-'
    assert xref == xpand.expand_seq(ref)

    ## add two dashes after "V--", one dash after "L"
    seq =  'MFVFL-VLLPLVSSQCV--N-LR-TTT-'
    xseq = 'MFVFL-VLLPLVSSQCV----N-L-R-TTT-'
    assert xseq == xpand.expand_seq(seq)

    ## add two dashes after "VTT", one dash after "L"
    seq =  'MFVFL-VLLPLVSSQCVTTN-LR-TTT-'
    xseq = 'MFVFL-VLLPLVSSQCVTT--N-L-R-TTT-'
    assert xseq == xpand.expand_seq(seq)

    ## Let's rebuild xpand in a different way
    xpand = tku.ExpandSeq.from_mstringpairs(mstringpairs)
    try:
        ## even though xpand is not fully initialized
        xref = xpand.expand_seq(ref)
        ## if exception not raised, then that's an error!
        assert False
    except RuntimeError:
        ## ought to raise an exception!
        assert True

    ## Let's finish initializing

    xpand.update_refseq(ref)
    ## add two dashes after "VTT", one dash after "L"
    seq =  'MFVFL-VLLPLVSSQCVTTN-LR-TTT-'
    xseq = 'MFVFL-VLLPLVSSQCVTT--N-L-R-TTT-'
    assert xseq == xpand.expand_seq(seq)
    

if __name__ == '__main__':
    for msp in mstringpairs:
        print('msp:',msp)

    test_extra_dashes()
    test_tweaks_from_identical_mstringpairs()
    test_tweaks_from_mstringpairs()
    test_expand()


    
