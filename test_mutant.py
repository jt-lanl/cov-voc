'''Module for testing mutant.py routines'''

import mutant

RefSeq = "ABC-D--EFG"
NewSeqList = [
    "BBC-E-WEFF",
    "BBC-EW-EFF",
    "BBC-E--EFF",
    "ABCWEXYEFG",
    "ABCWDXYEFG",
    "ABC-E--EFG",
    "ABCWD--EFH",
    "AB--D--EFG",
]
mut_dict = {
    'BBC-EW-EFF': '[A1B,D4E,+4W,G7F]!',
    'BBC-E--EFF': '[A1B,D4E,G7F]!',
    'ABCWEXYEFG': '[+3W,D4E,+4XY]!',
    'ABCWDXYEFG': '[+3W,+4XY]!',
    'ABC-E--EFG': '[D4E]!',
    'ABCWD--EFH': '[+3W,G7H]!',
    'AB--D--EFG': '[C3-]!',
}

mut_dict_xtra = {
    'BBC-E-WEFF': '[A1B,D4E,+4-W,G7F]!',
}

MM = mutant.MutationManager(RefSeq)

def test_ssm():
    ssmstring = '+4-W'
    ssm = mutant.SingleSiteMutation.from_string(ssmstring)
    assert str(ssm) == ssmstring
    
def test_hard_seq_to_mutations():
    for seq,mstring in mut_dict_xtra.items():
        #mut = MM.seq_to_mutation(seq)
        mut = MM.substr_to_mutation(seq)
        mut.exact = True
        assert mut.as_string() == mstring
    
def test_seq_to_mutations():
    for seq,mstring in mut_dict.items():
        mut = MM.seq_to_mutation(seq)
        assert mut.as_string() == mstring
    for seq,mstring in mut_dict.items():
        mut = MM.seq_to_mutation(seq)
        assert seq == MM.seq_from_mutation(mut)

    for seq,mstring in mut_dict.items():
        ## these are sequences which /cannot/ be reconstructed from mstring
        mut = MM.seq_to_mutation(seq)
        #assert seq != MM.seq_from_mutation(mut)

def seq_fits_mstring(mstring,seq,**kw):
    ## NEW way
    try:
        re_mpatt = MM.regex_from_mstring(mstring,compile=True,**kw)
        re_fits = bool(re_mpatt.match(seq))
    except ValueError:
        re_fits = False
    ## OLD way
    if False:
        mpatt = mutant.Mutation.from_mstring(mstring)
        fits = MM.seq_fits_pattern(mpatt,seq,**kw)
        if re_fits != fits:
            print(mstring,seq,MM.regex_from_mstring(mstring,**kw),'mpatt=',mpatt)
    return re_fits

def test_extended_mstrings():

    assert seq_fits_mstring('[A1B,D4E,+4W,G7*]', 'BBC-EW-EFF',extended=False) == False
    assert seq_fits_mstring('[A1BC,D4E,+4W,G7*]', 'BBC-EW-EFF',extended=False) == False

    def seq_fits_mstring_x(*args,**kw):
        return seq_fits_mstring(*args,**kw,extended=True)

    assert seq_fits_mstring_x('[A1B,D4E,+4W,G7F]', 'BBC-EW-EFF')
    assert seq_fits_mstring_x('[A1B,D4E,+4W,G7*]', 'BBC-EW-EFF')
    assert seq_fits_mstring_x('[A1BC,D4E,+4W,G7*]', 'BBC-EW-EFF')
    
    assert seq_fits_mstring_x('[D4*]', 'BBC-E-WEFF')            == True
    assert seq_fits_mstring_x('[D4*]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'BBC-EW-EFF')            == True
    assert seq_fits_mstring_x('[D4*]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'BBC-E--EFF')            == True
    assert seq_fits_mstring_x('[D4*]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring_x('[D4*]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring_x('[D4*]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'ABC-E--EFG')            == True
    assert seq_fits_mstring_x('[D4*]', 'ABC-E--EFG',exact=True) == True
    assert seq_fits_mstring_x('[D4*]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring_x('[D4*]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring_x('[D4*]', 'AB--D--EFG')            == False
    assert seq_fits_mstring_x('[D4*]', 'AB--D--EFG',exact=True) == False


def test_seq_to_mutation():
    for seq,mstring in mut_dict.items():
        mut = MM.seq_to_mutation(seq)
        assert mut.as_string() == mstring        
    

def test_hard_mstrings():

    seq = 'BBC-E-WEFF'
    mseq = MM.seq_to_mutation(seq)
    alt_mseq = MM.get_alt_mutation_ssmlist(mseq)
    print(f'{seq=}, {mseq=!s}, alt_mseq={alt_mseq}')
    
    assert seq_fits_mstring('[A1B,D4E,+4-W,G7F]', 'BBC-E-WEFF')            == True
    assert seq_fits_mstring('[A1B,D4E,+4-W,G7F]', 'BBC-E-WEFF',exact=True) == True
    
def test_mstrings():
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'BBC-EW-EFF')            == True
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'BBC-EW-EFF',exact=True) == True
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWEXYEFG')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[A1B,D4E,+4W,G7F]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'BBC-E-WEFF')            == True
    assert seq_fits_mstring('[D4E]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'BBC-EW-EFF')            == True
    assert seq_fits_mstring('[D4E]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'BBC-E--EFF')            == True
    assert seq_fits_mstring('[D4E]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring('[D4E]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[D4E]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'ABC-E--EFG')            == True
    assert seq_fits_mstring('[D4E]', 'ABC-E--EFG',exact=True) == True
    assert seq_fits_mstring('[D4E]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring('[D4E]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[D4E]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[D4E]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWEXYEFG')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWD--EFH')            == True
    assert seq_fits_mstring('[+3W,G7H]', 'ABCWD--EFH',exact=True) == True
    assert seq_fits_mstring('[+3W,G7H]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[+3W,G7H]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-E-WEFF')            == True
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-EW-EFF')            == True
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-E--EFF')            == True
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'BBC-E--EFF',exact=True) == True
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWEXYEFG')            == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[A1B,D4E,G7F]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWEXYEFG',exact=True) == True
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[+3W,D4E,+4XY]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWDXYEFG')            == True
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWDXYEFG',exact=True) == True
    assert seq_fits_mstring('[+3W,+4XY]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWD--EFH')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[+3W,+4XY]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[+3W,+4XY]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[A1A,D4.]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWDXYEFG')            == True
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'ABC-E--EFG')            == True
    assert seq_fits_mstring('[A1A,D4.]', 'ABC-E--EFG',exact=True) == True
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWD--EFH')            == True
    assert seq_fits_mstring('[A1A,D4.]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[A1A,D4.]', 'AB--D--EFG')            == True
    assert seq_fits_mstring('[A1A,D4.]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWEXYEFG')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWDXYEFG')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABC-E--EFG')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWD--EFH')            == True
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'ABCWD--EFH',exact=True) == True
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'AB--D--EFG')            == False
    assert seq_fits_mstring('[+3W,C3C,G7H]', 'AB--D--EFG',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-E-WEFF')            == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-E-WEFF',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-EW-EFF')            == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-EW-EFF',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-E--EFF')            == False
    assert seq_fits_mstring('[A1A,C3.]', 'BBC-E--EFF',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWEXYEFG')            == True
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWEXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWDXYEFG')            == True
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWDXYEFG',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'ABC-E--EFG')            == True
    assert seq_fits_mstring('[A1A,C3.]', 'ABC-E--EFG',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWD--EFH')            == True
    assert seq_fits_mstring('[A1A,C3.]', 'ABCWD--EFH',exact=True) == False
    assert seq_fits_mstring('[A1A,C3.]', 'AB--D--EFG')            == True
    assert seq_fits_mstring('[A1A,C3.]', 'AB--D--EFG',exact=True) == True

if __name__ == '__main__':

    test_hard_mstrings()
    test_seq_to_mutations()
    test_mstrings()
    test_seq_to_mutation()
    
