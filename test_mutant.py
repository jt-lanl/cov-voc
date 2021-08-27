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
    'BBC-E-WEFF': '[A1B,D4E,+4W,G7F]!',
}

MM = mutant.MutationManager(RefSeq)
    
def checkit(mstring,seq,**kw):
    mpatt = mutant.Mutation.from_mstring(mstring)
    return MM.seq_fits_pattern(mpatt,seq,**kw)

def test_get_mutations():
    for seq,mstring in mut_dict.items():
        mut = MM.get_mutation(seq)
        assert mut.as_string() == mstring
    for seq,mstring in mut_dict_xtra.items():
        mut = MM.get_mutation(seq)
        assert mut.as_string() == mstring
    for seq,mstring in mut_dict.items():
        mut = MM.get_mutation(seq)
        assert seq == MM.seq_from_mutation(mut)

    for seq,mstring in mut_dict.items():
        ## these are sequences which /cannot/ be reconstructed from mstring
        mut = MM.get_mutation(seq)
        #assert seq != MM.seq_from_mutation(mut)

    

def test_mstrings():
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-E-WEFF')            == True
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-E-WEFF',exact=True) == True
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-EW-EFF')            == True
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-EW-EFF',exact=True) == True
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-E--EFF')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWEXYEFG')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWDXYEFG')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABC-E--EFG')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWD--EFH')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'AB--D--EFG')            == False
    assert checkit('[A1B,D4E,+4W,G7F]', 'AB--D--EFG',exact=True) == False
    assert checkit('[D4E]', 'BBC-E-WEFF')            == True
    assert checkit('[D4E]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[D4E]', 'BBC-EW-EFF')            == True
    assert checkit('[D4E]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[D4E]', 'BBC-E--EFF')            == True
    assert checkit('[D4E]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[D4E]', 'ABCWEXYEFG')            == True
    assert checkit('[D4E]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[D4E]', 'ABCWDXYEFG')            == False
    assert checkit('[D4E]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[D4E]', 'ABC-E--EFG')            == True
    assert checkit('[D4E]', 'ABC-E--EFG',exact=True) == True
    assert checkit('[D4E]', 'ABCWD--EFH')            == False
    assert checkit('[D4E]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[D4E]', 'AB--D--EFG')            == False
    assert checkit('[D4E]', 'AB--D--EFG',exact=True) == False
    assert checkit('[+3W,G7H]', 'BBC-E-WEFF')            == False
    assert checkit('[+3W,G7H]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[+3W,G7H]', 'BBC-EW-EFF')            == False
    assert checkit('[+3W,G7H]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[+3W,G7H]', 'BBC-E--EFF')            == False
    assert checkit('[+3W,G7H]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[+3W,G7H]', 'ABCWEXYEFG')            == False
    assert checkit('[+3W,G7H]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[+3W,G7H]', 'ABCWDXYEFG')            == False
    assert checkit('[+3W,G7H]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[+3W,G7H]', 'ABC-E--EFG')            == False
    assert checkit('[+3W,G7H]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[+3W,G7H]', 'ABCWD--EFH')            == True
    assert checkit('[+3W,G7H]', 'ABCWD--EFH',exact=True) == True
    assert checkit('[+3W,G7H]', 'AB--D--EFG')            == False
    assert checkit('[+3W,G7H]', 'AB--D--EFG',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'BBC-E-WEFF')            == True
    assert checkit('[A1B,D4E,G7F]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'BBC-EW-EFF')            == True
    assert checkit('[A1B,D4E,G7F]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'BBC-E--EFF')            == True
    assert checkit('[A1B,D4E,G7F]', 'BBC-E--EFF',exact=True) == True
    assert checkit('[A1B,D4E,G7F]', 'ABCWEXYEFG')            == False
    assert checkit('[A1B,D4E,G7F]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'ABCWDXYEFG')            == False
    assert checkit('[A1B,D4E,G7F]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'ABC-E--EFG')            == False
    assert checkit('[A1B,D4E,G7F]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'ABCWD--EFH')            == False
    assert checkit('[A1B,D4E,G7F]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[A1B,D4E,G7F]', 'AB--D--EFG')            == False
    assert checkit('[A1B,D4E,G7F]', 'AB--D--EFG',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-E-WEFF')            == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-EW-EFF')            == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-E--EFF')            == False
    assert checkit('[+3W,D4E,+4XY]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'ABCWEXYEFG')            == True
    assert checkit('[+3W,D4E,+4XY]', 'ABCWEXYEFG',exact=True) == True
    assert checkit('[+3W,D4E,+4XY]', 'ABCWDXYEFG')            == False
    assert checkit('[+3W,D4E,+4XY]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'ABC-E--EFG')            == False
    assert checkit('[+3W,D4E,+4XY]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'ABCWD--EFH')            == False
    assert checkit('[+3W,D4E,+4XY]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[+3W,D4E,+4XY]', 'AB--D--EFG')            == False
    assert checkit('[+3W,D4E,+4XY]', 'AB--D--EFG',exact=True) == False
    assert checkit('[+3W,+4XY]', 'BBC-E-WEFF')            == False
    assert checkit('[+3W,+4XY]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[+3W,+4XY]', 'BBC-EW-EFF')            == False
    assert checkit('[+3W,+4XY]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[+3W,+4XY]', 'BBC-E--EFF')            == False
    assert checkit('[+3W,+4XY]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[+3W,+4XY]', 'ABCWEXYEFG')            == True
    assert checkit('[+3W,+4XY]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[+3W,+4XY]', 'ABCWDXYEFG')            == True
    assert checkit('[+3W,+4XY]', 'ABCWDXYEFG',exact=True) == True
    assert checkit('[+3W,+4XY]', 'ABC-E--EFG')            == False
    assert checkit('[+3W,+4XY]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[+3W,+4XY]', 'ABCWD--EFH')            == False
    assert checkit('[+3W,+4XY]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[+3W,+4XY]', 'AB--D--EFG')            == False
    assert checkit('[+3W,+4XY]', 'AB--D--EFG',exact=True) == False
    assert checkit('[D4*]', 'BBC-E-WEFF')            == True
    assert checkit('[D4*]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[D4*]', 'BBC-EW-EFF')            == True
    assert checkit('[D4*]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[D4*]', 'BBC-E--EFF')            == True
    assert checkit('[D4*]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[D4*]', 'ABCWEXYEFG')            == True
    assert checkit('[D4*]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[D4*]', 'ABCWDXYEFG')            == False
    assert checkit('[D4*]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[D4*]', 'ABC-E--EFG')            == True
    assert checkit('[D4*]', 'ABC-E--EFG',exact=True) == True
    assert checkit('[D4*]', 'ABCWD--EFH')            == False
    assert checkit('[D4*]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[D4*]', 'AB--D--EFG')            == False
    assert checkit('[D4*]', 'AB--D--EFG',exact=True) == False
    assert checkit('[A1A,D4.]', 'BBC-E-WEFF')            == False
    assert checkit('[A1A,D4.]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[A1A,D4.]', 'BBC-EW-EFF')            == False
    assert checkit('[A1A,D4.]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[A1A,D4.]', 'BBC-E--EFF')            == False
    assert checkit('[A1A,D4.]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[A1A,D4.]', 'ABCWEXYEFG')            == True
    assert checkit('[A1A,D4.]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[A1A,D4.]', 'ABCWDXYEFG')            == True
    assert checkit('[A1A,D4.]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[A1A,D4.]', 'ABC-E--EFG')            == True
    assert checkit('[A1A,D4.]', 'ABC-E--EFG',exact=True) == True
    assert checkit('[A1A,D4.]', 'ABCWD--EFH')            == True
    assert checkit('[A1A,D4.]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[A1A,D4.]', 'AB--D--EFG')            == True
    assert checkit('[A1A,D4.]', 'AB--D--EFG',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-E-WEFF')            == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-EW-EFF')            == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-E--EFF')            == False
    assert checkit('[+3W,C3C,G7H]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'ABCWEXYEFG')            == False
    assert checkit('[+3W,C3C,G7H]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'ABCWDXYEFG')            == False
    assert checkit('[+3W,C3C,G7H]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'ABC-E--EFG')            == False
    assert checkit('[+3W,C3C,G7H]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[+3W,C3C,G7H]', 'ABCWD--EFH')            == True
    assert checkit('[+3W,C3C,G7H]', 'ABCWD--EFH',exact=True) == True
    assert checkit('[+3W,C3C,G7H]', 'AB--D--EFG')            == False
    assert checkit('[+3W,C3C,G7H]', 'AB--D--EFG',exact=True) == False
    assert checkit('[A1A,C3.]', 'BBC-E-WEFF')            == False
    assert checkit('[A1A,C3.]', 'BBC-E-WEFF',exact=True) == False
    assert checkit('[A1A,C3.]', 'BBC-EW-EFF')            == False
    assert checkit('[A1A,C3.]', 'BBC-EW-EFF',exact=True) == False
    assert checkit('[A1A,C3.]', 'BBC-E--EFF')            == False
    assert checkit('[A1A,C3.]', 'BBC-E--EFF',exact=True) == False
    assert checkit('[A1A,C3.]', 'ABCWEXYEFG')            == True
    assert checkit('[A1A,C3.]', 'ABCWEXYEFG',exact=True) == False
    assert checkit('[A1A,C3.]', 'ABCWDXYEFG')            == True
    assert checkit('[A1A,C3.]', 'ABCWDXYEFG',exact=True) == False
    assert checkit('[A1A,C3.]', 'ABC-E--EFG')            == True
    assert checkit('[A1A,C3.]', 'ABC-E--EFG',exact=True) == False
    assert checkit('[A1A,C3.]', 'ABCWD--EFH')            == True
    assert checkit('[A1A,C3.]', 'ABCWD--EFH',exact=True) == False
    assert checkit('[A1A,C3.]', 'AB--D--EFG')            == True
    assert checkit('[A1A,C3.]', 'AB--D--EFG',exact=True) == True

if __name__ == '__main__':

    test_get_mutations()
    
