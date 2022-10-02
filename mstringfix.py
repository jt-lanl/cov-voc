'''module for "fixing" mstrings to account for equivalent alignments'''

import re
import warnings
import verbose as v

MSTRING_DEFAULT_FIXES = '''
[E156-,F157-,R158G] [E156G,F157-,R158-]
[+214AAG,D215Y] [D215A,+215AGY]
[Y144Y,Y145-] [Y144-,Y145Y]
[N211I,L212-] [N211-,L212I]
[L244Y,H245-] [L244-,H245Y]
[L24S,P25-,P26-,A27-] [L24-,P25-,P26-,A27S]
[S247S,Y248-,L249-,T250-,P251-,G252-,D253-,S254-] [S247-,Y248-,L249-,T250-,P251-,G252-,D253-,S254S]
[S247S,Y248-,L249-,T250-,P251-,G252-,D253-,S256-] [S247-,Y248-,L249-,T250-,P251-,G252-,D253-,S256S]
[Y248S,L249-,T250-,P251-,G252-,D253-,S254-] [Y248-,L249-,T250-,P251-,G252-,D253-]
[V213S,R214G,+214RGR] [+212SGR,V213G,R214R]
[Y144T,Y145S,+145N] [+143T,Y144S,Y145N]
[V213S,+213G,R214R,+214GR] [+212SGR,V213G,R214R]
[D215G,+215ARN] [+214GAR,D215N]
[L242Y,A243-,L244-,H245-] [L242-,A243-,L244-,H245Y]
[+214T,D215D,+215RD] [+214TDR,D215D]
[+209V,N211-,L212-] [I210V,N211-,L212I]
[+18R,T19T,T20T,R21-] [T19R,T20T,R21T]
'''

def mstring_nobrackets(mstring):
    '''clean up mstring, including removal of brackets'''
    mstring = re.sub(r'\s','',mstring) ## remove extra spaces between ssms
    mstring = re.sub(r'\s*ancestral\s*','',mstring) ## take out 'ancestral'
    mstring = re.sub(r'[\[\]]','',mstring) ## remove brackets if there
    return mstring

def mstring_brackets(mstring):
    '''clean up mstring and make sure it has brackets around it'''
    mstring = mstring_nobrackets(mstring) ## first clean-up and remove
    mstring = f'[{mstring}]' ## then add back
    return mstring

def get_mstring_pairs(lines):
    '''parse lines to get mstring pairs'''
    for line in lines:
        line = line.strip()
        line = re.sub('#.*','',line)
        if not line:
            continue

        match = re.search(r'\[(.*)\].+\[(.*)\]',line)
        ## although specified with bracket in string or file
        ## the mstrings themselves do not contain the brackets
        if not match:
            warnings.warn("Could not read line:",line)
            continue
        yield match[1],match[2]


def default_mstring_pairs():
    '''get good/bad mstring pairs from default string'''
    yield from get_mstring_pairs( MSTRING_DEFAULT_FIXES.splitlines() )

def read_mstring_pairs(filename):
    '''read good/bad mstring pairs from a file'''
    if not filename:
        yield from default_mstring_pairs()
        return
    with open(filename) as f_in:
        yield from get_mstring_pairs(f_in.readlines())

class _MStringPair:
    '''pair of mstrings (from->to), along with re.compile'd version of from_mstring'''
    def __init__(self,from_mstring,to_mstring):
        self.a = from_mstring
        self.b = to_mstring
        a_re = self.a
        a_re = re.sub(r'\+',r'\\+',a_re)
        a_re = re.sub(r'\-',r'\\-',a_re)
        self.a_re = re.compile(a_re)

    def __str__(self):
        return f'[{self.a}]->[{self.b}]'

class MStringFixer:
    '''provides the function fix() that fixes mstrings'''

    def __init__(self,tweakfile=None):
        self.mstring_pairs = [_MStringPair(a,b) for a,b in read_mstring_pairs(tweakfile)]

    def append(self,mstring_from,mstring_to):
        '''append pair mstring_from -> mstring_to to list of fixes'''
        self.mstring_pairs.append( _MStringPair(mstring_from,mstring_to) )

    def __str__(self):
        return '\n'.join(str(fix) for fix in self.mstring_pairs)

    def fix(self,mstring):
        '''fix an mstring by substituting standard alignment variants'''
        for pair in self.mstring_pairs:
            if pair.a_re.search(mstring):
                v.vprint(f'mstring={mstring}\n'
                         f'Replace: {pair.a} -> {pair.b}')
                mstring = pair.a_re.sub(pair.b,mstring)
        return mstring
