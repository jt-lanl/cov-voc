'''Table of colors, names, and pango lineages'''

import re
import warnings
from collections import Counter
import itertools

import covid
import colornames

## Roy G Biv
COLORS_DEFAULT = [
    'red', 'yellow', 'blue', 'violet',
    'orange', 'green', 'indigo'
]

## From Will Fischer via Paul Tol
## (this is third row of scheme in Fig.14)
COLORS_DEFAULT = [
    '#771155',
    '#AA4488',
    '#CC99BB',
    '#114477',
    '#4477AA',
    '#77AADD',
    '#117777',
    '#44AAAA',
    '#77CCCC',
    '#777711',
    '#AAAA44',
    '#DDDD77',
    '#774411',
    '#AA7744',
    '#DDAA77',
# comment out last triplet, too similar to first triplet
#    '#771122',
#    '#AA4455',
#    '#DD7788',
]

OTHER='$OTHER' #regexp that doesn't match anything (since it begins with $)

DefaultLineageTable = [

    (      'Black',                                  'Unassigned', r'(None|Unassigned|)'),
    ( 'WhiteSmoke',                                'Recombinants', r'X.*'),
    (   'DarkGray',                                   'Ancestral', r'(A)|(A\.2)|(A\.2\.[^5](\..*)?)|(A.[^2](\..*)?)|(A\.1[01234567](\..*)?)|(A\.[^1][0-9](\..*)?)|(B)|(B\.1\.1\.161)|(B\.1\.14)|(B\.1\.260)|(B\.1[0123589])|(B\.2[036789])|(B\.3(\.1)?)|(B.3[0-9])|(B\.4(\.[124567])?)|(B\.4[012345679])|(B\.5([01235678]?))|(B.6(\.[1234568])?)|(B\.6[01])'),
    (     'Yellow',                                       'D614G', r'$OTHER'),
    (     'Orange',                                       'Alpha', r'(B\.1\.1\.7(\..*)?)|(Q\..*)'),
    (  'LightPink',                                        'Beta', r'B\.1\.351(\..*)?'),
    (  'FireBrick',                                       'Gamma', r'(B\.1\.1\.28)|(P\.1(\..*)?)'),
    (   'DarkCyan',                                          'Mu', r'(B\.1\.621(\..*)?)|(BB.2)'),
    (       'Cyan',                                     'Epsilon', r'B\.1\.42[97](\..*)?'),
    ( 'DodgerBlue',                                        'Iota', r'B\.1\.526'),
    ( 'BlueViolet',                                       'Delta', r'(B\.1\.617\.2)|(AY\..*)'),
    (    'Magenta',                                     'BA.1/BD', r'(B.1.1.529)|(BA\.1(\..*)?)|(BD\..*)'),
    (    'Thistle',                                   'BA.1.1/BC', r'(BA\.1\.1(\..*)?)|(BC\..*)'),
    (  'RoyalBlue',                             'BA.2/B[HJPS]/DD', r'(BA\.2(\..*)?)|(B[HJPS]\..*)|(DD\..*)'),
    (       'Pink',                                'BA.2.12.1/BG', r'(BA\.2\.12\.1(\..*)?)|(BG\..*)'),
    (  'LimeGreen',                                     'BA.4/CS', r'(BA\.4(\..*)?)|(CS\..*)'),
    (  'Chocolate',                                   'BA.4.6/DC', r'(BA\.4\.6(\..*)?)|(DC\..*)'),
    ( 'DodgerBlue', 'BA.5/B[EFKTUVWZ]/C[CDEFGKLNRTU]/D[ABEFGHJQ]', r'(BA\.5(\..*)?)|(B[EFKTUVWZ]\..*)|(C[CDEFGKLNQRTU]\..*)|(D[ABEFGHJQ]\..*)'),
    (  'LightBlue',                                   'BF.7-like', r'(BF\.(7|(11))(\..*)?)|(CP\..*)|(BA.5.2.(6|(13)|(3[45]))(\..*)?)|(BA\.5\.1\.((18)|(27))(\..*)?)|(BA\.5\.10\.1(\..*)?)|(BE\.1\.2(\..*)?)|(DW\..*)'),
    (       'Gold',                                        'BQ.1', r'BQ\.1(\..*)?'),
    (  'Goldenrod',                                'BQ.1.1/C[WZ]', r'(BQ\.1\.1(\..*)?)|(C[WZ]\..*)|(D[KMNP]\..*)'),
    (  'Burlywood',                              'BQ.1.1.22-like', r'(BQ\.1\.1\.((22)|(13)|(20)|(24)|(10)|(21))(\..*)?)|(BQ\.1\.18(\..*)?)'),
    (       'Blue',                                'BA.2.3.20/CM', r'(BA\.2\.3\.20(\..*)?)|(CM\..*)'),
    (   'Cornsilk',                              'BA.5.2.18-like', r'(BA\.5\.2\.18(\..*)?)|(BF\.16(\..*)?)|(C[RQ]\..*)'),
    (    'DarkRed',                                         'XBB', r'XBB(\..*)?'),
    (    'Crimson',                                     'XBB.1.5', r'XBB\.1\.5(\..*)?'),
    ('ForestGreen',               'BA.2.75/B[LMNRY]/C[BHV]/D[SV]', r'(BA\.2\.75(\..*)?)|(B[LMNRY]\..*)|(C[BHV]\..*)|(D[SV]\..*)'),
    (  'LimeGreen',                                        'CH.1', r'CH\.1(\..*)?'),
    (      'Khaki',                                        'BN.1', r'BN\.1(\..*)?'),
    (   'DarkCyan',                                'BA.2.75.2/CA', r'(BA\.2\.75\.2(\..*)?)|(CA\..*)'),
    ('YellowGreen',                           'XBF/BM.1.1.1-like', r'(BM\.1\.1\.1)|(XBF(\..*)?)|(CJ(\..*)?)'),
    (      'Green',                                        'BR.2', r'BR\.2(\..*)?'),

]

##______________________________________________________________________


def rd_lineage_table(filename):
    '''read file with columns: color, name, pattern'''
    lineage_table=[]
    with open(filename) as ftable:
        for line in ftable:
            line = re.sub(r'#.*','',line)
            line = line.strip()
            if not line:
                continue
            try:
                color,name,pattern = line.split()
                lineage_table.append( (color,name,pattern) )
            except ValueError:
                warnings.warn(f'Bad line in lineage file: {line}')
    return lineage_table

def write_lineage_table_python(lineage_table):
    '''write out the lineage table as a python list of tuples'''
    print("DefaultLineageTable = [")
    print()
    maxlen_color = maxlen_name = 0
    for color,name,pattern in lineage_table:
        maxlen_color = max([maxlen_color,len(color)])
        maxlen_name  = max([maxlen_name,len(name)])
    fmt = "    (%%%ds %%%ds %%s)," % (maxlen_color+3, maxlen_name+3)
    for color,name,pattern in lineage_table:
        qcolor = f"'{color}',"
        qname = f"'{name}',"
        qpattern = f"r'{pattern}'"
        print(fmt % (qcolor,qname,qpattern))
    print()
    print("]")

class LineageTablePatterns:
    '''adds some functionality to the raw lineage table'''

    def __init__(self,table):
        self.patterns =  [patt
                          for color,name,patt in table]
        self.colors =    {patt: colornames.tohex(color)
                          for color,name,patt in table}
        self.names =     {patt: name
                          for color,name,patt in table}
        self.regexpatt = {patt: re.compile(r'\.('+patt+r')$')
                          for color,name,patt in table}

    def _match_generator(self,seqname,reverse=False):
        patternlist = self.patterns[::-1] if reverse else self.patterns
        return (patt for patt in patternlist
                if self.regexpatt[patt].search(seqname))

    def last_match(self,seqname,notfound=OTHER):
        '''return the last pattern found whose regexp matches the sequence name'''
        ## by "last" we mean first in the reversed list
        return next(self._match_generator(seqname,reverse=True), notfound)

    def first_match(self,seqname,notfound=OTHER):
        '''return the first pattern found whose regexp matches the sequence name'''
        return next(self._match_generator(seqname),notfound)

    def all_matches(self,seqname):
        '''return a list of all patterns whose regexp matches the sequence name'''
        return list(self._match_generator(seqname))

    def add_pattern(self,color,name,patt):
        '''add a (color,name,pattern) to the table'''
        self.patterns.append(patt)
        self.colors[patt] = colornames.tohex(color)
        self.names[patt] = name
        self.regexpatt[patt] = re.compile(r'\.('+patt+r')$')

    def del_pattern(self,patt):
        '''delete an entry in the table'''
        if patt not in self.patterns:
            warnings.warn(f'pattern {patt} not in lineage table, cannot delete')
            return
        self.patterns.remove(patt)
        del self.colors[patt]
        del self.names[patt]
        del self.regexpatt[patt]

def get_lineage_table(lineagetable_file=None,other=None):
    '''return a lineage table; specifically, a LineageTablePatterns class instance'''

    table = DefaultLineageTable
    if lineagetable_file:
        table = rd_lineage_table(lineagetable_file)

    if other:
        color,name,posn = other
        table.insert(int(posn), (color,name,OTHER) )

    if OTHER not in [items[2] for items in table]:
        table.insert(0,('Gainsboro','other',OTHER))

    return LineageTablePatterns(table)

def get_lineage_table_from_seqs(seqs,num_lineages=15,
                                skipnone=True,other=None):
    '''
    create a lineage table from the most frequent
    pango designations in the list of sequences
    '''
    lineage_counter = Counter()
    for s in seqs:
        lineage = covid.get_lineage_from_name(s.name)
        if skipnone and lineage in ['None','Unassigned','']:
            continue
        lineage_counter[lineage] += 1
    lineage_list = sorted(lineage_counter,
                          key=lineage_counter.get,
                          reverse=True)
    table = []
    colors = itertools.cycle(COLORS_DEFAULT)
    for name in lineage_list[:num_lineages]:
        pattern = re.sub(r'\.',r'\.',name)
        color = next(colors)
        table.append((color,name,pattern))

    ## reverse order of table so most frequent on top
    table = list(reversed(table))

    if other:
        color,name,posn = other
        table.insert(int(posn), (color,name,OTHER) )

    if OTHER not in [items[2] for items in table]:
        table.insert(0,('Gainsboro','other',OTHER))

    return LineageTablePatterns(table)

if __name__ == "__main__":

    ltable = rd_lineage_table('lineage-table-latest.txt')
    write_lineage_table_python(ltable)
