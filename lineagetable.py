'''Table of colors, names, and pango lineages'''

import re
import warnings
import colornames

OTHER='$OTHER' #regexp that doesn't match anything (since it begins with $)

DefaultLineageTable = [

    (      'Black',                                  'Unassigned', r'(None|Unassigned|)'),
    ( 'WhiteSmoke',                                'Recombinants', r'X.*'),
    (   'DarkGray',                                   'Ancestral', r'(A)|(A.1)|(A.11)|(A.12)|(A.15)|(A.16)|(A.17)|(A.2)|(A.2.2)|(A.2.3)|(A.2.4)|(A.21)|(A.22)|(A.23)|(A.23.1)|(A.24)|(A.25)|(A.26)|(A.27)|(A.28)|(A.29)|(A.3)|(A.30)|(A.4)|(A.5)|(A.6)|(A.7)|(A.9)|(B)|(B.1.1.161)|(B.1.14)|(B.1.260)|(B.10)|(B.11)|(B.12)|(B.13)|(B.15)|(B.18)|(B.19)|(B.20)|(B.23)|(B.26)|(B.27)|(B.28)|(B.29)|(B.3)|(B.3.1)|(B.30)|(B.31)|(B.32)|(B.33)|(B.34)|(B.35)|(B.36)|(B.37)|(B.38)|(B.39)|(B.4)|(B.4.1)|(B.4.2)|(B.4.4)|(B.4.5)|(B.4.6)|(B.4.7)|(B.40)|(B.41)|(B.42)|(B.43)|(B.44)|(B.45)|(B.46)|(B.47)|(B.49)|(B.5)|(B.50)|(B.51)|(B.52)|(B.53)|(B.55)|(B.56)|(B.57)|(B.58)|(B.6)|(B.6.1)|(B.6.2)|(B.6.3)|(B.6.4)|(B.6.5)|(B.6.6)|(B.6.8)|(B.60)|(B.61)'),
    (     'Yellow',                                       'D614G', r'$OTHER'),
    (     'Orange',                                       'Alpha', r'(B\.1\.1\.7(\..*)?)|(Q\..*)'),
    (  'LightPink',                                        'Beta', r'B\.1\.351(\..*)?'),
    (  'FireBrick',                                       'Gamma', r'(B\.1\.1\.28)|(P\.1(\..*)?)'),
    (   'DarkCyan',                                          'Mu', r'(B\.1\.621(\..*)?)|(BB.2)'),
    (       'Cyan',                                     'Epsilon', r'B\.1\.42[97](\..*)?'),
    ( 'DodgerBlue',                                        'Iota', r'B\.1\.526'),
    ( 'BlueViolet',                                       'Delta', r'(B\.1\.617\.2)|(AY\..*)'),
    (    'Magenta',                             'Omicron_BA.1/BD', r'(B.1.1.529)|(BA\.1(\..*)?)|(BD\..*)'),
    (    'Thistle',                           'Omicron_BA.1.1/BC', r'(BA\.1\.1(\..*)?)|(BC\..*)'),
    (  'RoyalBlue',                'Omicron_BA.2/B[HJSP]/CM/CS.1', r'(BA\.2(\..*)?)|(B[SPHJ]\..*)|(CM\.(.*))|(CS\.1(\..*)?)'),
    (  'Burlywood',                        'Omicron_BA.2.12.1/BG', r'(BA\.2\.12\.1(\..*)?)|(BG\..*)'),
    (  'LimeGreen',                          'Omicron_BA.4/CS.2+', r'(BA\.4(\..*)?)|(CS\.[2-9].*)'),
    (  'Chocolate',                              'Omicron_BA.4.6', r'BA\.4\.6(\..*)?'),
    ('YellowGreen', 'Omicron_BA.5/B[EFKTQUVWZ]/C[CDEFGKLNUPQRTU]', r'(BA\.5(\..*)?)|(B[EFKTQUVWZ]\..*)|(C[CDEFGKLNUPQRTU]\..*)'),
    (    'Thistle',             'Omicron_BA.2.75/B[LMNRY]/C[BJH]', r'(BA\.2\.75(\..*)?)|(B[LMNRY]\..*)|(C[BJH]\..*)'),
    (    'Magenta',                        'Omicron_BA.2.75.2/CA', r'(BA\.2\.75\.2(\..*)?)|(CA\..*)|'),
    (  'LightBlue',                                        'BF.7', r'BF\.7(\..*)?'),
    (       'Gold',                                        'BQ.1', r'BQ\.1(\..*)?'),
    (  'Goldenrod',                                      'BQ.1.1', r'BQ\.1\.1(\..*)?'),
    ('ForestGreen',                                         'XBB', r'XBB(\..*)?'),

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
        self.patterns =  [patt                                for color,name,patt in table]
        self.colors =    {patt: colornames.tohex(color)       for color,name,patt in table}
        self.names =     {patt: name                          for color,name,patt in table}
        self.regexpatt = {patt: re.compile(r'\.('+patt+r')$') for color,name,patt in table}

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



if __name__ == "__main__":

    ltable = rd_lineage_table('lineage-table-latest.txt')
    write_lineage_table_python(ltable)
