'''Table of colors, names, and pango lineages'''

import re
import warnings
import colornames

OTHER='$OTHER' #regexp that doesn't match anything (since it begins with $)

DefaultLineageTable=[

    ('Orange',     'Alpha',   r'(B\.1\.1\.7)|(Q\..*)'),
    ('ForestGreen','Lambda',  r'C\.37(\..*)?'),
    ('LightPink',  'Beta',    r'B\.1\.351(\..*)?'),
    ('LimeGreen',  'Mu',      r'B\.1\.621(\..*)?'),
    ('FireBrick',  'Gamma',   r'P\.1(\..*)?'),
    ('Cyan',       'Epsilon', r'B\.1\.42[97](\..*)?'),
    ('DodgerBlue', 'Iota',    r'B\.1\.526(\..*)?'),
    ('Gold',       'Eta',     r'B\.1\.525(\..*)?'),
    ('OliveDrab',  'Kappa',   r'B\.1\.617\.1(\..*)?'),
    ('Goldenrod',  'R.1',     r'R\.1(\..*)?'),

    #('Lavender',  'Delta_AY.25',    r'AY\.25(\..*)?'),
    #('Purple',    'Delta_AY.26',    r'AY\.26(\..*)?'),
    #('Lavender',  'Delta_AY.33',    r'AY\.33(\..*)?'),
    #('SkyBlue',   'Delta_AY.98.1',  r'AY\.98\.1(\..*)?'),
    #('Cornsilk',  'Delta_AY.47',    r'AY\.47(\..*)?'),
    #('HotPink',   'Delta_AY.35',    r'AY\.35(\..*)?'),
    #('Magenta',   'Delta_AY.4.2',   r'AY\.4\.2'),
    #('Tan',       'Delta_AY.4.2.1', r'AY\.4\.2\.1(\..*)?'),
    ('BlueViolet', 'Delta',          r'(B\.1\.617\.2)|(AY\..*)'),

    #('Yellow',       'C.1.2', r'C\.1\.2(\..*)?'),
    #('Khaki',        'B.1.1.318', r'B\.1\.1\.318(\..*)?'),
    #('DarkTurquoise','B.1.640', r'B\.1\.640(\..*)?'),
    ('Red',        'Omicron', r'(B\.1\.1\.529)|(BA\.1)|(BA\.[^12](\..*)?)'),
    ('Pink',       'Omicron_BA.1.1',r'BA\.1\..*'),
    ('Maroon',     'Omicron_BA.2',r'BA\.2(\..*)?'),
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

class LineageTablePatterns:
    '''adds some functionality to the raw lineage table'''

    def __init__(self,table):
        self.patterns =  [patt                                for color,name,patt in table]
        self.colors =    {patt: colornames.tohex(color)       for color,name,patt in table}
        self.names =     {patt: name                          for color,name,patt in table}
        self.regexpatt = {patt: re.compile(r'\.('+patt+r')$') for color,name,patt in table}

    def _match_generator(self,seqname):
        return (patt for patt in self.patterns
                if self.regexpatt[patt].search(seqname))

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
