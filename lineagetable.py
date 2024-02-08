'''Table of colors, names, and pango lineages'''

import re
import warnings
from collections import Counter
import functools
import itertools

import covid
import colornames
import defaultlineage

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
    fmt = "    (%%%ds %%%ds %%s)," % (maxlen_color+3,
                                      maxlen_name+3)
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

        #debug
        #for color,name,patt in table:
        #    print(color,name,patt)
        #    re.compile(r'\.('+patt+r')$')

        self.regexpatt = {patt: re.compile(r'\.('+patt+r')$')
                          for color,name,patt in table}

    def _match_generator(self,seqname,reverse=False):
        patternlist = self.patterns[::-1] if reverse else self.patterns
        return (patt for patt in patternlist
                if self.regexpatt[patt].search(seqname))

    @functools.lru_cache(maxsize=None)
    def last_match(self,seqname,notfound=OTHER):
        '''return the last pattern found 
        whose regexp matches the sequence name'''
        ## by "last" we mean first in the reversed list
        return next(self._match_generator(seqname,reverse=True),
                    notfound)

    def first_match(self,seqname,notfound=OTHER):
        '''return the first pattern found 
        whose regexp matches the sequence name'''
        return next(self._match_generator(seqname),notfound)

    def all_matches(self,seqname):
        '''return a list of all patterns 
        whose regexp matches the sequence name'''
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
            warnings.warn(f'pattern {patt} not in '
                          'lineage table, cannot delete')
            return
        self.patterns.remove(patt)
        del self.colors[patt]
        del self.names[patt]
        del self.regexpatt[patt]

def get_lineage_table(lineagetable_file=None,other=None):
    '''return a lineage table; specifically, 
    a LineageTablePatterns class instance'''

    table = defaultlineage.DefaultLineageTable
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
        lineage = covid.get_lineage(s.name)
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
