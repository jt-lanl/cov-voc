'''read the lineage_notes.txt file'''
## Available at:
## https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
import re
from typing import List
import verbose as v

def first_last_splits(name,reverse=False):
    '''
    For dot-deliminated name of form A.B.C...D,
    iteratively yield tuples that split into first and last parts;
    eg (A, B.C...D), (A.B, C...D), (A.B.C, ...D), etc
    '''
    tokens = name.split('.')
    nlist = range(len(tokens))
    if reverse:
        nlist = reversed(list(nlist))
    for ndx in range(len(tokens)):
        first,last = (".".join(tokens[:ndx]),
                      ".".join(tokens[ndx:]))
        yield first,last


class LineageNotes:
    '''parse the lineage_notes.txt file
    Usage: ln = LineageNotes.from_file(filename)
    Attributes:
       ln.lineages: list of (short) lineage names
       ln.fullname: dict for translating short names to long/full names
       ln.shortname: dict for translating long/full names to short names
    '''

    default_file = "lineage_notes.txt"

    def __init__(self,alias_of,describe=None,fix=False):
        self.fullname = alias_of
        self.shortname = {val:key for key,val in self.fullname.items()}
        self.describe = describe
        if fix:
            self.fix_inconsistencies()
            for bad in self.remove_inconsistencies():
                v.print(bad)

    @property
    def lineages(self):
        '''return all the lineages'''
        return set(self.fullname)

    @classmethod
    def from_file(cls,lineagefile,keepwithdrawn=False,fix=False):
        '''initialize LineageNotes object from file'''
        alias_of,describe = cls.read_lineage_file(lineagefile,
                                        keepwithdrawn=keepwithdrawn)
        return cls(alias_of,describe=describe,fix=fix)

    @staticmethod
    def read_lineage_file(lineagefile,keepwithdrawn=False):
        '''return a list of lineages, and dict of aliases'''
        ## parse lines of the form:
        ## GL.1	Alias of XAY.1.1.1.1, Europe, S:D420N, C19441T, from issue #2032
        re_alias = re.compile(r'[Aa]lias\s+of\s+(((B)|(X[A-Z]+))[\.0-9]+)')
        lineages = list()
        alias_of = dict()
        describe = dict()
        with open(lineagefile,'r') as fin:
            fin.readline() # skip header: "Lineage Description"
            for line in fin:
                if not line.strip():
                    continue
                name,description = line.split(None,1)
                if not keepwithdrawn and name.startswith('*'):
                    ## asterisk indicates withdrawn lineage
                    continue
                lineages.append(name)
                describe[name] = description.strip()
                if (mat := re_alias.search(description)):
                    v.vvvprint(f'Alias: {name} -> {mat[1]}')
                    alias_of[name] = mat[1]
                elif re.match(r'[A-Z][A-Z]\..*',name):
                    ## all two-letter names should have aliases
                    v.print('missing "Alias of"?',line)
                else:
                    alias_of[name] = name

        assert len(lineages) == len(alias_of)
        v.vprint(f'Lineages:  {len(alias_of)}')
        return alias_of,describe

    def report_size(self):
        '''return string with size of lineages and aliases'''
        sizelines = [f'lineages: {len(self.lineages)}',
                     f'aliases : {len(self.fullname)}']
        return "\n".join(sizelines)

    def get_fullname(self,short):
        '''more robust version of fullname.get(),
        if short is not itn ehte fullname dict,
        then infers fullname
        by trimming the last numbers of the end and trying again
        '''
        long = self.fullname.get(short,None)
        if not long:
            for first,last in first_last_splits(short):
                firstlong = self.fullname.get(first,None)
                if firstlong:
                    long = f'{firstlong}.{last}'
                    break
        return long or short

    def get_shortname(self,long):
        '''more robust version of shortname.get()'''
        long = self.get_fullname(long) or long
        short = self.shortname.get(long,None)
        if not short:
            for first,last in first_last_splits(long,reverse=True):
                shortfirst = self.shortname.get(first,None)
                if shortfirst:
                    short = f'{shortfirst}.{last}'
        return short

    def inconsistencies(self,remove=True,fix=False):
        '''check consistency of lineage notes file'''
        inconsistent=[]
        fixed=[]
        bad_aliases = list()
        fixed_aliases = dict()
        for alias,full in self.fullname.items():
            try:
                _,numbers = alias.split('.',1)
            except ValueError:
                ## if alias has no numbers (eg, is B or XBB), continue
                continue
            if not full.endswith(numbers):
                bad_aliases.append(alias)
                inconsistent.append(f'Inconsistent ending: {alias}->{full}')
                if not fix:
                    continue
                parent_alias = re.sub(r'\.\d+$','',alias)
                lastno_alias = re.sub(r'.*(\.\d+)$',r'\1',alias)
                if parent_alias in self.fullname:
                    parent_full = self.fullname[parent_alias]
                    new_full = parent_full + lastno_alias
                else:
                    # just fix the last number
                    new_full = re.sub(r'(.*)\.\d+$',r'\1'+lastno_alias,full)
                fixed_aliases[alias] = new_full
                fixed.append(f'Fix {alias}: {full} -> {new_full}')
        if fix:
            for alias,new_full in fixed_aliases.items():
                self.fullname[alias] = new_full
        if remove:
            ## remove inconsistent aliases/lineages
            for alias in bad_aliases:
                del self.fullname[alias]
        if fix:
            return fixed
        return inconsistent

    def remove_inconsistencies(self):
        '''remove inconsistent lineage definitions'''
        return self.inconsistencies(remove=True,fix=False)
    def fix_inconsistencies(self):
        '''fix inconsistent lineage definitions'''
        return self.inconsistencies(remove=False,fix=True)

    def parent_of(self,child):
        '''return the lineage name that is the parent of the input lineage'''
        full = self.fullname.get(child,child)
        if re.search(r"\.\d+$",full):
            fullparent = re.sub(r"\.\d+$","",full)
        else:
            fullparent = "Wuhan"
        parent = self.shortname.get(fullparent,fullparent)
        return parent

    def children_of(self, parent: str) -> List[str]:
        """return list of lineages that are immediate chlindren of parent"""
        children = []
        parent = self.fullname.get(parent, parent)
        pcount = parent.count(".")
        for lin in self.lineages:
            longlin = self.fullname.get(lin, lin)
            if longlin == parent:
                continue
            if longlin.startswith(parent) and longlin.count(".") == pcount + 1:
                children.append(lin)
        return children

    def allchildren_of(self,parent,excludeparent=False):
        '''return a list of lineages, all of which are children
        (or grandchildren, or great grandchildren, etc) of the parent
        NOTE: The parent is included, unless excludeparent=True
        '''
        children=[]
        parent = self.fullname.get(parent,parent)
        for lin in self.lineages:
            longlin = self.fullname.get(lin,lin)
            if excludeparent and longlin == parent:
                continue
            if longlin.startswith(parent):
                children.append(lin)
        return children

    def get_lineage_set(self,parent=None,excludeparent=False):
        '''return a set of lineages that have parent as an ancestor'''
        lineage_set = self.lineages
        if parent and parent in lineage_set:
            lineage_set = self.allchildren_of(parent,
                                              excludeparent=excludeparent)
        return set(lineage_set)

if __name__ == "__main__":
    import covid
    file = covid.find_seqfile(LineageNotes.default_file)
    lins = LineageNotes.from_file(file)
    print("Inconsistencies")
    for inc in lins.inconsistencies(remove=False):
        print(inc)
    print("Fixes")
    for eachfix in lins.fix_inconsistencies():
        print(eachfix)

    clade = "BA.2"
    while clade:
        mychildren = lins.allchildren_of(clade,
                                         excludeparent=True)
        print(f'clade={clade} says: I have {len(mychildren)} children!')
        clade = mychildren[0] if mychildren else None
