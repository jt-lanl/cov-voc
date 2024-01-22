'''read the lineage_notes.txt file'''
## Available at:
## https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
import re
from functools import cache
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
    for n in range(len(tokens)):
        first,last = (".".join(tokens[:n]),
                      ".".join(tokens[n:]))
        yield first,last


def read_lineage_file(lineagefile,keepwithdrawn=False):
    '''return a list of lineages, and dict of aliases'''
    ## parse lines of the form:
    ## GL.1	Alias of XAY.1.1.1.1, Europe, S:D420N, C19441T, from issue #2032
    lineages = list()
    fullname = dict()
    with open(lineagefile,'r') as fin:
        ## skip top two lines
        fin.readline()
        fin.readline()
        ## begin reading in earnest
        for line in fin:
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            name = tokens[0]
            if name[0]=='*' and not keepwithdrawn:
                ## asterisk indicates withdrawn lineage
                continue
            lineages.append(name)
            if tokens[1] == 'Alias':
                full = tokens[3]
                if full[-1] == ',':
                    ## strip trailing comma
                    full = full[:-1]
                fullname[name] = full
            elif re.match(r'[A-Z][A-Z]\..*',tokens[0]):
                v.vprint('missing "Alias of"?',line)
    return lineages,fullname

class LineageNotes:
    '''parse the lineage_notes.txt file'''
    ## lineages: list of (short) lineage names
    ## fullname: dict for translating short names to long/full names
    ## shortname: dict for translating long/full names to short names
    def __init__(self,lineagefile=None):
        if lineagefile:
            self.lineages, self.fullname = read_lineage_file(lineagefile)
            self.shortname = {val:key for key,val in self.fullname.items()}
        else:
            self.lineages=[]
            self.fullname=dict()
            self.shortname=dict()

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
            _,numbers = alias.split('.',1)
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
                    #just fix the last number
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
            self.lineages = [lin for lin in self.lineages
                             if lin not in bad_aliases]
        if fix:
            return fixed
        else:
            return inconsistent

    def remove_inconsistencies(self):
        return self.inconsistencies(remove=True,fix=False)
    def fix_inconsistencies(self):
        return self.inconsistencies(remove=False,fix=True)
    @cache
    def parent_of(self,child):
        '''return the lineage name that is the parent of the input lineage'''
        full = self.fullname.get(child,child)
        if re.search(r"\.\d+$",full):
            fullparent = re.sub(r"\.\d+$","",full)
        else:
            fullparent = "Wuhan"
        parent = self.shortname.get(fullparent,fullparent)
        return parent

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

    def restrict_to_clade(self,parent,excludeparent=False):
        '''return a LineageNotes object that restricts itself
        lineages that originate with specified parent'''
        children = self.allchildren_of(parent,excludeparent=excludeparent)
        new_lin_notes = LineageNotes()
        for lin in children:
            new_lin_notes.lineages.append(lin)
            new_lin_notes.shortname[lin] = self.shortname.get(lin,lin)
            new_lin_notes.fullname[lin] = self.fullname.get(lin,lin)
        return new_lin_notes        

if __name__ == "__main__":
    lins = LineageNotes("data/lineage_notes.txt")
    print("Inconsistencies")
    for inc in lins.inconsistencies(remove=False):
        print(inc)
    print("Fixes")
    for eachfix in lins.fix_inconsistencies():
        print(eachfix)
        
    if 0:
        parent = "BA.2"
        print(f'{parent=}')
        print("Children:\n",
              ",".join(lins.allchildren_of(parent,
                                           excludeparent=True)))

        newlins = lins.restrict_to_clade('BA.2')
        print("New lineage notes:\n",",".join(newlins.lineages))
    
