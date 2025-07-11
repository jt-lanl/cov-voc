"""read the lineage_notes.txt file"""

## Available at:
## https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
import re
from collections.abc import Iterator
import json
import warnings
import verbose as v

RE_INT = re.compile(r"^[0-9]+$")

def letters_numbers_split(name: str):
    """
    splta name like AB.7.1.5 into AB and 7.1.5
    assume letters all appear before first "."
    then all numbers after that
    """
    try:
        letters,numbers = name.split(".",1)
    except ValueError:
        letters = name
        numbers = None
    return letters,numbers
    
def first_last_splits(name: str, reverse: bool = False):
    """
    For dot-delimited name of form A.B.C...D,
    iteratively yield tuples that split into first and last parts;
    eg (A, B.C...D), (A.B, C...D), (A.B.C, ...D), etc
    """
    tokens = name.split(".")
    nlist = range(len(tokens))
    if reverse:
        nlist = reversed(list(nlist))
    for ndx in range(len(tokens)):
        first, last = (".".join(tokens[:ndx]), ".".join(tokens[ndx:]))
        yield first, last


class LineageNotes:
    """parse the lineage_notes.txt file
    Usage: ln = LineageNotes.from_file(notesfile,keyfile,fix=True)
    Attributes:
       ln.lineages: list of (short) lineage names
       ln.fullname: dict for translating short names to long/full names
       ln.shortname: dict for translating long/full names to short names
    """

    default_file = "lineage_notes.txt"
    default_key = "alias_key.json"

    def __init__(self, alias_of: dict[str, str],
                 describe: str = None,
                 fix: bool = False,
                 alias_key: dict[str,str] = None):
        self.fullname = alias_of
        self.shortname = {val: key for key, val in self.fullname.items()}
        self.describe = describe
        self.alias_key = alias_key
        if fix:
            self.fix_inconsistencies()
            for bad in self.remove_inconsistencies():
                v.print(bad)

    @property
    def lineages(self):
        """return all the lineages"""
        return set(self.fullname)

    @classmethod
    def from_file(cls, lineagefile: str,
                  aliaskeyfile: str, # = None,
                  keepwithdrawn: bool = False,
                  fix: bool = False):
        """initialize LineageNotes object from file"""
        if aliaskeyfile:
            alias_key = cls.read_alias_key_file(aliaskeyfile)
            v.vprint(f"alias keyfile = {aliaskeyfile}")
        else:
            alias_key=None
            warnings.warn(f"File alias_key.json NOT used")
        alias_of, describe = cls.read_lineage_file(lineagefile, keepwithdrawn=keepwithdrawn)
        return cls(alias_of, describe=describe, fix=fix, alias_key=alias_key)

    @staticmethod
    def read_alias_key_file(alias_key_file):
        """read alias_key.json file to get alias_key dict"""
        with open(alias_key_file) as fp:
            alias_key = json.load(fp)
        ## exclude dict elements where full is empty or short begins with "X..."
        alias_key = {short:full for short,full in alias_key.items()
                     if full and not short.startswith("X")}
        ## two special cases
        alias_key["A"]="A"
        alias_key["B"]="B"
        return alias_key

    @staticmethod
    def read_lineage_file(lineagefile: str, keepwithdrawn: bool = False):
        """return a list of lineages, and dict of aliases"""
        ## parse lines of the form:
        ## GL.1	Alias of XAY.1.1.1.1, Europe, S:D420N, C19441T, from issue #2032
        re_alias = re.compile(r"[Aa]lias\s+of\s+(((B)|(X[A-Z]+))[\.0-9]+)")
        lineages = list()
        alias_of = dict()
        describe = dict()
        with open(lineagefile, "r") as fin:
            fin.readline()  # skip header: "Lineage Description"
            for line in fin:
                if not line.strip():
                    continue
                name, description = line.split(None, 1)
                if not keepwithdrawn and name.startswith("*"):
                    ## asterisk indicates withdrawn lineage
                    continue
                if name in lineages:
                    v.print(f"Name {name} already appeared in lineage notes tile")
                else:
                    lineages.append(name)
                alias_of[name]=name ## default?
                describe[name] = description.strip()
                if mat := re_alias.search(description):
                    v.vvvprint(f"Alias: {name} -> {mat[1]}")
                    alias_of[name] = mat[1]
                elif re.match(r"[A-Z][A-Z]\..*", name):
                    ## all two-letter names should have aliases
                    v.print('missing "Alias of"?', line)
                else:
                    alias_of[name] = name

        v.vprint(f"Lineages:  {len(alias_of)} =?= {len(lineages)}")
        assert len(lineages) == len(alias_of)
        return alias_of, describe

    def report_size(self):
        """return string with size of lineages and aliases"""
        sizelines = [f"lineages: {len(self.lineages)}", f"aliases : {len(self.fullname)}"]
        return "\n".join(sizelines)

    def get_fullname(self, short: str):
        """more robust version of fullname.get(),
        if short is not itn ehte fullname dict,
        then infers fullname
        by trimming the last numbers of the end and trying again
        """
        long = self.fullname.get(short, None)
        if not long:
            for first, last in first_last_splits(short):
                firstlong = self.fullname.get(first, None)
                if firstlong:
                    long = f"{firstlong}.{last}"
                    break
        return long or short

    def get_shortname(self, long: str):
        """more robust version of shortname.get()"""
        long = self.get_fullname(long) or long
        short = self.shortname.get(long, None)
        if not short:
            for first, last in first_last_splits(long, reverse=True):
                shortfirst = self.shortname.get(first, None)
                if shortfirst:
                    short = f"{shortfirst}.{last}"
        return short

    def sortkey(self, lin: str) -> tuple:
        """key used for sorted() to sort lineage names based on fullnamnes"""
        return tuple(int(t) if RE_INT.match(t) else t for t in self.get_fullname(lin).split("."))

    def sort_lineages(self, lineagenames: Iterator[str]) -> list[str]:
        """return list of lineage names, sorted"""
        return sorted(lineagenames, key=self.sortkey)

    def alias_key_expand(self,alias):
        if not self.alias_key:
            return None
        try:
            letters,numbers = alias.split(".",1)
        except ValueError: ### Shouldn't happen if called after split in inconsistencies fcn
            letters = alias
            numbers = None
        if letters in self.alias_key:
            akfull = self.alias_key[letters]
            if numbers:
                akfull = akfull + "." + numbers
            return akfull
        return None # no evidence of inconsistency
        

    def inconsistencies(self, remove: bool = True, fix: bool = False) -> list[str]:
        """check consistency of lineage notes file"""
        inconsistent = []
        fixed = []
        bad_aliases = list()
        fixed_aliases = dict()
        for alias, full in self.fullname.items():

            ## First check if alias,full are consistent based on alias_key
            ## if alias_key is not available, then this part assumes consistency
            akfull = self.alias_key_expand(alias) or full
            if akfull != full:
                bad_aliases.append(alias)
                inconsistent.append(f"Invalid alias: {alias}: {full} -> {akfull}")
                if fix:
                    fixed_aliases[alias] = akfull
                    fixed.append(f"Fix {alias}: {full} -> {akfull}")
                full = akfull

            ## Second check, if tailing numbers are consistent
            letters,numbers = letters_numbers_split(alias)
            if not numbers:
                continue
            if not full.endswith(numbers):
                bad_aliases.append(alias)
                inconsistent.append(f"Inconsistent ending: {alias}->{full}")
                if not fix:
                    continue
                parent_alias = re.sub(r"\.\d+$", "", alias)
                lastno_alias = re.sub(r".*(\.\d+)$", r"\1", alias)
                if parent_alias in self.fullname:
                    parent_full = self.fullname[parent_alias]
                    new_full = parent_full + lastno_alias
                else:
                    # just fix the last number
                    new_full = re.sub(r"(.*)\.\d+$", r"\1" + lastno_alias, full)
                fixed_aliases[alias] = new_full
                fixed.append(f"Fix {alias}: {full} -> {new_full}")
        if fix:
            for alias, new_full in fixed_aliases.items():
                self.fullname[alias] = new_full
                self.shortname[new_full] = alias
        if remove:
            ## remove inconsistent aliases/lineages
            for alias in bad_aliases:
                del self.fullname[alias]
        if fix:
            return fixed
        return inconsistent

    def remove_inconsistencies(self) -> list[str]:
        """remove inconsistent lineage definitions"""
        return self.inconsistencies(remove=True, fix=False)

    def fix_inconsistencies(self) -> list[str]:
        """fix inconsistent lineage definitions"""
        return self.inconsistencies(remove=False, fix=True)

    def parent_of(self, child: str) -> str:
        """return the lineage name that is the parent of the input lineage"""
        full = self.fullname.get(child, child)
        if re.search(r"\.\d+$", full):
            fullparent = re.sub(r"\.\d+$", "", full)
        else:
            fullparent = "Wuhan"
        parent = self.shortname.get(fullparent, fullparent)
        return parent

    def children_of(self, parent: str) -> list[str]:
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

    def allchildren_of(self, parent: str, excludeparent: bool = False) -> list[str]:
        """return a list of lineages, all of which are children
        (or grandchildren, or great grandchildren, etc) of the parent
        NOTE: The parent is included, unless excludeparent=True
        """
        children = []
        longparent = self.fullname.get(parent, parent)
        for child in self.lineages:
            longchild = self.fullname.get(child, child)
            if longchild.startswith(longparent + ".") or (
                not excludeparent and longchild == longparent
            ):
                children.append(child)
        return children

    def get_lineage_set(self, parent: str = None, excludeparent: str = False):
        """return a set of lineages that have parent as an ancestor"""
        lineage_set = self.lineages
        if parent and parent in lineage_set:
            lineage_set = self.allchildren_of(parent, excludeparent=excludeparent)
        return set(lineage_set)


if __name__ == "__main__":
    import covid

    file = covid.find_seqfile(LineageNotes.default_file)
    akfile = covid.find_seqfile(LineageNotes.default_key)
    v.print("akfile=",akfile)
    lins = LineageNotes.from_file(file,akfile)
    print("Inconsistencies")
    for inc in lins.inconsistencies(remove=False):
        print(inc)
    print("Fixes")
    for eachfix in lins.fix_inconsistencies():
        print(eachfix)

    clade = "BA.2"
    while clade:
        mychildren = lins.allchildren_of(clade, excludeparent=True)
        print(f"clade={clade} says: I have {len(mychildren)} children!")
        clade = mychildren[0] if mychildren else None
