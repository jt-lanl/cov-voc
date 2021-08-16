'''
mutantx: this version does insertions

Mutation and SingleSiteMutation classes,
where Mutation is a list of SingleSiteMutations, with auxiliary information such as an 'exact' flag.
A SingleSiteMutation has three components:
  site = integer site at which mutation occurs
  ref  = character amino acid at that site in reference sequence
         + indicates insertion after that site
  mut  = character amino acid at that site in mutated sequence
         . indicates any character
         _ indicates ref
         * indicates any character except ref
         multicharacter strings indicate either:
           if ref == +, then multiple characters are inserted
           otherwise, indicates any character in the multicharacter string

Note that there are two kinds of mutantion strings: prescriptive and descriptive
Prescriptive: allows us create a mutation sequence from a reference sequence
Descriptive: we can assess whether a sequence is described by the mutation

eg, [A222V,D614G] is prescriptive (and being prescriptive, is also descriptive)
and [A222*,D614G] is not prescriptive, but is descriptive
    [+222FF] indicates that "FF" is inserted after site 222; this is prescriptive
    [A222DG] indicates that A mutates to either D or G

'''

import sys
import re
import itertools as it
from collections import defaultdict
from functools import lru_cache
import warnings

class SiteIndexTranslator():
    '''provide functions to translate between string index and site number
       refseq: ABC--D
       site:   123334
       index:  012345
    '''
    def __init__(self,refseq):
        # self.refseq = refseq ## needed?
        self.site = list(it.accumulate(int(b!='-') for b in refseq))
        self.ndx = [None]*(2+max(self.site))
        for n,p in enumerate(self.site):
            if self.ndx[p] is None:
                self.ndx[p] = n
        ## final index is len(refseq)=len(self.site)
        self.ndx[-1] = len(refseq)
        self.topsite = max(self.site) 

    def site_from_index(self,n):
        '''return site number associated with index'''
        return self.site[n]

    def index_from_site(self,site):
        '''return index associated with site number'''
        return self.ndx[site] if site <= self.topsite else self.ndx[-1]

    def indices_from_site(self,site):
        ''' return range of indices for a single site '''
        return range(self.ndx[site],self.ndx[site+1])

    def indices_from_sitelist(self,sitelist):
        '''return list of indices from list of sites'''
        ndxlist = []
        for site in sitelist:
            ndxlist.extend( self.indices_from_site(site) )
        return ndxlist


class SingleSiteMutation():
    ''' eg, D614G is a single site mutation '''
    def __init__(self,mstring=None):
        if isinstance(mstring,str):
            self.init_from_mstring(mstring)
        elif isinstance(mstring,tuple):
            self.init_from_ref_site_mut(*mstring)

    def init_from_ref_site_mut(self,ref,site,mut):
        self.ref = ref
        self.site = site
        self.mut = mut
        self.mstring = None
        
    def init_from_mstring(self,mstring):
        mstring = mstring.strip()
        m = re.match(r"(.)(\d+)(.[A-Z-_]*)",mstring)  ## what's that second '.' doing?
        if not m:
            raise RuntimeError(f"Mutation string /{mstring}/ invalid")
        self.mstring = mstring
        self.ref = m[1]
        self.site = int(m[2])
        self.mut = m[3]

        if self.mut == self.ref:
            print(self,self.mut,"==",self.ref,file=sys.stderr)

        if "_" in self.mut:
            self.mut = re.sub("_",self.ref,self.mut)

    def adjust_site(self,offset):
        self.site += offset
        self.mstring = self.ref+str(self.site)+self.mut
        return self

    def __eq__(self,other):
        '''returns boolean: is self == other?'''
        return set(other) == set(self)


    def __eq__(self,other):
        return self.ref == other.ref and self.site == other.site and self.mut == other.mut

    def __lt__(self,other):
        return self.site < other.site

    def __hash__(self):
        ''' enables the making of sets of SSM's '''
        return hash(tuple((self.ref,self.site,self.mut)))

    def __str__(self):
        return self.mstring or self.ref+str(self.site)+self.mut



                            
class Mutation(list):
    ''' Mutation is a list of SingleSiteMutations '''

    @classmethod
    def merge(cls,*muts):
        '''merge two or more mutations into a single mutation'''
        mergedmut = set()
        for mut in muts:
            mergedmut.update(mut)
        mergedmut = cls(list(mergedmut))
        return mergedmut.sort()

    def __init__(self,ssms=None):
        super().__init__(self)
        try:
            ## default translators (will be replaced when refseq is specified)
            self.site_from_index = lambda ndx: ndx+1
            self.index_from_site = lambda pos: pos-1

            ## initialization according to what ssms is
            if isinstance(ssms,list):
                self.extend(ssms)
            elif isinstance(ssms,str):
                self.init_from_line(ssms)
            elif isinstance(ssms,tuple):
                self.init_from_sequences(*ssms)
            else:
                ## just initialize as empty list
                pass
        except Exception as exception:
            raise RuntimeError(f"Unable to initialze with ssms={ssms}") from exception

    def init_from_line(self,line):
        ''' initialize Mutation from mstring '''
        mat = re.match(r".*\[(.*)\](!?).*",line)
        if mat:
            self.exact = bool(mat[2])
            for m in mat[1].split(","):
                if m.strip():
                    self.append( SingleSiteMutation(m) )
        return self

    def init_from_sequences(self,refseq,seq,exact=False):
        ''' find all sites where seq differs from refseq,
            convert into a Mutation
        '''
        ## note multisite insertions are not yet supported
        self.exact = exact
        self.set_reference_sequence(refseq)
        muts = []
        for n,(r,c) in enumerate(zip(refseq,seq)):
            if c != r:
                site = self.site_from_index(n)
                muts.append(SingleSiteMutation(f"{r}{site}{c}"))

        def flush_ssmlist(ssmlist):
            allsites = [t.site for t in ssmlist]
            if not all(s == allsites[0] for s in allsites):
                raise RuntimeError(f"assertion failed {allsites} not all equal")
            site = ssmlist[-1].site
            mstr = "".join(t.mut for t in ssmlist)
            self.append(SingleSiteMutation(f"+{site}{mstr}"))
            ssmlist.clear()

        ## Merge [-67A,-67B,-67C] -> [+67ABC]
        tmp_ssmlist=[]
        for ssm in muts:
            if ssm.ref == '-':
                if tmp_ssmlist and tmp_ssmlist[-1].site != ssm.site:
                    flush_ssmlist(tmp_ssmlist)
                    #tmp_ssmlist=[]
                tmp_ssmlist.append(ssm)
            else:
                if tmp_ssmlist:
                    flush_ssmlist(tmp_ssmlist)
                    #tmp_ssmlist=[]
                self.append(ssm)
        if tmp_ssmlist:
            flush_ssmlist(tmp_ssmlist)
            #tmp_ssmlist=[]

        return self

    def set_reference_sequence(self,refseq):
        '''use refseq to set up site_from_index and index_from_site methods'''
        T = SiteIndexTranslator(refseq)
        self.site_from_index = T.site_from_index
        self.index_from_site = T.index_from_site
        return self

    def sort(self):
        '''return a sorted copy off self'''
        return type(self)(sorted(self,key=lambda m: m.site))

    def inconsistent(self):
        ''' return an inconsistent pair of ssm's if there is one '''
        ## inconsistent if site is the same but mutation is different
        ## [T95I] and [G142D] are consistent because sites are different
        ## [T95I] and [T95T] are inconsistent, but
        ## [T95I] and [T95.] are unequal but consistent
        ## don't (yet) handle [T95_] as synonym for [T95T]
        for ssma in self:
            for ssmb in self:
                if ssma.site != ssmb.site:
                    continue
                if ssma.ref != ssmb.ref:
                    ## not only inconsistent, but one of then doesn't match refseq
                    return (str(ssma),str(ssmb))
                if ssma.mut == "." or ssmb.mut == ".":
                    continue
                if ssma.mut == "*" or ssmb.mut == "*":
                    continue
                if ssma.mut != ssmb.mut:
                    return (str(ssma),str(ssmb))
        return None

    def check(self):
        ''' check consistency '''
        return not bool(self.inconsistent())

    def checkref(self,refseq):
        ''' return True if ssm.ref = refseq at all mutation sites '''
        self.set_reference_sequence(refseq)
        check=True
        for ssm in self:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == '+':
                continue
            if ssm.ref != refseq[ndx]:
                warnings.warn(f"Mismatch for {ssm}: {ssm.ref} != {refseq[ndx]}")
                check=False
        return check

    def contains(self,other):
        '''returns boolean: does self contain other?'''
        return set(other) <= set(self)

    def contained_in(self,other):
        '''returns boolean: is self contained in other?'''
        return set(self) <= set(other)



    def pattern(self,refseq,exact=False):
        '''sort of a regex-pattern, but lite'''
        if exact is not False:
            ## we ignore 'exact'
            raise RuntimeError("mutantx.Mutation().pattern exact")
        self.set_reference_sequence(refseq)
        mutseq = list(refseq)
        if not self.exact:
            mutseq = [c if c=="-" else "." for c in refseq]
        for ssm in self:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == '+':
                mutseq[ndx] += ssm.mut
                for n in range(ndx+1,ndx+1+len(ssm.mut)):
                    if mutseq[n] == '-':
                        mutseq[n] = ""
            elif len(ssm.mut) > 1:
                mutseq[ndx] = "&"
            else:
                mutseq[ndx] = ssm.mut
        mutseq = "".join(mutseq)
        return mutseq

    def regex_pattern(self,refseq,exact=None,nodash=False):
        '''return pattern that can be used as regex to search for mutation'''
        ## exact: needs to match refseq at all sites not in mutation
        ## otherwise: only needs to match at mutation sites
        if exact is None:
            exact = self.exact

        pattern = list(refseq)
        if not exact:
            ## all .'s except -'s where the -'s are in refseq
            pattern = ['.' if r != '-' else '-' for r in refseq]

        for ssm in self:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == "+":
                pattern[ndx] += ssm.mut
                for n in range(ndx+1,ndx+1+len(ssm.mut)):
                    if exact:
                        assert pattern[n] == '-'
                    pattern[n] = ""
            elif ssm.mut == "*":
                pattern[ndx] = "[^"+ssm.ref+"]"
            elif ssm.mut == "_":
                pattern[ndx] = ssm.ref
            elif len(ssm.mut) > 1:
                pattern[ndx] = "[" + ssm.mut + "]"
            else:
                try:
                    pattern[ndx] = ssm.mut
                except IndexError as exception:
                    print("pattern=",pattern)
                    print("ndx=",ndx)
                    print("ssm:",ssm)
                    raise IndexError from exception

        pattern = "".join(pattern)
        if nodash:
            pattern = re.sub("-","",pattern)
        return pattern

    def mutate_sequence(self,refseq,prescriptive=True):
        '''return sequence that is mutated version of refseq'''
        self.set_reference_sequence(refseq)
        mutseq = list(refseq)
        for ssm in self:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == '+':
                mutseq[ndx] += ssm.mut
                for n in range(ndx+1,ndx+1+len(ssm.mut)):
                    if mutseq[n] == '-':
                        mutseq[n] = ""
            elif ssm.mut == ".":
                if prescriptive:
                    warnings.warn(f"{ssm} is not prescriptive, using {ssm.ref}")
                    mutseq[ndx] = ssm.ref
                else:
                    mutseq[ndx] = "."
            elif ssm.mut == "*":
                if prescriptive:
                    warnings.warn(f"{ssm} is not prescriptive, using X")
                    mutseq[ndx] = "X"
                else:
                    mutseq[ndx] = "*"
            elif len(ssm.mut) > 1:
                if prescriptive:
                    warnings.warn(f"{ssm} is not prescriptive, using {ssm.mut[0]}")
                    mutseq[ndx] = ssm.mut[0]
                else:
                    mutseq[ndx] = "&"
            else:
                mutseq[ndx] = ssm.mut
        mutseq = "".join(mutseq)
        return mutseq

    def __str__(self):
        return "[" + ",".join(str(ssm) for ssm in self) + "]" #+ ("!" if self.exact else "")

class MutationMaker():
    ''' Keeps track of a single refseq and SiteIndexTranslator for multiple mutations '''
    def __init__(self,refseq):
        self.refseq = refseq
        self.T = SiteIndexTranslator(refseq)

    def site_from_index(self,ndx):
        return self.T.site_from_index(ndx)

    def index_from_site(self,site):
        return self.T.index_from_site(site)

    def get_hamming(self,seq):
        '''return hamming distance from refseq to seq'''
        return sum(bool(r!=c) for r,c in zip(self.refseq,seq))

    @lru_cache(maxsize=None)
    def get_mutation(self,seq):
        '''convert sequence into a mutation list'''

        ## Scan mutations and put into dict of lists so that
        ## dict is indexed by tuple (ssm.ref,ssm.site) so that 
        ## items like [-67A,-67B,-67C] end up on the same list
        mutdict = defaultdict(list)
        for n,(r,c) in enumerate(zip(self.refseq,seq)):
            if c != r:
                site = self.site_from_index(n)
                mutdict[(r,site)].append(c)

        mutlist = []
        for (r,site),clist in mutdict.items():
            mstr = "".join(clist)
            rr = "+" if r == "-" else r
            mutlist.append(SingleSiteMutation((rr,site,mstr)))

        return Mutation(mutlist).sort()

if __name__ == "__main__":

    ma = Mutation("[S13S,S494P,N501N,D614G,P681H,T716I,S982S,D1118D]")
    mb = Mutation().init_from_line("[S13S,L141-,G142-,V143-,A222A,L452R,D614G,T859T,D950D]")
    mc = Mutation.merge(ma,mb)
    print("A:",ma)
    print("B:",mb)
    print("C inconsistent:",mc.inconsistent())
    print("C:",mc.check(),mc)
    print("C:",mc.sort())

    RefSeq = "ABC-D--EFG"
    NewSeqList = [
        "ABC-E--EFG",
        "ABCWD--EFH",
        "BBC-E--EFF",
        "BBC-E-WEFF",
        "BBC-EW-EFF",
        "ABCWEXYEFG",
        "ABCWDXYEFG",
    ]
    for newseq in NewSeqList:
        mu = Mutation().init_from_sequences(RefSeq,newseq)
        rpatt = mu.regex_pattern(RefSeq)
        rpattx = mu.regex_pattern(RefSeq,exact=True)
        print(RefSeq,newseq,rpatt,rpattx,mu)

    print()
    MM = MutationMaker(RefSeq)
    for newseq in NewSeqList:
        mu = MM.get_mutation(newseq)
        print(RefSeq,newseq,mu)
        
