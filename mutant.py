'''
MutationManager, Mutation, and SingleSiteMutation classes
Also, the SiteIndexTranslator class

A Mutation is a list of SingleSiteMutations, with auxiliary information such as an 'exact' flag.
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
        ## if site too large, just return last index
        return self.ndx[site] if site <= self.topsite else self.ndx[-1]

    def indices_from_site(self,site):
        ''' return range of indices for a single site '''
        return range(self.index_from_site(site),
                     self.index_from_site(site+1))

    def indices_from_sitelist(self,sitelist):
        '''return list of indices from list of sites'''
        ndxlist = []
        for site in sitelist:
            ndxlist.extend( self.indices_from_site(site) )
        return ndxlist

    def extra_chars_dict(self):
        '''return a dict with the number of extra characters at a given site'''
        xtrachars = dict()
        for site in range(1,1+self.topsite):
            nchars = self.index_from_site(site+1) - self.index_from_site(site)
            if nchars > 1:
                xtrachars[site] = nchars-1
        return xtrachars

class SingleSiteMutation():
    ''' eg, D614G is a single site mutation, so is +614ABC '''
    def __init__(self,initializer=None):
        if isinstance(initializer,str):
            self.init_from_mstring(initializer)
        elif isinstance(initializer,tuple):
            self.init_from_ref_site_mut(*initializer)

    @classmethod
    def from_ref_site_mut(cls,ref,site,mut):
        ssm = SingleSiteMutation()
        ssm.ref = ref
        ssm.site = site
        ssm.mut = mut
        return ssm
        
    def init_from_ref_site_mut(self,ref,site,mut):
        '''initialize from ref char, site number, and mutation char'''
        self.ref = ref
        self.site = site
        self.mut = mut
        return self

    def init_from_mstring(self,mstring):
        '''initialize from mstring; eg "D614G" '''
        mstring = mstring.strip()
        m = re.match(r"(.)(\d+)(.[A-Z-_]*)",mstring)  ## what's that second '.' doing?
        if not m:
            raise RuntimeError(f"SingleSiteMutation string /{mstring}/ invalid")
        self.ref = m[1]
        self.site = int(m[2])
        self.mut = m[3]

        if self.mut == self.ref:
            ## this is really ok, do I need to complain about it?
            ## print(self,self.mut,"==",self.ref,file=sys.stderr)
            pass

        if "_" in self.mut:
            ## recognize that ssm.mut might have multiple characters
            self.mut = re.sub("_",self.ref,self.mut)

        return self

    def as_tuple(self):
        '''return a tuple version of the object, useful for hashing, etc'''
        return (self.ref,self.site,self.mut)

    def __eq__(self,other):
        '''
        equality test for single site mutations;
        '''
        if self.as_tuple() == other.as_tuple():
            return True
        return False

    def xxx(self,other):
        '''the rest of the __eq__ function, 
        in case you want [R256.] and [R256A] considered equal
        '''
        if self.ref != other.ref:
            return False
        if self.site == other.site:
            if self.mut == other.mut:
                return True
            if self.mut == "." or other.mut == ".":
                return True
            if self.mut == "*" and other.mut != other.ref:
                return True
            if self.mut != self.ref and other.mut == "*":
                return True
        return False

    def __lt__(self,other):
        ''' less-than function enables sorting of SSM's '''
        ## enables sorting lists of SSM's by site number
        ## but if two ssm's have the same site number, put the "+" after
        ## eg sorted([A7B,+3FF,D3G,C4A]) --> [D3G,+3FF,C4A,A7B]
        return (self.site,self.ref == '+') < (other.site,other.ref == '+')

    def __hash__(self):
        ''' enables the making of sets of SSM's, ability to use SSM's as dictionary keys, etc '''
        return hash(self.as_tuple())

    def __str__(self):
        '''mstring is, eg, 'D614G' '''
        return self.ref+str(self.site)+self.mut

class Mutation(list):
    ''' Mutation is a list of SingleSiteMutations '''

    ## Several variants of a mutation: mpatt vs mseq
    ## mpatt is a "pattern"
    ##      mpatt may or may not include wildcards
    ##      mpatt is pattern against which many distinct sequences might match
    ##      exact=False  (fullmatch=False)
    ## mseq is descriptive of a sequence;
    ##      mseq should not include wildcards, and
    ##      mseq implicitly asserts that sites not mentioned match a refseq[*]
    ##      exact=True   (fullmatch=True)
    ##
    ## [*] note, however, that 'refseq' is not part of a Mutation; thus, pattern matching
    ##       (exact pattern matching, in any case) cannot be done by Mutation alone,
    ##       the class MutationManager is designed for those tasks
    ##
    ## does a mutation need to know if it is an mpatt or an mseq?
    ##      yes: because when it reads an mstring that string might include '!' indicating exact
    ##

    @classmethod
    def merge(cls,*muts):
        '''merge two or more mutations into a single mutation'''
        ## does this ever actually get used???
        ## and what does it mean to merge T95I with T95. (or T95K)
        ## here, if any of muts are exact, then merged mut is exact -- does that make sense?
        mergedmut = set()
        exact = False
        for mut in muts:
            mergedmut.update(mut)
            exact = exact | mut.exact
            print("merge:",mut,mut.exact,exact)
        mergedmut = cls(sorted(mergedmut))
        mergedmut.exact = exact
        return mergedmut

    def __init__(self,ssms=None,exact=False):
        super().__init__(self)
        self.exact = exact
        try:
            ## initialization according to what ssms is
            if isinstance(ssms,list):
                self.extend(ssms)
            elif isinstance(ssms,set):
                self.extend(sorted(ssms))
            elif isinstance(ssms,str):
                self.init_from_mstring(ssms)
            elif isinstance(ssms,tuple):
                self.init_from_sequences(*ssms)
            elif ssms:
                raise ValueError(f"Invalid argument initializaing {type(self)}")
        except Exception as exception:
            raise RuntimeError(f"Unable to initialze with ssms={ssms}") from exception

    @classmethod
    def from_mstring(cls,mstring,exact=None):
        ''' initialize Mutation from mstring, eg [A222V,D614G] '''
        return cls().init_from_mstring(mstring,exact=exact)
        
    def init_from_mstring(self,mstring,exact=None):
        ''' initialize Mutation from mstring, eg [A222V,D614G] '''
        mat = re.match(r".*\[(.*)\](!?).*",mstring)
        if not mat:
            raise ValueError(f"Invalid mstring: {mstring}")
        if exact is None:
            exact = bool(mat[2])
        self.exact = exact
        ssmlist = []
        for m in mat[1].split(","):
            if m.strip():
                ssmlist.append( SingleSiteMutation(m) )
        self.extend(sorted(ssmlist))
        return self

    def init_from_sequences(self,refseq,seq,exact=False):
        ''' find all sites where seq differs from refseq,
            convert into a Mutation
        '''
        warnings.warn("not efficient, use MutationManager")
        MM = MutationManager(refseq)
        self.extend( MM.get_ssmlist(seq) )
        self.exact = exact
        return self

    #def __eq__(self,other):
    #    '''returns boolean: is self == other?'''
    #    ## note, if we are careful that Mutation is sorted, then 'sorted' here not necessary
    #    return sorted(other) == sorted(self)

    def __str__(self):
        return "[" + ",".join(str(ssm) for ssm in self) + "]" # + ("!" if self.exact else "")

    def as_string(self,exact=None):
        '''instead of str(mut), use mut.as_string() to incorporate the exact ('!') flag'''
        if exact is None:
            exact = self.exact
        return self.__str__() + ("!" if exact else "")

    ## note, no hash function, just too expensive to hash

    def sitelist(self):
        '''
        return list of relevant sites, with mutliple copies for insertions
        eg, [A1B,C5.,+7AB] --> 1,5,7,7,7
        '''

        sites=[]
        for ssm in self:
            sites.append(ssm.site)
            if ssm.ref == "+": ## then append a few more sites
                for _ in range(len(ssm.mut)-1):
                    sites.append(ssm.site)
        return sites

    def inconsistent(self):
        ''' return an inconsistent pair of ssm's if there is one '''
        ## on further reflection, the only time sites should be identical
        ## is with insertions; eg R214A,+214AGY
        for ssma,ssmb in it.combinations(self, 2):
            if ssma.site != ssmb.site:
                continue
            if (ssma.ref=='+') + (ssmb.ref=='+') == 1:
                ## ok only if exactly one of the two is a '+'
                continue
            return (str(ssma),str(ssmb))
        return None

    def check(self):
        ''' check consistency '''
        return not bool(self.inconsistent())

    ### routines below are potentially deprecated, to be
    ### replaced with MutationManager methods

    def checkref(self,refseq):
        ''' return True if ssm.ref = refseq at all mutation sites '''
        T = SiteIndexTranslator(refseq)
        for ssm in self:
            ndx = T.index_from_site(ssm.site)
            if ssm.ref == '+':
                continue
            if ssm.ref != refseq[ndx]:
                warnings.warn(f"Mismatch for {ssm}: {ssm.ref} != {refseq[ndx]}")
                return False
        return True

    def pattern(self,refseq,exact=None):
        '''
        returns a string of same length as refseq,
        but with substitutions, wildcards, etc noted
        note, it's NOT a regex pattern matcher
        '''
        raise RuntimeError("Use MutationManager.pattern_from_mutation instead")


    def regex_pattern(self,refseq,exact=None,nodash=False):
        '''return pattern that can be used as regex to search for mutation'''
        ## exact: needs to match refseq at all sites not in mutation
        ## otherwise: only needs to match at mutation sites
        raise RuntimeError("Not using regex's anymore, are we?")

    def mutate_sequence(self,refseq,prescriptive=True):
        '''return sequence that is mutated version of refseq'''
        raise RuntimeError("Use MutationManager.seq_from_mutation")

    
class MutationManager(SiteIndexTranslator):
    '''
    Object that helps manage Mutation() mstrings, 
    because it knows context; namely refseq
    '''
    def __init__(self,refseq):
        super().__init__(refseq)
        self.refseq = refseq

    def refval(self,site):
        '''return the value in the reference sequence at the specified site'''
        return self.refseq[self.index_from_site(site)]

    @lru_cache(maxsize=None)
    def get_ssmlist(self,seq):
        '''convert sequence into a list of ssm's '''
        ## Scan mutations and put into dict of lists so that
        ## dict is indexed by tuple (ssm.ref,ssm.site) so that
        ## items like [-67A,-67B,-67C] end up on the same list
        mutdict = defaultdict(list)
        for n,(r,c) in enumerate(zip(self.refseq,seq)):
            if c != r:
                site = self.site_from_index(n)
                mutdict[(r,site)].append(c)
        ## now unwind the mutdict and produce a list of ssms
        mutlist = []
        for (r,site),clist in mutdict.items():
            mstr = "".join(clist)
            rr = "+" if r == "-" else r
            mutlist.append(SingleSiteMutation.from_ref_site_mut(rr,site,mstr))
        return mutlist ## should already be sorted? 

    def get_mutation(self,seq,exact=True):
        '''convert sequence into a mutation list'''
        ## since obtained from seq, assume mseq has exact=True
        mutlist = self.get_ssmlist(seq)
        mut = Mutation(mutlist)
        mut.exact = exact
        return mut

    def get_alt_mutation_ssmlist(self,ssmlist):
        '''convert list of ssms into an alt mutation structure'''
        ## alt mutation structure is a dictionary of sites
        ## [A1B,C3D,+3YZ] -> {1: 'B', 3: 'DYZ'}
        ## [A1B,+3YZ] -> {1: 'B', 3:'CYZ'}
        ## [A1*] -> {1: '*A'}
        ## alt structure designed for describing sequences;
        ## if ssmlist is meant to be a pattern, then more likely to be problematic
        alt=dict()
        for ssm in sorted(ssmlist):
            if ssm.ref == "+":
                alt[ssm.site] = alt.get(ssm.site,self.refval(ssm.site)) + ssm.mut
            elif ssm.mut == "*":
                alt[ssm.site] = "*" + ssm.ref
            else:
                alt[ssm.site] = ssm.mut
        return alt

    def get_alt_mutation(self,seq):
        '''convert sequence into mutation dict'''
        mutlist = self.get_ssmlist(seq)
        return self.get_alt_mutation_ssmlist(mutlist)

    def pattern_from_mutation(self,mut,exact=None):
        '''
        returns a string of same length as refseq,
        but with substitutions, wildcards, etc
        note, it's NOT a regex pattern matcher
        '''
        ## very similar to seq_from_mutation,
        ## can we factor out common elements
        if exact is None:
            exact = mut.exact
        mutseq = list(self.refseq)
        if not exact:
            mutseq = [c if c=="-" else "." for c in self.refseq]
        for ssm in mut:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == '+':
                mutseq[ndx] += ssm.mut
                for n in range(ndx+1,ndx+1+len(ssm.mut)):
                    if mutseq[n] == '-':
                        mutseq[n] = ""
            else:
                if ssm.ref != self.refseq[ndx]:
                    raise RuntimeError(f"Error in mut={mut}, ssm={ssm}: {ssm.ref}!={self.refseq[ndx]}")
                if len(ssm.mut) > 1:
                    mutseq[ndx] = "&"  ## indicates multiple choices at this site
                else:
                    mutseq[ndx] = ssm.mut
        mutseq = "".join(mutseq)
        return mutseq

    def seq_from_mutation(self,mut,prescriptive=True):
        '''return sequence that is mutated version of refseq'''
        ## note: insertions at the head of the sequence will
        ## trigger an error; eg, [+0ABC]
        mutseq = list(self.refseq)
        for ssm in mut:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref == '+':
                mutseq[ndx] += ssm.mut
                for n in range(ndx+1,ndx+1+len(ssm.mut)):
                    if mutseq[n] == '-':
                        mutseq[n] = ""
                    else:
                        warnings.warn(f"ssm=[{ssm}] has too many insertion characters")
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
                try:
                    mutseq[ndx] = ssm.mut
                except IndexError as e:
                    print("ssm,site,topsite,ndx,len(mutseq)=",
                          ssm,ssm.site,self.topsite,
                          ndx,len(mutseq))
                    raise IndexError(e)
        mutseq = "".join(mutseq)
        return mutseq

    ## This is one very intricate function!!! so many special cases!!
    ## any way to offload some of the cognitive budget? not to mention raw compute time
    ## if mpatt has no wild cards, this could go faster (separate routine for that?)
    def seq_fits_pattern(self,mpatt,seq,exact=None):
        '''return True if the sequence fits the pattern'''
        #if not all( m==sm for m,sm in zip(mpatt,sorted(mpatt)) ):
        #    raise RuntimeError(f"unsorted mpatt: {mpatt}")
        if exact is None:
            exact = mpatt.exact
        mseq = self.get_mutation(seq)
        if mseq == mpatt: ## easy case: if seq and patt are identical, then True
            return True
        alt_mseq = self.get_alt_mutation_ssmlist(mseq)
        for ssm in mpatt:
            refval = self.refval(ssm.site)
            seqval = alt_mseq.get(ssm.site,refval)
            if ssm.mut == "_":
                ## may never happen, get's converted when read into SSM
                if seqval != ssm.ref:
                    return False
            elif ssm.mut == "*":
                if seqval[0] == ssm.ref or seqval[0] == ".":
                    return False
            elif ssm.mut == ".":
                pass
            elif ssm.ref == "+":
                if ssm.mut != seqval[1:]:
                    return False
            else:
                if ssm.mut != seqval[0]:
                    return False
        if not exact:
            return True
        ## made it this far; now in exact mode
        alt_mpatt = self.get_alt_mutation_ssmlist(mpatt)
        for site,seqval in alt_mseq.items():
            pval = alt_mpatt.get(site,None)
            if not pval:
                return False
            ref = self.refval(site)
            seqval = seqval or ref
            if pval[0] == "*" and seqval != ref: ## ref = pval[1] in this case
                continue
            if pval[0] == ".":
                continue
            if pval != seqval:
                return False
        return True

    def filter_seqs_by_pattern(self,mpatt,seqs,exact=None):
        '''filter an iterable of seqs that match the pattern'''
        ## if input mpatt is string, this converts to Mutation object
        ## also, makes sure the ssm's in mpatt are sorted
        mpatt = Mutation(sorted(Mutation(mpatt))) 
        for s in seqs:
            if self.seq_fits_pattern(mpatt,s.seq,exact=exact):
                yield s
        

if __name__ == "__main__":

    ma = Mutation.from_mstring("[S13S,S494P,N501N,D614G,P681H,T716I,S982S,D1118D]")
    mb = Mutation.from_mstring("[S13S,L141-,G142-,V143-,A222A,L452R,D614G,T859T,D950D]")
    mc = Mutation.merge(ma,mb)
    print("A:",ma)
    print("B:",mb)
    print("C:",mc)
    print("C inconsistent:",mc.inconsistent())
    print("C:",mc.check(),mc)
    print("C:",mc)
    print("C:",
          "Note, this is a list of ssms, not a Mutation:",
          sorted(mc))

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
    MM = MutationManager(RefSeq)
    for newseq in NewSeqList:
        mu = MM.get_mutation(newseq)
        fpatt = MM.pattern_from_mutation(mu)
        print(RefSeq,newseq,fpatt,mu)

    print()
    MM = MutationManager(RefSeq)
    for newseq in NewSeqList:
        mu = MM.get_mutation(newseq)
        print(RefSeq,newseq,mu)
    print()
    for patt in [
            '[A1B,D4E,+4W,G7F]',
            '[D4E]',
            '[+3W,G7H]',
            '[A1B,D4E,G7F]',
            '[+3W,D4E,+4XY]',
            '[+3W,+4XY]',
            '[D4*]',
            '[A1_,D4.]',
            '[+3W,C3_,G7H]',
            '[A1_,C3.]',
            ]:
        thepatt = Mutation().from_mstring(patt)
        for newseq in NewSeqList:
            fit = MM.seq_fits_pattern(thepatt,newseq)
            fullfit = MM.seq_fits_pattern(thepatt,newseq,exact=True)
            #print(thepatt,newseq,fit,fullfit)
            #print(f"assert MM.seq_fits_pattern(Mutation().init_from_mstring('{thepatt}'),'{newseq}')            == {fit}")
            #print(f"assert MM.seq_fits_pattern(Mutation().init_from_mstring('{thepatt}'),'{newseq}',exact=True) == {fullfit}")

    def checkit(mstring,seq,**kw):
        mpatt = Mutation.from_mstring(mstring)
        return MM.seq_fits_pattern(mpatt,seq,**kw)
    
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
