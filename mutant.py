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
         x inidcates a blank (only makes sense if ref == +, and means no insertion)
         multicharacter strings indicate either:
           if ref == +, then multiple characters are inserted
           otherwise, indicates any character in the multicharacter string

'''
import re
import itertools as it
from functools import cache
from collections import defaultdict
import warnings

class SiteIndexTranslator():
    '''provide functions to translate between string index and site number
       refseq: ABC--D
       site:   123334
       index:  012345
    '''
    def __init__(self,refseq):
        self.refseq = refseq
        self.site = list(it.accumulate(int(b!='-') for b in refseq))
        if self.site:
            self.ndx = [None]*(2+max(self.site))
        else:
            self.ndx = [None]*2
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

    def refval(self,site):
        '''return the value in the reference sequence at the specified site'''
        return self.refseq[self.index_from_site(site)]

    def extra_chars_dict(self):
        '''return a dict with the number of extra characters at a given site'''
        xtrachars = dict()
        for site in range(1,1+self.topsite):
            ## = len(self.indices_from-site(site))
            nchars = self.index_from_site(site+1) - self.index_from_site(site)
            if nchars > 1:
                xtrachars[site] = nchars-1
        return xtrachars

class SingleSiteMutation():
    ''' eg, D614G is a single site mutation, so is +614ABC '''
    def __init__(self,ref,site,mut):
        self.ref=ref
        self.site=int(site)
        self.mut=mut

    @classmethod
    def from_string(cls,ssmstring):
        '''Initialize from an ssm string (eg "D614G" or "+144YGS" '''
        ## that '.' in mut-string: eg, 'x'for +123x'
        mat = re.match(r"(\D+)(\d+)(.[A-Z-_]*)",ssmstring.strip())
        if not mat:
            raise RuntimeError(f"Invalid ssm={ssmstring}")

        ref = mat[1]
        site = int(mat[2])
        mut = mat[3]

        if not re.match(r'([A-Z+])|(ins)',ref):
            raise RuntimeError(f'Invalid ref={ref} in ssm={ssmstring}')
        if ref =="ins":
            ref = "+"
        if "_" in mut:
            ## D614_ -> D614D
            ## recognize that ssm.mut might have multiple characters
            mut = re.sub("_",ref,mut)
        if mut == "x":
            if ref == "+":
                ## +123x means no insertion at 123
                ## (Q: would +123- be a better notation?)
                mut = ""
            else:
                ## D123x is not a valid formulation;
                ## use D123- for deletion at site 123
                raise RuntimeError(f"Invalid ssm={ssmstring}")

        return cls(ref,site,mut)

    def as_tuple(self):
        '''return a tuple version of the object, useful for hashing, etc'''
        return (self.ref,self.site,self.mut)

    def __eq__(self,other):
        '''
        equality test for single site mutations;
        '''
        return self.as_tuple() == other.as_tuple()

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
        '''ssmstring is, eg, 'D614G' '''
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

    @staticmethod
    def parse_mstring(mstring,exact=None):
        '''from mstring, return an ssmlist and an exact flag'''
        mat = re.match(r".*\[(.*)\](!?).*",mstring)
        if not mat:
            mstring = f'[{mstring}]'
            mat = re.match(r".*\[(.*)\](!?).*",mstring)
        if not mat:
            raise ValueError(f"Invalid mstring: {mstring}")
        if exact is None:
            exact = bool(mat[2])
        ssmlist = []
        for m in mat[1].split(","):
            if m.strip():
                ssmlist.append( SingleSiteMutation.from_string(m) )
        return ssmlist,exact


    def __init__(self,ssms=None,exact=False):
        super().__init__()
        if isinstance(ssms,str):
            cls=type(self)
            ssms,exact = cls.parse_mstring(ssms,exact=exact)
        self.extend(sorted(ssms or []))
        self.exact = exact

    @classmethod
    def from_mstring(cls,mstring,exact=None):
        ''' return a Mutation object from mstring, eg [A222V,D614G] '''
        ssmlist,exact = cls.parse_mstring(mstring,exact)
        return cls(ssmlist,exact=exact)

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
        return f'{self}{"!" if exact else ""}'

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
        xlator = SiteIndexTranslator(refseq)
        for ssm in self:
            ndx = xlator.index_from_site(ssm.site)
            if ssm.ref == '+':
                continue
            if ssm.ref != refseq[ndx]:
                warnings.warn(f"Mismatch for {ssm}: {ssm.ref} != {refseq[ndx]}")
                return False
        return True

    def relative_to(self,baseline):
        '''write the current mutation 'relative to' a baseline mutation'''
        ssms_add = [ssm for ssm in self if ssm not in baseline]
        ssms_cut = [ssm for ssm in baseline if ssm not in self]
        mutstr = ""
        if ssms_add:
            mutstr += "+" + str(Mutation(ssms_add))
        if ssms_cut:
            mutstr += "-" + str(Mutation(ssms_cut))
        return mutstr

ALL_DASHES = re.compile(r'-')
TRAILING_DASHES = re.compile(r'-+$')

class MutationManager(SiteIndexTranslator):
    '''
    Object that helps manage Mutation() mstrings,
    because it knows context; namely refseq
    (which it inherits from SiteIndexTranslator
    '''
    ## A lot of deprecated functions here now; the new/good/recommended ones are:
    ##    substr_to_mutation: uses substr_to_ssms, 
    ##    seq_to_mutation: uses subseq_to_ssms, get_ndxlists
    ##    get_mutation  (from seq): uses seq_to_mutation (indeed, kind of redundant with it)
    ##    regex_from_mstring: uses regex_from_mutation
    ##    regex_from_mutation

    def substr_to_ssms(self,gseq,ndxlo=0,insertwithdashes=True,includeidentityssms=True):
        '''
        given a gappy substring gseq (such as 'R---T-T-')
        and the index associated with the first character
        yield ssm's that would make up an mstring
        call: Mutation(substr_to_ssms(...)) to get Mutation object
        '''
        re_dashes = TRAILING_DASHES if insertwithdashes else ALL_DASHES
        ndxhi = ndxlo + len(gseq)
        while ndxlo < ndxhi:
            site = self.site_from_index(ndxlo)
            ndxlolo = self.index_from_site(site)
            if ndxlolo == ndxlo:
                ## substitution (or deletion) character
                subchr,gseq = gseq[0],gseq[1:]
                if includeidentityssms or self.refseq[ndxlo] != subchr:
                    yield SingleSiteMutation(self.refseq[ndxlo],site,subchr)
                ndxlo = ndxlo + 1
            elif ndxlolo < ndxlo:
                ## insertion string
                ndx_next = self.index_from_site(site+1)
                insstr,gseq = gseq[:ndx_next-ndxlo],gseq[ndx_next-ndxlo:]
                insstr = re_dashes.sub('',insstr) ## remove (trailing) dashes
                if insstr:
                    yield SingleSiteMutation("+",site,insstr)
                ndxlo = ndx_next
            else: # ndxlolo > ndxlo:
                raise RuntimeError(f'{ndxlolo=} > {ndxlo=}: should never happen')

    def substr_to_mutation(self,gseq,ndxlo=0,insertwithdashes=True):
        '''return mutation object based on gapped subsequence'''
        fullseq = bool(ndxlo==0 and len(gseq) == len(self.refseq))
        mut = Mutation(self.substr_to_ssms(gseq,ndxlo,
                                           includeidentityssms = not fullseq,
                                           insertwithdashes=insertwithdashes))
        mut.exact = fullseq
        #v.print('substr_to_mut: mut=',str(mut))
        return mut

    @cache # pylint: disable=method-cache-max-size-none
    def subseq_to_ssms(self,gseq,ndxlo,insertwithdashes=True):
        '''cached version of substr_to_ssms'''
        return list(self.substr_to_ssms(gseq,ndxlo,
                                        insertwithdashes=insertwithdashes,
                                        includeidentityssms=False))

    @cache # pylint: disable=method-cache-max-size-none
    def get_ndxlists(self,nchunks):
        '''return two lists, for site-aligned ndxlo and ndxhi used for slicing sequence into chunks'''
        sites = [1+n*self.topsite//nchunks for n in range(nchunks)]
        ndxlo_list = [self.index_from_site(site) for site in sites]
        ndxhi_list = ndxlo_list[1:] + [len(self.refseq)]
        return ndxlo_list,ndxhi_list
        
    def seq_to_mutation(self,seq,exact=True):
        '''convert sequence to Mutation object, using chunked strategy for speed'''
        ## Using chunks means subseq_to_ssms only works with small subseq lengths,
        ## which makes the arguments more cache-able
        nchunks = min(64,self.topsite//2)
        los,his = self.get_ndxlists(nchunks)
        ssmlist = []
        for ndxlo,ndxhi in zip(los,his):
            ssms = self.subseq_to_ssms(seq[ndxlo:ndxhi],ndxlo)
            #print(ndxlo,ndxhi,seq[ndxlo:ndxhi],[str(ssm) for ssm in ssms])
            ssmlist.extend(ssms)
        mut = Mutation(ssmlist)
        mut.exact = exact
        return mut

    #TODO: deprecated
    @cache # pylint: disable=method-cache-max-size-none
    def get_partial_mutdict(self,subseq,ndxlo,ndxhi):
        '''cached helper function used by alt_get_ssmlist'''
        partial_mutdict = defaultdict(list)
        for n,(r,c) in enumerate(zip(self.refseq[ndxlo:ndxhi],subseq)):
            if c != r:
                site = self.site_from_index(ndxlo+n)
                rr = "+" if r == "-" else r
                partial_mutdict[(rr,site)].append(c)
        return partial_mutdict

    #TODO: deprecated
    def get_mutdict(self,seq):
        '''helper function used by alt_get_ssmlist'''
        nchunks=64
        lolist = [len(self.refseq)*n//nchunks     for n in range(nchunks)]
        hilist = [len(self.refseq)*(n+1)//nchunks for n in range(nchunks)]
        mutdict = defaultdict(list)
        for lo,hi in zip(lolist,hilist):
            tmp_mutdict = self.get_partial_mutdict(seq[lo:hi],lo,hi)
            for (rr,site),clist in tmp_mutdict.items():
                mutdict[(rr,site)].extend(clist)
        return mutdict

    #TODO: deprecated
    def alt_get_ssmlist(self,seq):
        '''alternative get_ssmlist: seq -> list of ssm's '''
        ## Can be faster because get_partial_mutdict can be more cache-able
        ## because the subseq is likely to be a repeated argument (especially
        ## when ndxhi-ndlxo is small
        mutlist = []
        mutdict = self.get_mutdict(seq)
        for (rr,site),clist in mutdict.items():
            mstr = "".join(clist)
            mutlist.append(SingleSiteMutation(rr,site,mstr))
        return mutlist

    #TODO: deprecated
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
            mutlist.append(SingleSiteMutation(rr,site,mstr))
        return mutlist ## should already be sorted?

    # Sort-of deprecated, still widely used, but seq_to_mutation does the same thing
    def get_mutation(self,seq,exact=True):
        '''convert sequence into a mutation list'''
        ## since obtained from seq, assume mseq has exact=True
        ## OLD deprecated method
        ## mutlist = self.alt_get_ssmlist(seq)
        ## mut = Mutation(mutlist)
        mut = self.seq_to_mutation(seq)
        mut.exact = exact
        return mut

    #TODO deprecated
    def get_alt_mutation_ssmlist(self,ssmlist):
        '''convert list of ssms into an alt mutation structure'''
        ## alt mutation structure is a dictionary of sites
        ## [A1B,C3D,+3YZ] -> {1: 'B', 3: 'DYZ'}
        ## [A1B,+3YZ] -> {1: 'B', 3:'CYZ'}
        ## [A1*] -> {1: '*A'}
        ## alt structure designed for describing sequences;
        ## if ssmlist is meant to be a pattern, then more likely to be problematic
        alt=dict()
        for ssm in ssmlist:
            if ssm.ref == "+":
                alt[ssm.site] = alt.get(ssm.site,self.refval(ssm.site)) + ssm.mut
            elif ssm.mut == "*":
                alt[ssm.site] = "*" + ssm.ref
            else:
                alt[ssm.site] = ssm.mut
        return alt

    #TODO deprecated
    def get_alt_mutation(self,seq):
        '''convert sequence into mutation dict'''
        mutlist = self.alt_get_ssmlist(seq)
        return self.get_alt_mutation_ssmlist(mutlist)


    def regex_from_mstring(self,mstring,**kw):
        '''return a string that can be used as a simple regex
        that corresopnds to the mstring'''
        if isinstance(mstring,Mutation):
            raise ValueError(f'mstring argument is Mutation, not string; '
                             f'use regex_from_mutation(...) instead')
        mut = Mutation.from_mstring(mstring)
        return self.regex_from_mutation(mut,**kw)

    def regex_from_mutation(self,mut,exact=None,compile=False):
        '''return a string that can be used as a simple regex'''
        ## similar to seq_from_mutation, especially if exact=True
        ## but doen't do special characters
        if isinstance(mut,str):
            raise ValueError(f'mut is a string, not a Mutation; '
                             f'use regex_from_mstring(...) instead')
        if exact is None:
            exact = mut.exact
        regex_as_list = list(self.refseq) if exact else ["."] * len(self.refseq)
        for ssm in mut:
            ndxlo = self.index_from_site(ssm.site)
            if ssm.ref == "+":
                ndxhi = self.index_from_site(ssm.site+1)
                if len(ssm.mut) >= ndxhi-ndxlo:
                    raise ValueError(f'Insertion {ssm} too long; '
                                     f'only room for {ndxhi-ndxlo-1} characters')
                regex_as_list[ndxlo+1:ndxlo+1+len(ssm.mut)] = ssm.mut
                continue
            if ssm.ref != self.refseq[ndxlo]:
                raise ValueError(f'{ssm=!s} ref {ssm.ref} disagrees with '
                                 f'ref={self.refseq[ndxlo]} at site {ndxlo}')
            if len(ssm.mut) != 1:
                raise ValueError(f'{ssm=!s} invalid: {len(ssm.mut)=} '
                                 f'must have exactly one character')
            ## possible extension: if ssm.mut == "*", then regex_as_list[ndxlo] = f'[^{ssm.ref}]'
            if ssm.mut == "*":
                ## asterisk means anything except ref
                regex_as_list[ndxlo] = f'[^{ssm.ref}]'
                continue
            regex_as_list[ndxlo] = ssm.mut

        regex = "".join(regex_as_list)
        if compile:
            regex = re.compile(regex)
        return regex

    ## Still used for spikevariants, but there's probably a better way
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
                    raise RuntimeError(f"Error in mut={mut}, ssm={ssm}: "
                                       f"{ssm.ref}!={self.refseq[ndx]}")
                if len(ssm.mut) > 1:
                    mutseq[ndx] = "&"  ## indicates multiple choices at this site
                else:
                    mutseq[ndx] = ssm.mut
        mutseq = "".join(mutseq)
        return mutseq

    ## Are you sure you wouldn't prefer regex_from_mutation(..,exact=True) ?
    def seq_from_mutation(self,mut,prescriptive=True):
        '''return sequence that is mutated version of refseq'''
        ## note: insertions at the head of the sequence will
        ## trigger an error; eg, [+0ABC]
        mutseq = list(self.refseq)
        for ssm in mut:
            ndx = self.index_from_site(ssm.site)
            if ssm.ref not in ['+', self.refseq[ndx]]:
                ## complain but for now, don't exit
                warnings.warn(f'at site {ssm.site}, ref={self.refseq[ndx]} != {ssm}')
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
                except IndexError as err:
                    errmsg = \
                        f'ssm={ssm},site={ssm.site},topsite={self.topsite},' \
                        f'ndx={ndx},len(mutseq)={len(mutseq)}'
                    raise IndexError(errmsg) from err
        mutseq = "".join(mutseq)
        return mutseq

    ## Now much faster to use regex_from_mutation, followed by re.match(...)
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

    ## Now much faster to use regex_from_mutation, followed by re.match(...)
    def filter_seqs_by_pattern(self,mpatt,seqs,exact=None):
        '''filter an iterable of seqs that match the pattern'''
        ## if input mpatt is string, this converts to Mutation object
        ## also, makes sure the ssm's in mpatt are sorted
        mpatt = Mutation(sorted(Mutation(mpatt)))
        for s in seqs:
            if self.seq_fits_pattern(mpatt,s.seq,exact=exact):
                yield s
