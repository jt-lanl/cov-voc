'''
Mutation and SingleSiteMutation classes,
where Mutation is a list of SingleSiteMutations.
A SingleSiteMutation has three components:
  site = integer site at which mutation occurs
  ref  = character amino acid at that site in reference sequence
  mut  = character amino acid at that site in mutated sequence ('.' == any, '!' == any but ref)
'''

import sys
import re

class SingleSiteMutation():
    ''' eg, D614G is a single site mutation '''
    def __init__(self,mstring):
        mstring = mstring.strip()
        m = re.match(r"(.)(\d+)(.)",mstring)
        if not m:
            raise RuntimeError(f"Mutation string /{mstring}/ invalid")
        self.mstring = mstring
        self.ref = m[1]
        self.site = int(m[2])
        self.mut = m[3]

    def __eq__(self,other):
        return self.ref == other.ref and self.site == other.site and self.mut == other.mut

    def __hash__(self):
        ''' enables the making of sets of SSM's '''
        return hash(tuple((self.ref,self.site,self.mut)))
        
    def __str__(self):
        return self.mstring

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
        try:
            if isinstance(ssms,list):
                self.extend(ssms)
            elif isinstance(ssms,str):
                self.init_from_line(ssms)
            elif isinstance(ssms,tuple):
                self.init_from_sequences(*ssms)
            else:
                ## just initialize as empty list
                pass
        except:
            raise RuntimeError(f"Unable to initialze with ssms={ssms}")

    def init_from_line(self,line):
        muts = line.strip().strip("[]")
        for m in muts.split(","):
            if m.strip():
                self.append( SingleSiteMutation(m) )
        return self

    def init_from_sequences(self,refseq,seq):
        ''' find all sites where seq differs from refseq, 
        convert into a Mutation '''
        for n,(r,c) in enumerate(zip(refseq,seq)):
            if c != r:
                self.append(SingleSiteMutation(f"{r}{n+1}{c}"))
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
        for ssma in self:
            for ssmb in self:
                if ssma.site != ssmb.site:
                    continue
                if ssma.ref != ssmb.ref:
                    ## not only inconsistent, but one of then doesn't match refseq
                    return((str(ssma),str(ssmb)))
                if ssma.mut == "." or ssmb.mut == ".":
                    continue
                if ssma.mut == "!" or ssmb.mut == "!":
                    continue
                if ssma.mut != ssmb.mut:
                    return((str(ssma),str(ssmb)))
        return None

    def check(self):
        return False if self.inconsistent() else True
    
    def checkref(self,refseq,verbose=False):
        ''' return True if ssm.ref = refseq at all mutation sites '''
        check=True
        for ssm in self:
            if ssm.ref != refseq[ssm.site-1]:
                if verbose:
                    print(f"Mismatch for {ssm}: {ssm.ref} != {refseq[ssm.site-1]}",
                          file=sys.stderr,flush=True)
                check=False
        return check

    def pattern(self,refseq,exact=False):
        '''return pattern to describe mutation'''
        ## exact: needs to match refseq at all sites not in mutation
        ## otherwise: only needs to match at mutation sites
        pattern = list(refseq) if exact else list("." * len(refseq))
        for ssm in self:
            pattern[ssm.site-1] = ssm.mut
        return "".join(pattern)
    
    def regex_pattern(self,refseq,exact=False):
        '''return pattern that can be used as regex to search for mutation'''
        ## exact: needs to match refseq at all sites not in mutation
        ## otherwise: only needs to match at mutation sites
        pattern = list(refseq) if exact else list("." * len(refseq))
        for ssm in self:
            if ssm.mut == "!":
                pattern[ssm.site-1] = "[^"+ssm.ref+"]"
            else:
                pattern[ssm.site-1] = ssm.mut
        return "".join(pattern)

    def mutate_sequence(self,refseq):
        '''return sequence that is mutated version of refseq'''
        if any( ssm.mut == "!" for ssm in self ):
            print("mutated sequence includes !'s",file=sys.stderr,flush=True)
        return self.pattern(refseq,exact=True)

    def __str__(self):
        return "[" + ",".join(str(ssm) for ssm in self) + "]"
        
## TO-DO: an out-of-class function to read a file with mutations (also: colors, names?)

if __name__ == "__main__":

    ma = Mutation("[S13S,S494P,N501N,D614G,P681H,T716I,S982S,D1118D]")
    mb = Mutation().init_from_line("[S13S,L141-,G142-,V143-,A222A,L452R,D614G,T859T,D950D]")
    mc = Mutation.merge(ma,mb)
    print("A:",ma)
    print("B:",mb)
    print("C inconsistent:",mc.inconsistent())
    print("C:",mc.check(),mc)
    print("C:",mc.sort())

    
