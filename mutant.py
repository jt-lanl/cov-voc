'''
Mutation and SingleSiteMutation classes,
where Mutation is a list of SingleSiteMutations.
A SingleSiteMutation has three components:
  site = integer site at which mutation occurs
  ref  = character amino acid at that site in reference sequence
  mut  = character amino acid at that site in mutated sequence
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

    def __str__(self):
        return self.mstring

class Mutation(list):
    ''' Mutation is a list of SingleSiteMutations '''

    def __init__(self,ssms=None):
        if isinstance(ssms,list):
            self.extend(ssms)
        elif isinstance(ssms,str):
            self.init_from_line(ssms)
        elif isinstance(ssms,tuple):
            self.init_from_sequences(*ssms)
        else:
            ## just initialize as empty list
            pass

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
        '''return pattern that can be used as regex to search for mutation'''
        ## exact: needs to match refseq at all sites not in mutation
        ## otherwise: only needs to match at mutation sites
        pattern = list(refseq) if exact else list("." * len(refseq))
        for ssm in self:
            pattern[ssm.site-1] = ssm.mut
        return "".join(pattern)
    
    def mutate_sequence(self,refseq):
        '''return squence that is mutated versino of refseq'''
        return self.pattern(refseq,exact=True)

    def __str__(self):
        return "[" + ",".join(str(ssm) for ssm in self) + "]"
        


