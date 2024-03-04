'''module for manipulating sequence tweaks'''

import re
import itertools as it
import verbose as v
import xopen
import intlist

class IndexTweak():
    '''pair of inconsistent subsequnces,
    along with their location in the full sequences
    and the count of how many sequences each appeared in'''
    ## alternative to ca,cb could make "associated_info"
    ## attribute which is a dict; could do counts (ca,cb)
    ## but could also do mstringpairs (ma,mb)

    def __init__(self,ndxlo,ndxhi,sa,sb,ca=None,cb=None):
        ndxlo = int(ndxlo)
        ndxhi = int(ndxhi)
        assert ndxhi-ndxlo == len(sa) == len(sb)
        self.ndxlo = ndxlo
        self.ndxhi = ndxhi
        self.sa = sa
        self.sb = sb
        self.ca = ca
        self.cb = cb

    @classmethod
    def from_string(cls,line):
        '''initialize from string'''
        args = line.strip().split()
        return cls(*args)

    @classmethod
    def tweaks_from_file(cls,file):
        '''read IndexTweak's from a file'''
        if not file:
            return []
        tweaklist=[]
        with xopen.xopen(file) as fpin:
            for fullline in fpin:
                ## Walrus Operator Ahead!!
                if not (line := re.sub(r'#.*','',fullline).strip()):
                    ## ignore empty lines
                    continue
                try:
                    tweak = cls.from_string(line)
                    tweaklist.append(tweak)
                except ValueError:
                    v.print(f'In file={file}, invalid line: [{fullline}]')
        return tweaklist

    def is_trim(self):
        '''return True if sa and sb cannot be trimmed'''
        return self.sa[0]!=self.sb[0] and self.sa[-1]!=self.sb[-1]

    def __str__(self):
        '''simple output string, does not include counts'''
        if hasattr(self,'ma'):
            ## if mstrings are available, prefer them for printing out
            return f'{self.ma} {self.mb}'
        return (f'{self.ndxlo} {self.ndxhi} '
                f'{self.sa} {self.sb}')

    def _as_tuple(self):
        return (self.ndxlo,self.ndxhi,self.sa,self.sb)

    def __hash__(self):
        return hash(self._as_tuple())

    def __eq__(self,other):
        return self._as_tuple() == other._as_tuple()

    def apply_to_seq(self,seq):
        '''return (seq,True) if seq tweaked; (seq,False) if not'''
        if seq[self.ndxlo:self.ndxhi] == self.sa:
            return (seq[:self.ndxlo] + self.sb + seq[self.ndxhi:],True)
        return seq,False

    def overlaps_with(self,other):
        '''Do the two tweaks regions [ndxlo:ndxhi] overlap?
        If so, return a pair of slices,
        one for self and othe for other, so that:
        self.sa[selfoverlap] corresponds to other.sa[otheroverlap]
        '''
        ndxlo = max(self.ndxlo,other.ndxlo)
        ndxhi = min(self.ndxhi,other.ndxhi)
        if ndxlo >= ndxhi:
            return False
        return (slice(ndxlo-self.ndxlo,ndxhi-self.ndxlo),
                slice(ndxlo-other.ndxlo,ndxhi-other.ndxlo))

    def interacts_with(self,other):
        '''two overlapping tweaks interact with each other if
        they agree in the overlap regaion on sa, but
        they disagree in the oerlap regsion on sb.
        (Because in that case, applying sa->sb for the two
        tweaks will lead to a different result, depending
        on the order in which they are applied.)
        '''
        if not (overlaps := self.overlaps_with(other)):
            return False
        s_olap,o_olap = overlaps
        if self.sa[s_olap] != other.sa[o_olap]:
            return False
        if self.sb[s_olap] == other.sb[o_olap]:
            return False
        return (self.sa[s_olap],self.sb[s_olap],other.sb[o_olap])

    def contains(self,other):
        '''self contains other if:
        1/ other's interval is within self's
        2/ other's sa/sb agrees with self's over that interval
        note that a tweak contains itself
        '''
        if not (self.ndxlo <= other.ndxlo and
                self.ndxhi >= other.ndxhi):
            return False
        selfslice = slice(other.ndxlo-self.ndxlo,other.ndxhi-self.ndxlo)
        if not (self.sa[selfslice] == other.sa and
                self.sb[selfslice] == other.sb):
            return False
        return True

    @staticmethod
    def get_minimal(tweaklist):
        '''from an input list of tweaks, ouptut a subset
        of them that are "minimal" meaning that they do
        not contain smaller tweaks inside themselves'''
        ## this is similar to "is_trim()" but works better
        ## in --bysite mode
        losers = set()
        for ts,to in it.combinations(tweaklist,2):
            if ts not in losers and ts.contains(to):
                losers.add(ts)
            if to not in losers and to.contains(ts):
                losers.add(to)
        return [tweak for tweak in tweaklist
                if tweak not in losers]

    def viz(self,mmgr):
        '''print a kind of viz-ualization based on actual
        site numbers'''
        #TODO: conext-y version that shows all indices
        #      over the given range of sites
        sites = [mmgr.site_from_index(ndx) for
                 ndx in range(self.ndxlo,self.ndxhi)]

        ca = f'(count={self.ca})' if self.ca else ''
        cb = f'(count={self.cb})' if self.cb else ''

        ## if ma,mb attributes available
        ma = getattr(self,'ma','')
        mb = getattr(self,'mb','')

        ref = mmgr.refseq[self.ndxlo:self.ndxhi]

        lines = [f'    {line}'
                 for line in intlist.write_numbers_vertically(sites)]
        lines.append( f' R: {ref} Reference' )
        lines.append( f' B: {self.sb} {mb} {cb}' )
        lines.append( f' A: {self.sa} {ma} {ca}' )
        return "\n".join(lines)
