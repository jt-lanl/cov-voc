'''module for manipulating sequence tweaks'''

import intlist

class IndexTweak():
    '''pair of inconsistent subsequnces,
    along with their location in the full sequences
    and the count of how many sequences each appeared in'''
    ## alternative to ca,cb could make "associated_info"
    ## attribute which is a dict; could do counts (ca,cb)
    ## but could also do mstringpairs (ma,mb)
    def __init__(self,*args):
        if len(args) == 6:
            ndxlo,ndxhi,sa,sb,ca,cb = args
            ndxlo = int(ndxlo)
            ndxhi = int(ndxhi)
        elif len(args) == 4:
            ndxlo,ndxhi,sa,sb = args
            ndxlo = int(ndxlo)
            ndxhi = int(ndxhi)
            ca = cb = None
        elif len(args) == 3:
            ndxlo,sa,sb = args
            ndxlo = int(ndxlo)
            ndxhi = ndxlo + len(sa)
            ca = cb = None
        else:
            raise ValueError('Invalid initialization')
        assert ndxhi-ndxlo == len(sa) == len(sb)
        self.ndxlo = ndxlo
        self.ndxhi = ndxhi
        self.sa = sa
        self.sb = sb
        self.ca = ca
        self.cb = cb

    def is_minimal(self):
        '''return True if sa and sb cannot be trimmed'''
        return self.sa[0]!=self.sb[0] and self.sa[-1]!=self.sb[-1]

    def __str__(self):
        '''simple output string, does not include counts'''
        if hasattr(self,'ma'):
            ## if mstrings are available, prefer them for printing out
            return f'{self.ma} {self.mb}'
        return (f'{self.ndxlo} {self.ndxhi} '
                f'{self.sa} {self.sb}')

    def __hash__(self):
        return hash((self.ndxlo,self.ndxhi,self.sa,self.sb))

    def applytweak(self,seq):
        '''return (seq,True) if seq tweaked; (seq,False) if not'''
        if seq[self.ndxlo:self.ndxhi] == self.sa:
            return (seq[:self.ndxlo] + self.sb + seq[self.ndxhi:],True)
        return seq,False

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
