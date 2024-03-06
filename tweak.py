'''module for manipulating sequence tweaks'''

import re
import itertools as it
from functools import lru_cache
import warnings
import verbose as v
import xopen
import intlist
import sequtil
import mutant
import mstringfix

_DEDASH = re.compile("-")
@lru_cache(maxsize=None)
def de_gap(seq):
    '''remove '-'s from sequence string'''
    return _DEDASH.sub("",seq)

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

    @classmethod
    def tweaks_from_filelist(cls,filelist):
        '''return tweaklist by reading multiple files'''
        tweaklist = []
        for file in filelist or []:
            tweaklist.extend( cls.tweaks_from_file(file) )
        return tweaklist

    def is_trim(self):
        '''return True if sa and sb cannot be trimmed'''
        return self.sa[0]!=self.sb[0] and self.sa[-1]!=self.sb[-1]

    def __str__(self):
        '''simple output string, does not include counts'''
        if hasattr(self,'ma') and self.ma is not None:
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

########## SITE-BASED TWEAKS (mstring pairs)

def get_mstring_pairs(mfiles,mpairs):
    '''
    Return a list of 2-tuples of m-strings;
    Each pair has a from_mstring and a to_mstring
    '''
    mstringpairs = []

    ## first: read the file(s) indicated by '-M file'
    ## A file named '.' indicates that one should use the defaults
    ## (which are currently listed in mstringfix.py)
    ## Note that '-M' can be invoked more than once
    ##   That is: '-M file1 -M file2' is okay
    ##   Invalid: '-M file1 file2'
    for mfile in mfiles or []:
        if mfile == '.':
            mfile = None
        mstringpairs.extend( mstringfix.read_mstring_pairs(mfile) )

    ## second: read from the strings on the command line
    ## these are of the from -m from_mstring to_mstring
    ## and multiple invocations of '-m' are allowed
    for mspair in mpairs or []:
        assert len(mspair)==2
        mstringpairs.append(mspair)

    ## ensure mstrings have brackets around them
    mstringpairs = [ (mstringfix.mstring_brackets(a),
                      mstringfix.mstring_brackets(b))
                     for a,b in mstringpairs ]

    ## are there any duplicates? remove them
    mstringpairs = list(dict.fromkeys(mstringpairs))

    ## Now we have all of our mstrings
    for ma,mb in mstringpairs:
        v.vvprint(f'{ma} -> {mb}')

    return mstringpairs

def get_extra_chars(mstringpairs):
    '''
    Characterize extra chars:
    xtras[site] = number of extra chars after site
    Add xtra chars for all the +nnnABC mstring component
    '''
    xtras = dict()
    for mspair in mstringpairs:
        for mstring in mspair:
            for ssm in mutant.Mutation.from_mstring(mstring):
                if ssm.ref == "+":
                    xtras[ssm.site] = max( [xtras.get(ssm.site,0),
                                            len(ssm.mut)] )
    return xtras

def add_extra_dashes(seqs,xxtras):
    '''xxtras is dict keyed by indices of s.seq strings;
       for each of those strings, we expand by extra dashes
    '''
    for s in seqs:
        sseq = list(s.seq)
        for ndx in xxtras:
            sseq[ndx] += "-"*xxtras[ndx]
        s.seq = "".join(sseq)
        yield s

def mstrings_to_ndx_seqs(mut_mgr,mstring_a,mstring_b):
    '''
    convert mstrings into short sequence-alignment fragments
    along with the indices of where those sequences start/end
    '''
    mut_a = mutant.Mutation.from_mstring(mstring_a)
    mut_b = mutant.Mutation.from_mstring(mstring_b)

    sites = sorted(set(ssm.site for ssm in it.chain(mut_a,mut_b)))
    lo,hi = sites[0],sites[-1]+1
    ndxlo = min(mut_mgr.indices_from_site(lo))
    ndxhi = max(mut_mgr.indices_from_site(hi-1))+1

    seq_r = mut_mgr.refseq
    seq_a = mut_mgr.seq_from_mutation(mut_a)
    seq_b = mut_mgr.seq_from_mutation(mut_b)

    seq_r = seq_r[ndxlo:ndxhi]
    seq_a = seq_a[ndxlo:ndxhi]
    seq_b = seq_b[ndxlo:ndxhi]

    v.vvprint(f'{mstring_a}->{mstring_b}:')
    v.vvprint(f'   r: {seq_r}')
    v.vvprint(f'   a: {seq_a}')
    v.vvprint(f'   b: {seq_b}')

    ## None of these should happen
    if seq_r == seq_a:
        warnings.warn(f"Edit {mstring_a} is same as ref sequence; "
                      "Are you sure you wan to change it !?")
    if seq_a == seq_b:
        v.vprint(f"Edit {mstring_a}->{mstring_b} will do nothing!")
    if de_gap(seq_a) != de_gap(seq_b):
        v.print(".".join(a+b for a,b in zip(de_gap(seq_a),de_gap(seq_b))))
        v.print(f'   r: {seq_r}')
        v.print(f'   a: {seq_a}')
        v.print(f'   b: {seq_b}')
        warnings.warn(f"Edit {mstring_a}->{mstring_b} "
                      "will change actual sequence!"
                      " not just the alignment\n"
                      f"  {seq_a}->{seq_b}")
        ## this shouldn't happen...but if it does, then don't do any replacing
        seq_b = seq_a

    return ndxlo,ndxhi,seq_a,seq_b

def tweak_from_mstringpair(mut_mgr,mstring_a,mstring_b):
    '''
    convert mstring pair into IndexTweak object
    '''
    tweak_tuple = mstrings_to_ndx_seqs(mut_mgr,mstring_a,mstring_b)
    tweak = IndexTweak(*tweak_tuple)
    ## add ma,mb to tweak object ... in an unclean way
    tweak.ma = mstring_a
    tweak.mb = mstring_b
    v.vprint('tweak:',str(tweak))
    return tweak

def add_needed_dashes(seqs,mstringpairs):
    '''check if mstring pairs will require extra room in the sequences;
    and if so, add the needed dashes to the sequences'''

    if not mstringpairs:
        return seqs

    ## Identify where extra columns are needed (theoretical)
    xtras = get_extra_chars(mstringpairs)
    v.vprint("xtras:",xtras)

    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    mut_mgr = mutant.MutationManager(first.seq)

    ## Identify actual columns that need to be added
    xxtras = dict() ## xx is ndx based instead of site based
    for site in sorted(xtras):
        ## determine how many extra columns there already are
        mm_xtras = len(mut_mgr.indices_from_site(site))-1
        need_xtras = xtras[site] - mm_xtras
        if need_xtras > 0:
            ## we'll need to add some columns
            v.vprint("xtras:",site,mm_xtras,"->",xtras[site])
            ## ndx is the index of the column to whch new columns should be added
            ## should we use ndx+xtras[site] instead so we are adding columns
            ## "after" the ones that are alrady there??? just thinking out loud
            ndx = mut_mgr.index_from_site(site) ## hmm,adding from left side??
            xxtras[ndx] = need_xtras

    if len(xxtras):
        seqs = add_extra_dashes(seqs,xxtras)

    return seqs,xxtras
