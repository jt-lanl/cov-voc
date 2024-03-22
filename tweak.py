'''module for manipulating sequence tweaks'''

import re
import itertools as it
from functools import cache
import warnings
import verbose as v
import xopen
import intlist
import sequtil
import mutant
import mstringfix

_DEDASH = re.compile("-")
@cache
def de_gap(seq):
    '''remove '-'s from sequence string'''
    return _DEDASH.sub("",seq)

class IndexTweak(): # pylint: disable=too-many-instance-attributes
    '''pair of inconsistent subsequnces,
    along with their location in the full sequences
    and the count of how many sequences each appeared in'''
    ## alternative to ca,cb could make "associated_info"
    ## attribute which is a dict; could do counts (ca,cb)
    ## but could also do mstringpairs (ma,mb)

    def __init__(self,ndxlo,ndxhi,sa,sb,ca=None,cb=None): # pylint: disable=too-many-arguments
        ndxlo = int(ndxlo)
        ndxhi = int(ndxhi)
        assert ndxhi-ndxlo == len(sa) == len(sb)
        self.ndxlo = ndxlo
        self.ndxhi = ndxhi
        self.sa = sa
        self.sb = sb
        self.ca = ca
        self.cb = cb
        self.ma = self.mb = None

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

    def swap_ab(self):
        '''swap sa <-> sb; also: ca <-> cb and ma <-> mb'''
        self.sa,self.sb = self.sb,self.sa
        self.ca,self.cb = self.cb,self.ca
        self.ma,self.mb = self.mb,self.ma

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
        2/ other's sa/sb agrees with self's sa/sb (or sb/sa) over that interval
        note that a tweak contains itself
        '''
        if not (self.ndxlo <= other.ndxlo and
                self.ndxhi >= other.ndxhi):
            return False
        selfslice = slice(other.ndxlo-self.ndxlo,other.ndxhi-self.ndxlo)
        if (self.sa[selfslice] == other.sa and
            self.sb[selfslice] == other.sb):
            return True
        if (self.sa[selfslice] == other.sb and
            self.sb[selfslice] == other.sa):
            return True
        return False

    def update_mstringpair(self,xlator):
        '''update attributes ma/mb with mstrings associated with sa/sb'''
        self.ma = str(substr_to_mut(xlator,self.sa,self.ndxlo))
        self.mb = str(substr_to_mut(xlator,self.sb,self.ndxlo))
        return self.sa,self.sb

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

    def viz(self,xlator,showcontext=False):
        '''print a kind of viz-ualization based on actual site numbers'''

        # counts (blank if zero or None)
        ca = f'(count={self.ca})' if self.ca else ''
        cb = f'(count={self.cb})' if self.cb else ''

        ## if ma,mb attributes available (blank if not)
        ma = getattr(self,'ma','')
        mb = getattr(self,'mb','')

        sites = [xlator.site_from_index(ndx) for
                 ndx in range(self.ndxlo,self.ndxhi)]

        ssb = self.sb
        ssa = self.sa
        ref = xlator.refseq[self.ndxlo:self.ndxhi]
        ctx_ndxlo = self.ndxlo
        if showcontext:
            ctx_ndxlo = xlator.index_from_site(min(sites))
            if ctx_ndxlo < self.ndxlo:
                sites = [xlator.site_from_index(ndx) for
                         ndx in range(ctx_ndxlo,self.ndxhi)]
                ssb = " "*(self.ndxlo-ctx_ndxlo) + self.sb
                ssa = " "*(self.ndxlo-ctx_ndxlo) + self.sa
                ref = xlator.refseq[ctx_ndxlo:self.ndxhi]

        ndx_lines = [f'n   {line}'
                     for line in intlist.write_numbers_vertically(range(ctx_ndxlo,self.ndxhi))]
        site_lines = [f's   {line}'
                      for line in intlist.write_numbers_vertically(sites)]
        lines = sum([ndx_lines,site_lines],[])
        lines.append( f' R: {ref} Reference' )
        lines.append( f' B: {ssb} {mb} {cb}' )
        lines.append( f' A: {ssa} {ma} {ca}' )
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
    #seq_a = mut_mgr.seq_from_mutation(mut_a)
    #seq_b = mut_mgr.seq_from_mutation(mut_b)
    seq_a = mut_mgr.regex_from_mutation(mut_a,exact=True)
    seq_b = mut_mgr.regex_from_mutation(mut_b,exact=True)

    if seq_a == seq_b:
        warnings.warn(f'Tweak {mstring_a}->{mstring_b} will not change anything!')
        return 0,0,"",""


    while seq_a[ndxlo] == seq_b[ndxlo]:
        ## this can occur if one of the mstrings is of the form [+251V],
        ## so we want to start not at ndxlo associated with site=251
        ## but with ndxlo one higher than that.
        ndxlo += 1
    while seq_a[ndxhi-1] == seq_b[ndxhi-1]:
        ## this can occur if one of the mstrings ends with an insertion: [...,+251V]
        ndxhi -= 1

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
    return tweak

def expansion_needed(mut_mgr,mstringpairs):
    '''check if the mstring pairs will require extra room in the sequences'''
    if not mstringpairs:
        return False

    ## Identify where extra columns are needed (based just on mstrings)
    columns_needed = get_extra_chars(mstringpairs)
    v.vprint("columns needed:",columns_needed)

    ## Identify actual columns that need to be added/expanded
    xpand = dict() ## xx is ndx based instead of site based
    for site in columns_needed:
        ## determine how many columns there already are
        columns_available = len(mut_mgr.indices_from_site(site))-1
        need_xtras = columns_needed[site] - columns_available
        if need_xtras <= 0:
            continue
        ## we'll need to add some columns
        v.vprint("columns expanded:",site,columns_available,"->",columns_needed[site])
        ## ndx is the index of the column to whch new columns should be added
        ## should we use ndx+xtras[site] instead so we are adding columns
        ## "after" the ones that are alrady there??? just thinking out loud
        ndx = mut_mgr.index_from_site(site) ## hmm,adding from left side??
        xpand[ndx] = need_xtras

    return xpand

def expand_seq(xpand,seq):
    '''based on xpand dict, eqpand sequence string'''
    seq_as_list = list(seq)
    for ndx in xpand:
        seq_as_list[ndx] += "-"*xpand[ndx]
    return "".join(seq_as_list)

def expand_seqs(xpand,seqs):
    '''baed on xpand dict, expand a list/generator of seqs'''
    for s in seqs:
        s.seq = expand_seq(xpand,s.seq)
        yield s

def add_needed_dashes(seqs,mstringpairs):
    '''check if mstring pairs will require extra room in the sequences;
    and if so, add the needed dashes to the sequences'''

    first,seqs = sequtil.get_first_item(seqs,keepfirst=True)
    mut_mgr = mutant.MutationManager(first.seq)
    xpand = expansion_needed(mut_mgr,mstringpairs)
    if xpand:
        seqs = expand_seqs(xpand,seqs)
    return seqs,xpand


ALL_DASHES = re.compile(r'-')
TRAILING_DASHES = re.compile(r'-+$')

def substr_to_ssms(xlator,gseq,ndxlo=0,insertwithdashes=True):
    '''
    given a gappy substring gaseq (such as 'R---T-T-')
    and the index associated with the first character
    yield ssm's that would make up an mstring
    call: Mutation(substr_to_ssms(...)) to get Mutation object
    '''
    re_dashes = TRAILING_DASHES if insertwithdashes else ALL_DASHES
    ndxhi = ndxlo + len(gseq)
    while ndxlo < ndxhi:
        site = xlator.site_from_index(ndxlo)
        ndxlolo = xlator.index_from_site(site)
        if ndxlolo == ndxlo:
            ## substitution (or deletion) character
            subchr,gseq = gseq[0],gseq[1:]
            yield mutant.SingleSiteMutation(xlator.refseq[ndxlo],site,subchr)
            ndxlo = ndxlo + 1
        elif ndxlolo < ndxlo:
            ## insertion string
            ndx_next = xlator.index_from_site(site+1)
            insstr,gseq = gseq[:ndx_next-ndxlo],gseq[ndx_next-ndxlo:]
            insstr = re_dashes.sub('',insstr) ## remove (trailing) dashes
            if insstr:
                yield mutant.SingleSiteMutation("+",site,insstr)
            ndxlo = ndx_next
        else: # ndxlolo > ndxlo:
            raise RuntimeError(f'{ndxlolo=} > {ndxlo=}: should never happen')

def substr_to_mut(xlator,gseq,ndxlo=0,insertwithdashes=True):
    '''return mutation object based on gapped subsequence'''
    mut = mutant.Mutation(substr_to_ssms(xlator,gseq,ndxlo,insertwithdashes=insertwithdashes))
    #v.print('substr_to_mut: mut=',str(mut))
    return mut
