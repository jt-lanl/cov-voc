'''
library of routines for spike variants
'''
import re
from collections import Counter
from functools import lru_cache
import warnings

import intlist
import mutant
import colornames
try:
    from defaultspikevars import \
        get_mstrings, get_colors, get_names
except ImportError:
    warnings.warn("Cannot import defaultspikevars.py -- well, maybe we won't need it")

## VOC (variant of concern) is basically a Mutation but with also a color and a name

class VOC(mutant.Mutation):
    '''VOC = Variant Of Concern = Mutation() with color and name '''
    def __init__(self,m,c,n):
        super().__init__(m)
        self.color=c
        self.name=n

    def __hash__(self):
        return hash(tuple((self.as_string(),self.color,self.name)))

def read_colormut_fp(colormutfileptr):
    ''' yields VOC's from lines of the fileptr '''
    for line in colormutfileptr:
        line = re.sub("#.*","",line).strip()
        if not line:
            #ignore empty and commented-out lines
            continue

        ## Match: Color [Mutation]! Name, with "!" optional
        m = re.match(r'(\S+)\s+(\[.*\])(!?)\s*(\S*).*',line)
        if not m:
            warnings.warn(f"No match: {line}")
            continue
        color = m[1].strip()
        try:
            color = colornames.tohex(color)
        except KeyError:
            raise RuntimeError(f"Invalid color: {color}")

        mstring = m[2]+m[3]
        name = m[4]
        yield VOC(mstring,color,name)


class SpikeVariants():
    '''Variants of Spike protein'''
    #mostly, this is a list of VOCS
    ## Hmmm, wonder if this could be a child of the MutationManager class???

    OTHERNAME="other"
    OTHERCOLOR='#dddddd'

    def __init__(self, vocs=None, refseq=None):
        if vocs:
            self.init_from_vocs(vocs,refseq=refseq)

    @property
    def mstrings(self):
        return [v.as_string() for v in self.vocs]

    @property
    def colors(self):
        return [v.color for v in self.vocs]

    @property
    def names(self):
        return [v.name for v in self.vocs]

    @classmethod
    def default(cls,refseq=None):
        return cls().init_from_defaults(refseq=refseq)

    @classmethod
    def from_colormut(cls,colormutfile,refseq=None):
        return cls().init_from_colormut(colormutfile,refseq=refseq)

    def init_from_defaults(self,refseq=None):
        mstrings = get_mstrings()
        colors = get_colors()
        names = get_names()
        self.init_from_vocs((VOC(m,c,n)
                             for m,c,n in zip(mstrings,colors,names)),
                            refseq=refseq)
        return self

    def init_from_colormut(self,colormutfile,refseq=None):
        '''initialize from color muation table file'''
        with open(colormutfile) as fp:
            return self.init_from_fp(fp,refseq=refseq)

    def init_from_fp(self,fp,refseq=None):
        '''initialize from file pointer'''
        return self.init_from_vocs(read_colormut_fp(fp),
                                   refseq=refseq)

    def init_from_vocs(self,vocs,refseq=None):
        '''initialize from list of VOC items'''
        self.vocs = list(vocs)

        ## union of sites that appear in all the different mstrings
        sites=set()
        for voc in self.vocs:
            sites.update(ssm.site for ssm in voc)
        self.sites = list(sorted(sites))

        # nb, if no insertions then xtrachars=0 for all sites
        xtrachars = {site: 0 for site in sites}
        for voc in self.vocs:
            for ssm in voc:
                if ssm.ref == '+':
                    xtrachars[ssm.site] = max([xtrachars[ssm.site],
                                               len(ssm.mut)])
        self.xtrachars = xtrachars ## a handy data structure to keep around

        ## if refseq=None, this will still make a /plausible/ refseq
        ## if we later run set_refseq on a real refseq, then we can get a better one
        self.set_refseq(refseq)
        
        return self

    def set_refseq(self,refseq):
        '''continue initialization using refseq'''

        if not refseq:
            ## then make a pseudo-refseq, consistent with actual refseq
            ## at all sites specified in all the mutations
            ## and "x" everywhere else
            warnings.warn("I don't trust this not-specifying of refseq")

            refseq = []
            for losite,hisite in zip([0]+self.sites[:-1],self.sites):
                refseq.extend(['x']*(hisite-losite) +
                              ['-']*self.xtrachars[hisite])

            T = mutant.SiteIndexTranslator(refseq)
            for voc in self.vocs:
                for ssm in voc:
                    ndx = T.index_from_site(ssm.site)
                    if refseq[ndx] in '+x':
                        refseq[ndx] = ssm.ref
            refseq = "".join(refseq)

        self.refseq = refseq
        self.MM = mutant.MutationManager(refseq)
        self.master = self.shorten(refseq)
        return self

    @lru_cache(maxsize=None)
    def vocmatch(self,seq):
        '''return list of voc patterns that match the sequence'''
        ## ideally, length of that list is 1 or 0
        return [voc for voc in self.vocs
                if self.MM.seq_fits_pattern(voc,seq)]

    def shorten(self,seq):
        shortseq = []
        for site in self.sites:
            ndx = self.MM.index_from_site(site)
            shortseq.append(seq[ndx:ndx+1+self.xtrachars[site]])
        return "".join(shortseq)

    def flatpattern(self,voc):
        flatpatt = self.MM.pattern_from_mutation(voc)
        flatpatt = self.relpattern(flatpatt,self.refseq)
        flatpatt = self.shorten(flatpatt)
        return flatpatt

    def relpattern(self,seq,refseq,dittochar='_'):
        ## note, self not actually needed; this routine could stand alone
        ## ORRR ... we could use self's idea of what refseq is
        ## except that sometimes want to use master instead of refseq
        rseq = [dittochar if (s==r and r not in "-") else s
               for s,r in zip(seq,refseq)]
        return "".join(rseq)


    def check(self):
        Ls = len(self.ssites())
        if self.refseq:
            Lm = len(self.master)
            if Ls != Lm:
                raise RuntimeError(f"ssites:{Ls} != {Lm}=len({self.master})")
            for v in self.vocs:
                shortpatt = self.shorten(v.pattern(self.refseq))
                Lp = len(shortpatt)
                if Lp != Lm:
                    raise RuntimeError(f"{str(v)} => {shortpatt}; len={Lp} != {Lm}")

    def checkmaster(self,refseq):
        '''broken!'''  ## Have we decided what 'master' even is, yet?
        if len(self.refseq) != len(refseq):
            print(f"checkmaster: unequal lengths {len(self.refseq)} != {len(refseq)}")
            return
        for n,(m,r) in enumerate(zip(self.refseq,refseq)):
            if m!=r and m!='x':
                print(n,m,r)
                raise RuntimeError("checkmaster fail")

    def ssites(self):
        ''' list self.sites but with repeats for insertions '''
        ss = []
        for site in self.sites:
            ss.extend([site]*(1+self.xtrachars[site]))
        return ss

    def less_exact(self):
        '''
        redefine the "exact" matches to only be exact
        over the common sites in the list of vocs
        '''
        allsites = self.sites
        for voc in self.vocs:
            if not voc.exact:
                continue
            vocsites = set(ssm.site for ssm in voc)
            for site in set(allsites)-vocsites:
                ref = self.refseq[self.MM.index_from_site(site)]
                voc.append( mutant.SingleSiteMutation((ref,site,ref)) )
            voc.exact=False
            voc.sort()



    def pprint(self,**kwxtra):
        print("\n".join( intlist.write_numbers_vertically(self.ssites()) ),**kwxtra)
        print(self.master,"Master",**kwxtra)
        for v in self.vocs:
            shortpatt = self.shorten(v.pattern(self.refseq))
            #print(shortpatt,v.name,"shortpatt",**kwxtra)
            newshortpatt = self.relpattern(shortpatt,self.master)
            print(newshortpatt,v.name,str(v),**kwxtra)

    ## Q: should all of this key view stuff get moved to mkkeyfile.py ?

    def key_view1(self):
        ''' pattern variant '''

        namelen = max(len(v.name) for v in self.vocs)
        fmt = "%%s %%-%ds %%s" % namelen
        lines = intlist.write_numbers_vertically(self.ssites())
        white = "#FFFFFF"
        tic = " "
        for line in lines:
            yield fmt % (white,tic,line)
        yield fmt % (white,tic,self.master)
        yield fmt % (self.OTHERCOLOR,self.OTHERNAME,
                     "."*len(self.master))
        for v in self.vocs[::-1]:
            flatpatt = self.flatpattern(v)
            yield fmt % (v.color,v.name,flatpatt)

    def key_view2(self):
        ''' mutation string variant '''
        namelen = max(len(v.name) for v in self.vocs)
        fmt = "%%s %%-%ds %%s" % (namelen,)
        yield fmt % (self.OTHERCOLOR,self.OTHERNAME,"")
        for v in self.vocs[::-1]:
            yield fmt % (v.color,v.name,v.as_string())

    def key_view3(self,seqs):
        ''' most-common sequence m-string for each variant '''
        namelen = max(len(v.name) for v in self.vocs)
        fmt = "%%s %%-%ds %%s" % (namelen,)
        yield fmt % (self.OTHERCOLOR,self.OTHERNAME,"")
        seq_counter = dict()
        for v in self.vocs:
            seq_counter[v] = Counter()
        for s in seqs:
            for v in self.vocmatch(s.seq):
                seq_counter[v][s.seq] += 1
        for v in self.vocs[::-1]:
            z = seq_counter[v].most_common(1)
            if len(z):
                [(seq,_)] = z[:1]
                mstring = str(self.MM.get_mutation(seq))
            else:
                mstring = ""
                continue ## don't put empty mstrings into key
            yield fmt % (v.color,v.name,mstring)

    def key_view(self,view,seqs=None):
        '''get lines to be printed for keyfile'''
        assert view in [1,2,3]
        if view == 1:
            keylines = self.key_view1()
        elif view == 2:
            keylines = self.key_view2()
        elif view == 3:
            keylines = self.key_view3(seqs)
        return keylines

    def key_print(self,view,seqs=None,**kwxtra):
        '''print a keyfile '''
        for line in self.key_view(view,seqs=seqs):
            print(line,**kwxtra)

    def pyprint(self,fileptr):
        ''' writes python code to define functions:
            get_mstrings(), get_colors(), get_names()
            based on the current stuctures in SpikeVariants
        '''
        def fprint(*p,**kw):
            print(*p,file=fileptr,**kw)

        def fprint_item(name,array,quoted=True):
            fprint(f"def get_{name}():")
            fprint(f"    {name} = [")
            for item in array:
                if quoted:
                    fprint(f"        '{item}',")
                else:
                    fprint(f"        {item},")
            fprint(f"    ]")
            fprint(f"    return {name}")
            fprint()

        fprint("# Default spike variants")
        fprint_item("mstrings",self.mstrings,quoted=True)
        fprint_item("colors",self.colors,quoted=True)
        fprint_item("names",self.names,quoted=True)

        fprint("if __name__ == '__main__':")
        fprint("    import spikevariants")
        fprint("    mstrings = get_mstrings()")
        fprint("    colors = get_colors()")
        fprint("    names = get_names()")
        fprint("    vocs = [spikevariants.VOC(m,c,n)")
        fprint("            for m,c,n in zip(mstrings,colors,names)]")
        fprint("    svar = spikevariants.SpikeVariants().init_from_vocs(vocs)")
        fprint("    svar.pprint()")


if __name__ == "__main__":

    import sys

    sv = SpikeVariants.default()
    sv.set_refseq(None)
    sv.check() ## doesn't make sense without refseq!!
    sv.pprint(file=sys.stderr)

    sv = SpikeVariants.from_colormut("color-mutation-table-v4.txt")
    sv.set_refseq(None)
    sv.check()
    sv.pprint(file=sys.stderr)

    sv.key_print(2,file=sys.stderr)

    sv.pyprint(sys.stdout)
