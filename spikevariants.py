import re
from collections import namedtuple
import warnings

import intlist
import mutant
import colornames
from defaultspikevars import \
    get_sites, get_master, get_mutants, get_colors, get_names


def get_regex_pattern(master,pattern):
    ## warning, vaguely redundant with mutant.Mutation.regex_pattern()
    ## bigger warning, returns a COMPILED regex, not just a string
    r=""
    for a,b in zip(master,pattern):
        if b == "*":
            r += "[^"+a+"]"
        else:
            r += b
    return re.compile(r)

## VOC (variant of concern)
#VOC = namedtuple('VOC',['pattern','color','name','re_pattern'])

class VOC():
    def __init__(self,p,c,n):
        self.pattern=p
        self.color=c
        self.name=n
        self.re_pattern=None

    def get_re_pattern(self,master):
        self.re_pattern=get_regex_pattern(master,self.pattern)
        return self.re_pattern

    def get_best_mstring(self,sites,master):
        # "best" == shortest
        ex_string = self.get_mstring(sites,master,exact=True)
        un_string = self.get_mstring(sites,master,exact=False)
        if len(un_string) <= len(ex_string):
            return un_string
        else:
            return ex_string
        

    def get_mstring(self,sites,master,exact=False):
        muts = []
        for n,m,p in zip(sites,master,self.pattern):
            if p == '.' and not exact:
                continue
            if p == m and exact:
                continue
            muts.append( "%s%d%s" % (m,n,p) )
        mstring = "[" + ",".join(muts) + "]"
        if exact:
            mstring += "!"
        return mstring

    def relpattern(self,master):
        s = ""
        for p,m in zip(self.pattern,master):
            s += p if p != m else "_"
        return s
            
            
        

class SpikeVariants():
    def __init__(self, sites=None, master=None, vocs=None):
        self.sites = sites
        self.master = master
        self.vocs = vocs

    @property
    def mutants(self):
        return [v.pattern for v in self.vocs]

    @property
    def colors(self):
        return [v.color for v in self.vocs]

    @property
    def names(self):
        return [v.name for v in self.vocs]

    def init_from_defaults(self):
        self.sites = get_sites()
        self.master = get_master()
        patterns = get_mutants()
        colors = get_colors()
        names = get_names()
        self.vocs = [VOC(p,c,n)
                     for p,c,n in zip(patterns,colors,names)]
        for v in self.vocs:
            v.re_pattern = get_regex_pattern(self.master,v.pattern)
        return self

        
    def init_from_colormut(self,colormutfile,refseq=None):

        mutants=[]
        colors=[]
        exact=[]
        names=[]

        with open(colormutfile) as f:
            for line in f:
                line = re.sub("#.*","",line).strip()
                if not line:
                    #ignore empty and commented-out lines
                    continue

                ## Match: Color [Mutation]! Name, with "!" optional 
                m = re.match("(\S+)\s+(\[.*\])(!?)\s*(\S*).*",line)
                if not m:
                    warnings.warn(f"No match: {line}")
                color = m[1].strip()
                try:
                    color = colornames.tohex(color)
                except KeyError:
                    raise RuntimeError(f"Invalid color: {color}")

                mutants.append( mutant.Mutation(m[2]) )
                colors.append(color)
                exact.append(bool(m[3]))
                names.append(m[4])

        ## Get list of all sites
        ## union of sites that appear in all the different mutants
        sites=set()
        for mut in mutants:
            sites.update(ssm.site for ssm in mut)
        sites = list(sorted(sites))

        if not refseq:
            ## then make a pseudo-refseq, consistent with actual refseq
            ## at all sites specified in all the mutations
            ## and "x" everywhere else
            refval = dict()
            for mut in mutants:
                for ssm in mut:
                    refval[ssm.site] = ssm.ref
            refseq = ["x"] * (max(sites)+1)
            for n in sites:
                refseq[n-1] = refval[n]
            refseq = "".join(refseq)

        mutant_patterns = []
        for mut,xact in zip(mutants,exact):
            if not mut.checkref(refseq,verbose=True):
                warnings.warn(f"Mismatch with refseq in mutant: {mut}")
            full_pattern = mut.pattern(refseq,exact=xact)
            pattern = ''.join(full_pattern[n-1] for n in sites)
            mutant_patterns.append( pattern )

                
        master = "".join(refseq[n-1] for n in sites)

        vocs = [ VOC(p,c,n)
                 for p,c,n in zip(mutant_patterns,colors,names) ]
        for v in vocs:
            v.get_re_pattern(master)

        #vocs.insert(0, VOC(master,"#eeeeee","Ancestral",get_regex_pattern(master,master)))

        self.sites = sites
        self.master = master
        self.vocs = vocs
        
        return self

    def append_other(self,master):
        ## really don't like this function
        vother = VOC("." * len(master), "#000000", "other")
        vother.get_re_pattern(self.master)
        self.vocs.append( vother )
        return self

    def check(self):
        L = len(self.sites)
        assert( len(self.master) == L )
        for m in self.mutants:
            if m != "other":
                assert( len(m) == L )
        assert( len(self.mutants) == len(self.colors) )
        assert( len(self.names) == len(self.colors) )

    def checkmaster(self,refseq):
        for n,s in enumerate(self.sites):
            assert ( self.master[n] == refseq[s-1] )
            

    def pprint(self,**kwxtra):
        print("\n".join( intlist.write_numbers_vertically(self.sites) ),**kwxtra)
        print(self.master,"Master",**kwxtra)
        fmt="%%-%ds" % (len(self.master),)
        for v in self.vocs:
            p,c,n = v.pattern, v.color, v.name
            r = v.relpattern(self.master)
            name = n if n else ""
            print(fmt % r,c,name,**kwxtra)

    def key_print(self,**kwxtra):
        ''' pattern variant '''
        namelen = max(len(v.name) for v in self.vocs)
        fmt = "%%-%ds" % namelen
        lines = intlist.write_numbers_vertically(self.sites)
        white = "#FFFFFF"
        for line in lines:
            print(white, fmt % ("'",), line, **kwxtra)
        print(white, fmt % ("'",), self.master, **kwxtra)
        for v in self.vocs[::-1]:
            print(v.color,fmt % (v.name,), v.relpattern(self.master), **kwxtra)

    def key_print2(self,**kwxtra):
        ''' mutation string variant '''
        mstringlist = [v.get_best_mstring(self.sites,self.master)
                       for v in self.vocs]
        strlen = max(len(mstr) for mstr in mstringlist)
        namelen = max(len(v.name) for v in self.vocs)
        fmt = "%%-%ds %%-%ds" % (namelen,strlen)
        for v in self.vocs[::-1]:
            m = v.get_best_mstring(self.sites,self.master)
            print(v.color,fmt % (v.name,m),**kwxtra)

    def pyprint(self,fileptr):
        ''' writes python code to define functions:
            get_sites() get_master(), get_mutants(), get_colors(), get_names()
        based on the current stuctures in SpikeVariants'''
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
        fprint()
        fprint_item("sites",intlist.format_intlist(self.sites),quoted=False)
        fprint("def get_master():")
        fprint("    master = \\")
        fprint(f"        '{self.master}'")
        fprint("    return master")
        fprint()           
        fprint_item("mutants",self.mutants,quoted=True)
        fprint_item("colors",self.colors,quoted=True)
        fprint_item("names",self.names,quoted=True)

if __name__ == "__main__":

    import sys
    import readseq

    sv = SpikeVariants().init_from_defaults()
    sv.check()
    sv.pprint(file=sys.stderr)

    seqs = readseq.read_seqfile("Data/wuhan.fasta.gz")
    sv = SpikeVariants().init_from_colormut("color-mutation-table-v2.txt")
    sv.check()
    sv.pprint(file=sys.stderr)

    sv.key_print2(file=sys.stderr)
    
    sv.pyprint(sys.stdout)

        
