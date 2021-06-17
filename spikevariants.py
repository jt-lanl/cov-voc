import re
from collections import namedtuple
import warnings

import intlist
import mutant
import colornames

### Note: now that we've made a distinction between pattern and re_pattern, we can use _'s in pattern instead of letters
### which will be much more readable

## Defaults:

def get_sites():
    sites = [
           5,   13,   18,   20,   26,   52,   67,   69,   70,   75,   76,   80,   95,
          98,  138,  141,  142,  143,  144,  152,  153,  154,  157,  180,  184,  189,
         190,  215,  222,  242,  243,  244,  246,  247,  248,  249,  250,  251,  252,
         253,  262,  272,  367,  417,  439,  452,  477,  478,  484,  490,  494,  501,
         570,  613,  614,  653,  655,  675,  677,  681,  701,  716,  732,  769,  772,
         796,  859,  888,  950,  982, 1027, 1071, 1118, 1176, 1219,
    ]
    return sites

def get_master():
    master = \
        'LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQDAHQQPATTGVDTFDSTQDVG'
    return master

def get_mutants_colors_names():
    mcnlist = [
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQDAHQQPATTGVDTFDSTQDVG', '#eeeeee', 'Ancestral'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#DCDCDC', 'G=D614G'),
        ('FSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#DCDCDC', 'L5F_G'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDVLALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#D3D3D3', 'GV=A222V'),
        ('LSFTPQAHVGTDTSDLGVYWMEFEGLRDVLALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#D3D3D3', 'L18F_GV'),
        ('FSFTPQAHVGTDTSDLGVYWMEFEGLRDVLALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#D3D3D3', 'L5F_GV'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLNTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#B0C4DE', 'S477N'),
        ('LSLTPQAHVGTDTFDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#E6E6FA', 'S98F'),
        ('LSLTPQAHVGTDTSDLGVYWTEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#E0FFFF', 'M153T'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFPNAQGAHQQPATTGVDTFDSTQDVG', '#CD5C5C', 'S494P'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTKFSNAQGAHQQPATTGVDTFDSTQDVG', '#FF7F50', 'E484K'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNRSTEFSNAQGAHQQPATTGVDTFDSTQDVG', '#F08080', 'L452R'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQHATTGVDTFDSTQDVG', '#FF00FF', 'P681H'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQHPATTGVDTFDSTQDVG', '#FF00FF', 'Q677H'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHHQPATTGVDTFDSTQDVG', '#FF00FF', 'Q675H'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQQRATTGVDTFDSTQDVG', '#FF00FF', 'P681R'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQRPATTGVDTFDSTQDVG', '#FF00FF', 'Q677R'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHRQPATTGVDTFDSTQDVG', '#FF00FF', 'Q675R'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDVLALRSYLTPGDAPVKNLSTEFSNAQGAHHQPATTGVDTFDSTQDVG', '#FF00FF', 'Q675H_GV'),
        ('FSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQHPATTGVDTFDSTQDVG', '#FF00FF', 'Q677H'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSNAQGAHQPPATTGVDTFDSTQDVG', '#9400D3', 'Q677P'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSTAQGAHQQPATTGVDTFDSTQDVG', '#228B22', 'N501T'),
        ('LSLTPQAHVGTDTSDLGVYWMEFEGLRDALALRSYLTPGDAPVKNLSTEFSYAQGAHQQPATTGVDTFDSTQDVG', '#006400', 'N501Y'),
        ('.......--.........-................................YD.G....H.I.......A..H..', '#FFA500', 'B.1.1.7'),
        ('.I.................C.........................R........G....................', '#00008B', 'B.1.429/7'),
        ('...........A...............G.---...........N....K..Y..G.....V..............', '#DDA0DD', 'B.1.351'),
        ('..FNS.........Y...........S................T....K..Y..G.Y.............I..F.', '#B22222', 'P.1'),
        ('...........................................K....K..N..G..................F.', '#FF0000', 'P.2'),
        ('F...........I..........................G........K.....G.....V..............', '#800080', 'B.1.526'),
        ('......................L...................F..........HD....R...............', '#D2691E', 'A23.1'),
        ('.....RV--.........-.............................K.....G...H........L.......', '#8FBC8F', 'B.1.525'),
        ('...............................................K......G....H..A............', '#4169E1', 'B.1.1.519'),
        ('................S......V..............................G...H................', '#5F9EA0', 'B.1.234'),
        ('...........G......-...S......................R........G...........N.H......', '#7CFC00', 'B.1.526.1'),
        ('..F..........................................R.....Y..DVY........Y........V', '#7FFFD4', 'A.27'),
        ('....................T...S.............................G...H................', '#00FFFF', 'B.1.1.284'),
        ('...................L............................K.....G........V...........', '#8A2BE2', 'R.1'),
        ('.........VI.....................-------N.....Q...S....G...........N........', '#008000', 'B.1.1.1'),
        ('............I...D....K.......................R..Q.....G....R...........H...', '#8B008B', 'B.1.617.1'),
        ('..................................................PY..G....H.I.......S..D..', '#556B2F', 'B.1.575'),
        ('.......--................F..................K.........G.........I..........', '#FFD700', 'B.1.258.17'),
        ('............................V...........SL............G....................', '#4B0082', 'B.1.177'),
        ('...................R........................K.........G....R...............', '#FF69B4', 'B.1.466.2'),
        ('.S.............---...........................R........G...........N.D......', '#FFC0CB', 'A.2.5.2'),
        ('other', '#000000', 'other'),
    ]
    return mcnlist

def get_mutants():
    return [mcn[0] for mcn in get_mutants_colors_names()]

def get_colors():
    return [mcn[1] for mcn in get_mutants_colors_names()]

def get_names():
    return [mcn[2] for mcn in get_mutants_colors_names()]

def get_regex_pattern(master,pattern):
    ## warning, vaguely redundant with mutant.Mutation.regex_pattern()
    ## bigger warning, returns a COMPILED regex, not just a string
    r=""
    for a,b in zip(master,pattern):
        if b == "!":
            r += "[^"+a+"]"
        else:
            r += b
    return re.compile(r)

## VOC (variant of concern)
VOC = namedtuple('VOC',['pattern','color','name','re_pattern'])

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
        self.vocs = [VOC(p,c,n,get_regex_pattern(p))
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
                    color = colornames.random_hex()

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

        vocs = [ VOC(pattern,color,name,get_regex_pattern(master,pattern))
                 for pattern,color,name in zip(mutant_patterns,colors,names) ]


        #vocs.insert(0, VOC(master,"#eeeeee","Ancestral",get_regex_pattern(master,master)))

        self.sites = sites
        self.master = master
        self.vocs = vocs
        
        return self

    def append_other(self):
        ## really don't like this function
        self.vocs.append( VOC("other", "#000000", "other", re.compile("$-") ) )
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
        for m,c,n in zip(self.mutants,self.colors,self.names):
            name = n if n else ""
            print(fmt % m,c,name,**kwxtra)        

    def pyprint(self,fileptr):
        ''' writes python code to define functions:
            get_sites() get_mutants(), get_colors(), get_names(), get_master()
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

        fprint_item("sites",intlist.format_intlist(self.sites),quoted=False)
        #fprint_item("mutants",self.mutants,quoted=True)
        #fprint_item("colors",self.colors,quoted=True)
        #fprint_item("names",self.names,quoted=True)
        fprint("def get_master():")
        fprint("    master = \\")
        fprint(f"        '{self.master}'")
        fprint("    return master")
        fprint()           

        fprint("def get_mutants_colors_names():")
        fprint("    mcnlist = [")
        for m,c,n in zip(self.mutants,self.colors,self.names):
            fprint(f"        ('{m}', '{c}', '{n}'),")
        fprint("    ]")
        fprint("    return mcnlist")
        fprint()
        fprint("def get_mutants():")
        fprint("    return [mcn[0] for mcn in get_mutants_colors_names()]")
        fprint()
        fprint("def get_colors():")
        fprint("    return [mcn[1] for mcn in get_mutants_colors_names()]")
        fprint()
        fprint("def get_names():")
        fprint("    return [mcn[2] for mcn in get_mutants_colors_names()]")

if __name__ == "__main__":

    import sys
    import readseq

    sv = SpikeVariants().init_from_defaults()
    sv.check()
    sv.pprint(file=sys.stderr)

    seqs = readseq.read_seqfile("Data/wuhan.fasta")
    sv = SpikeVariants().init_from_colormut("color-mut-Apr18.txt",seqs[0].seq,includeother=False).append_other()
    sv.check()
    sv.pprint(file=sys.stderr)
    
    sv.pyprint(sys.stdout)

        
