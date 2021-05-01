import re

import intlist
import mutant
import colornames

## Defaults:

def get_sites():
    return [13,20,26,52,67,69,70,80,95,138,144,152,153,157,215,222,242,
            243,244,253,367,417,439,452,453,477,478,484,501,570,613,614,
            655,675,677,681,701,716,732,888,982,1118,1176]

def get_master():
    master = \
        'STPQAHVDTDYWMFDALALDVKNLYSTENAQDHQQPATTFSDV'
    return master

def get_mutants():
    mutants = [
        'S.PQA.....YWMF.ALALDVKNLYSTENAQGH.QPAT.FSD.', ##  G clade, 
        'STPQA...TDYWMF.VLALDVKNLYSTENAQGH...ATTFSDV', ##  GV clade, 
        'STPQAHVDTDYWMFDALALDVKNLYNTENAQGHQQPATTFSDV', ##  S477N, 
        'STPQA..DTDYWMFDALALDVKNLFSTENAQGHQQPATTFSDV', ##  Y453F , Denmark
        'STPQA..DTDYWMFDALALDVKKLYSTENAQGHQQPATTFSDV', ##  N439K, 
        'STPQAHVDTDYWTFDALALDVKNLYSTENAQGH...ATTFSD.', ##  M153T, Japan
        'STPQAHVDTDYWMFDALALDVKNLY.TENAQGHQHPATTFSDV', ##  Q677H, 
        'STPQAHVDTDYWMFDALALDVKNLYSTENAQGHQQHA.TFSDV', ##  P681HR+716-any(I), 
        'STPQAHVDTDYWMFDALALDVKNLYSTENAQGHQQRA.TFSDV', ##  P681HR, 
        'STPQAHVDTDYWMFDALALDVKNLYSTENAQGHQPPATTFSDV', ##  Q677P, 
        'STPQA..DTDY.MFDALALDVKNLYSTKNAQGH...ATTFSD.', ##  E484K, 
        'STPQA..DTDYWMFDALALDVKNLYSTETAQGH...ATTFSDV', ##  N501T, 
        'STPQA..DTDYWMFDALALDVKNLYSTEYAQGH...ATTFSDV', ##  N501Y, 
        'STPQA..DTD.WMFDALALDVKNRYSTENAQGH...ATTFSDV', ##  L452R, 
        'STPQA--.T.-WMFDALALDVKNLYSTEYDQGHQ.HAITFAHV', ##  B.1.1.7 UK, England
        'STPQAHVATDYWMFGA---DVNNLYSTKYAQGHQQPVTTFSDV', ##  B.1.352 ZA, South Africa
        'IT.QAHVDTDYCMFDALALDVKNRYSTENAQGHQQPATTFSDV', ##  B.1.429 US CA, California
        'SNSQAHVDTYYWMFDALALDVTNLYSTKYAQGYQQPATTFSDF', ##  B.1.1.248 BR, Brazil
        'ST.QAHVDTDYWMLDALALDFKNLYSTENAHDHQQRATTFSDV', ##  A.23.1 UG, Uganda
        'STPQAHVDIDYWMFDALALGVKNLY.T.NAQGHQQP.TTFSDV', ##  B.1.526 US NY, New York
        #   * was R below where that asterisk is
        'STP.V--DTD-WMFDALALDVKNLYSTKNAQGHQHPATTLSDV', ##  A67V, Nigeria
        'STPQAHVDTDYWMFDALALDVKNLYSKENAQGHQQHATAFSDV', ##  T478K, Mexico
        'STPQAHVGTD-WMSDALALDVKNRYSTENAQGHQQPATTFSDV', ##  Alt-NY, Alt-L452R
        'other',
        ]
    return mutants

def get_colors():
    colors = [
        '#e2e2e2', #lighter gray
        '#cccccc', #light gray
        '#07a0b0', ## darker cyan?
        '#ffcc00', #'gold',
        '#f3e600', #lightercream',
        '#70ff30', #lightgreen',
        '#9900ff', #purple',
        '#ff33ff', #magenta
        '#ff33ff', #magenta
        '#1f77b4', #blue 
        '#ee0000', #red0
        '#00b300', #darkgreen
        '#009000', #darkergreen
        '#5f57b4', #blue 
        #'#b38600', #darkercream',
        '#ff7f0e', #orange
        '#e5ccff', #light-violet 'lavender',
        '#0000b4', #darker blue 
        '#aa0000', #darkred
        '#ffbf00', #lightbrown',
        '#aa00ff', #purple',  ## again!
        '#70db70', #lightgreen',
        '#5f5fff', #babyblue
        '#60ef20', #lightgreen',
        #'#27cedf', #light-cyan ...        'teal',
        #'#e5ccff', #light-violet 'lavender',
        #'#cc9900', #medbrown',
        #'#00b300', #darkgreen
        #'#e6e600', #yellow',
        '#000000', #black
        ]
    return colors

def get_names():
    return [m for m in get_mutants()]

class SpikeVariants():
    def __init__(self, sites=None, master=None, mutants=None, colors=None, names=None):
        self.sites = sites
        self.master = master
        self.mutants = mutants
        self.colors = colors
        self.names = names

    def init_from_defaults(self):
        self.sites = get_sites()
        self.master = get_master()
        self.mutants = get_mutants()
        self.colors = get_colors()
        self.names = get_names()
        return self

    def set_sites(self,sites):
        self.sites = sites
    def set_master(self,master):
        self.master = master
    def set_mutants(self,mutants):
        self.mutants = mutants
    def set_colors(self,colors):
        self.colors = colors

    def check(self):
        L = len(self.sites)
        assert( len(self.master) == L )
        for m in self.mutants:
            if m != "other":
                assert( len(m) == L )
        assert( len(self.mutants) == len(self.colors) )
        assert( len(self.names) == len(self.colors) )

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

        def fprint_item(name,array,quoted=False):
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

        fprint_item("sites",intlist.format_intlist(self.sites)) #,intro="sites = [")self.sites)
        fprint_item("mutants",self.mutants,quoted=True)
        fprint_item("colors",self.colors)
        fprint_item("names",self.names,quoted=True)
        fprint("def get_master():")
        fprint("    master = \\")
        fprint(f"        '{self.master}'")
        fprint("    return master")




## sv_fromfile takes a filename and reference sequence, and returns
## sitelist array, master string, mutants array of strings, and colors array of hex strings
## which are stored in a SpikeVariants class structure

def sv_fromfile(colormutfile,refseq,includeother=True):

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

            ## Match: Color [Mutation]! Name, with "!" optional and Name optional
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
        sites.update(m.site for m in mut)
    sites = list(sorted(sites))
    print("sites =",sites)

    master = "".join(refseq[n-1] for n in sites)
    mutant_patterns = [master]
    for mut,xact in zip(mutants,exact):
        full_pattern = mut.pattern(refseq,exact=xact)
        pattern = ''.join(full_pattern[n-1] for n in sites)
        mutant_patterns.append( pattern )
        
    colors = ["#eeeeee"] + colors
    names = ["Ancestral"] + names

    spv = SpikeVariants(sites = sites,
                        master = master,
                        mutants = mutant_patterns,
                        colors = colors,
                        names = names)

    if includeother:
        spv.mutants.append("other")
        spv.colors.append("#000000")
        spv.names.append("other")
        
    return spv


if __name__ == "__main__":

    import sys

    sv = SpikeVariants().init_from_defaults()
    sv.check()
    sv.pprint(file=sys.stderr)

    sv.pyprint(sys.stdout)

        
