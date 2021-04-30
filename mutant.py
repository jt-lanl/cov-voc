import re

class SingleSiteMutation():
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

def parse_mutline(line):
    '''Convert string of form '[W152R,N439K,D614G,P681R]' 
    into a list of ssm's (single-site mutations)'''
    line = line.strip()
    muts = line.strip("[]").split(",")
    ssmlist = [SingleSiteMutation(m.strip()) for m in muts]
    return ssmlist
    
def rd_mutfile(filename):
    mutants=[]
    with open(filename) as f:
        for line in f:
            line = re.sub("#.*","",line)
            line = line.strip()
            if not line:
                continue
            ssms = parse_mutline(line)
            mutants.append(ssms)
    return mutants

def name_from_ssms(ssms):
    return "[" + ",".join(str(ssm) for ssm in ssms) + "]"

def mkmutname(refseq,seq):
    ssmlist = []
    for n,(r,c) in enumerate(zip(refseq,seq)):
        if c != r:
            ssmlist.append( f"{r}{n+1}{c}" )
    return "[" + ",".join(ssmlist) + "]"
        

