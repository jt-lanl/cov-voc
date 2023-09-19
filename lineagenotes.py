'''read the lineage_notes.txt file'''
## Available at:
## https://github.com/cov-lineages/pango-designation/blob/master/lineage_notes.txt
import re
import verbose as v

def read_lineage_file(lineagefile):
    '''return a list of lineages, and dict of aliases'''
    ## parse lines of the form:
    ## GL.1	Alias of XAY.1.1.1.1, Europe, S:D420N, C19441T, from issue #2032
    lineages = list()
    fullname = dict()
    with open(lineagefile,'r') as fin:
        ## skip top two lines
        fin.readline()
        fin.readline()
        ## begin reading in earnest
        for line in fin:
            line = line.strip()
            if not line:
                continue
            tokens = line.split()
            name = tokens[0]
            lineages.append(name)
            if tokens[1] == 'Alias':
                full = tokens[3]
                if full[-1] == ',':
                    ## strip trailing comma
                    full = full[:-1]
                fullname[name] = full
            elif re.match(r'[A-Z][A-Z]\..*',tokens[0]):
                v.vprint('missing "Alias of"?',line)
    return lineages,fullname

class LineageNotes:
    '''parse the lineage_notes.txt file'''
    def __init__(self,lineagefile):
        self.lineages, self.fullname = read_lineage_file(lineagefile)

    def report_size(self):
        '''return string with size of lineages and aliases'''
        sizelines = []
        sizelines.append( f'lineages: {len(self.lineages)}' )
        sizelines.append( f'aliases : {len(self.fullname)}' )
        return "\n".join(sizelines)

    def inconsistencies(self,remove=True):
        '''check consistency of lineage notes file'''
        inconsistent=[]
        bad_aliases = list()
        for alias,full in self.fullname.items():
            _,numbers = alias.split('.',1)
            if not full.endswith(numbers):
                inconsistent.append(f'Inconsistent ending: {alias}->{full}')
                bad_aliases.append(alias)
        if remove:
            ## remove inconsistent aliases/lineages
            for alias in bad_aliases:
                del self.fullname[alias]
            self.lineages = [lin for lin in self.lineages
                             if lin not in bad_aliases]

        return inconsistent
