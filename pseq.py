'''Processed sequence'''

import warnings

import covidfast as covid
import mutant

def parse_seq_name(name):
    try:
        tokens = name.split('.',6)
        lineage = tokens[6]
        isl = tokens[5]
        isodate = tokens[4]
    except IndexError:
        raise RuntimeError(f'Bad name: {name}')
    
    return isodate,isl,lineage


class ProcessedSequence():
    def __init__(self,MM,s):
        #self.s = s
        self.seq = s.seq
        self.name = s.name
        isodate,isl,lineage = parse_seq_name(s.name)
        self.ISL = isl #covid.get_isl(s.name)
        self.lineage = lineage #covid.get_lineage_from_name(s.name)
        self.date = covid.date_fromiso(isodate)
        self.mut = MM.get_mutation(s.seq)
        self.altmut = MM.get_alt_mutation(s.seq)
        #if not all(self.altmut.values()):
        #    warnings.warn(f"init empty: {self.altmut}")

    def fits_pattern_inclusive(self,mpatt):
        '''return True if the sequence is consistent with pattern'''

        #if not all(self.altmut.values()):
        #    warnings.warn(f"incl empty: {self.altmut}")

        #if mpatt == self.mut: ## easy case: if seq and patt are identical, then True
        #    return True
        for ssm in mpatt:
            seqval = self.altmut.get(ssm.site, ssm.ref)
            #if ssm.mut == "_":
            #    ## may never happen, get's converted when read into SSM
            #    if seqval != ssm.ref:
            #        return False
            if ssm.mut == ".":
                pass
            elif ssm.mut == "*":
                if seqval[0] == ssm.ref or seqval[0] == ".":
                    return False
            elif ssm.ref == "+":
                if ssm.mut != seqval[1:]:
                    return False
            else:
                if ssm.mut != seqval[0]:
                    return False
        return True
        
    def fits_pattern(self,mpatt,alt_mpatt,exact=None,assume_inclusive=False):
        '''return True if the sequence fits the pattern'''
        #if mpatt == self.mut: ## easy case: if seq and patt are identical, then True
        #    ## easy, but not fast; quicker to skip this step and do it the "hard way"
        #    return True
        if not assume_inclusive:
            if not self.fits_pattern_inclusive(mpatt):
                return False
        if exact is None:
            exact = mpatt.exact
        if not exact:
            return True
        ## made it this far; now in exact mode
        for site,seqval in self.altmut.items():
            if not seqval:
                ## can be blank becuse altmut is a defaultdict(str)
                continue
            pval = alt_mpatt.get(site,None)
            if not pval:
                ## site in mutated seq not in mpatt
                return False
            if pval[0]=='.':
                continue
            if pval[0] == "*" and seqval != pval[1]:
                continue
            if pval != seqval:
                return False
        return True

def filter_pseqs_by_pattern(MM,mpatt,pseqs,**kw):
    '''filter an iterable of processed seqs that match the pattern'''
    ## if input mpatt is string, this converts to Mutation object
    ## also, makes sure the ssm's in mpatt are sorted
    mpatt = mutant.Mutation(sorted(mutant.Mutation(mpatt)))
    alt_mpatt = MM.get_alt_mutation_ssmlist(mpatt)
    #print("mpatt:",mpatt.as_string())
    #print("altmp:",alt_mpatt)
    for ps in pseqs:
        if ps.fits_pattern(mpatt,alt_mpatt,**kw):
            yield ps
 
