'''
definitions and subroutines for both pangocommonforms and commonformschange
'''
from collections import Counter,defaultdict
import verbose as v
import sequtil
import mutant
import covid

def mostcommonchar(clist):
    '''return the most common item in the list'''
    [(c,_)] = Counter(clist).most_common(1)
    return c

def consensus(seqlist):
    '''create a consesnsus sequence from the sequence list'''
    return "".join(mostcommonchar(clist)
                   for clist in sequtil.gen_columns_seqlist(seqlist))

def get_input_sequences(args):
    '''read input file and return firstseq and seqlist'''
    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqlist = sequtil.get_first_item(seqs,keepfirst=False)
    firstseq = first.seq

    seqlist = list(seqlist)
    v.vprint_only_summary('Invalid date:','skipped sequences')

    return firstseq,seqlist

class LineagePartition:
    '''partition sequences by lineage, provide:
    list of lineages, sorted by number of sequenes in that lineage;
    count of lineages, number of sequences for each lineage;
    individual sequence lists for each lineage;
    '''

    def __init__(self,fullseqlist,bylineage=True):
        self.sequences = defaultdict(list) ## dict of lists
        for s in fullseqlist:
            lin = covid.get_lineage_from_name(s.name) if bylineage else "N/A"
            lin = lin or "None"
            self.sequences[lin].append(s)
        self.counts = {lin: len(self.sequences[lin])
                       for lin in self.sequences}
        self.lineages = sorted(self.counts,key=self.counts.get,reverse=True)
        ## format lineage strings so they line up
        maxlinlen = max(len(lin) for lin in self.lineages+["Lineage"])
        self.fmt = "%%-%ds" % (maxlinlen,)

    def format(self,lin):
        '''return a fixed-width formatted lineage name'''
        if lin == "EMPTY" or not lin:
            lin=""
        return self.fmt % lin

def get_baseline_mutation(lin_baseline,mut_manager,linpart,protein='Spike'):
    '''from baseline lineage (eg "XBB.1.5") return the baseline mutation'''

    ##     For Other proteins, use most common form of the given baseline pango type
    if not lin_baseline:
        base_mut = mutant.Mutation("[]")
    elif protein.lower() == 'spike':
        ## for Spike, use hardcoded baseline mutation
        base_mut = mutant.Mutation(covid.get_baseline_mstring(lin_baseline))
    else:
        v.vprint('Will obtain baseline from most common',lin_baseline)
        if lin_baseline not in linpart.lineages:
            v.vprint(f'Baseline {lin_baseline} not in data!')
            v.vprint('Lineages:',linpart.lineages)
            base_mut = mutant.Mutation("[]")
        else:
            cntr = Counter(s.seq for s in linpart.sequences[lin_baseline])
            base_seq = cntr.most_common(1)[0][0]
            base_mut = mut_manager.get_mutation(base_seq)
    if lin_baseline:
        print()
        print(f"Baseline {lin_baseline}: {str(base_mut)}")

    return base_mut
