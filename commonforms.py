'''
definitions and subroutines for both pangocommonforms and commonformschange
'''
from collections import Counter,defaultdict
import verbose as v
import sequtil
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

    def __init__(self,fullseqlist,nopango=False):
        self.sequences = defaultdict(list) ## dict of lists
        for s in fullseqlist:
            lin = covid.get_lineage_from_name(s.name) if not nopango else "N/A"
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
