'''
definitions and subroutines for both pangocommonforms and commonformschange
'''

import sys
from collections import Counter,defaultdict
import verbose as v
import sequtil
import mutant
import covid

def commonforms_args(argparser):
    '''args for arparse common to pangocommonforms and commonformschange'''
    covid.corona_args(argparser)
    paa = argparser.add_argument_group("Common Forms Options").add_argument
    paa("--npatterns","-n",type=int,default=0,
        help="How many of the most common patterns per lineage (0=all)")
    paa("--mincount","-m",type=int,default=10,
        help="Show only patterns with at least this many counts")
    paa("--protein",default="Spike",
        help="Protein name to be used in the header")
    paa("--baseline",default="XBB.1.5",type=str.upper,
        choices=tuple(map(str.upper,covid.BASELINE_MSTRINGS)),
        help="Use this sequence as basline for mutation strings")
    paa("--lineagebaseline",action="store_true",
        help="Use each lineage most common form as mstring baseline for that lineage")
    paa("--bylineage",action="store_true",
        help="Partition sequences by pango lineage")
    paa("--notbylineage",action="store_false",dest='bylineage',
        help="Do not partition sequences by pango lineges")

def commonforms_fixargs(args):
    '''after args are parsed, do some checks and fixes'''
    if args.lineagebaseline or args.baseline == 'WUHAN':
        args.baseline = None
    if not args.bylineage and args.lineagebaseline:
        v.print('Warning: use --bylineage if you also want --lineagebaseline.')
    args = covid.corona_fixargs(args)
    return args

def mostcommonchar(clist):
    '''return the most common item in the list'''
    [(item,_)] = Counter(clist).most_common(1)
    return item

def consensus(seqlist):
    '''create a consesnsus sequence from the sequence list'''
    return "".join(mostcommonchar(clist)
                   for clist in sequtil.gen_columns_seqlist(seqlist))

def get_input_sequences(args,minseqs=1):
    '''read input file and return firstseq and seqlist'''
    seqs = covid.read_filter_seqfile(args)
    seqs = sequtil.checkseqlengths(seqs)

    first,seqlist = sequtil.get_first_item(seqs,keepfirst=False)
    firstseq = first.seq

    seqlist = list(seqlist)
    v.vprint_only_summary('Invalid date:','skipped sequences')

    if len(seqlist) < minseqs:
        v.print(args)
        v.print(f'Only {len(seqlist)} sequences -- aborting.')
        sys.exit(1) ## avoid 'RuntimeError' because of its longwinded traceback

    return firstseq,seqlist

class LineagePartition:
    '''partition sequences by lineage, provide:
    list of lineages, sorted by number of sequenes in that lineage;
    count of lineages, number of sequences for each lineage;
    individual sequence lists for each lineage;
    '''

    def __init__(self,fullseqlist,bylineage=True,restrict_to=None):
        self.sequences = defaultdict(list) ## dict of lists
        for s in fullseqlist:
            lin = covid.get_lineage(s) if bylineage else "N/A"
            lin = lin or "None"
            if restrict_to and lin not in restrict_to:
                continue
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
        base_mut = mutant.Mutation()
    elif protein.lower() == 'spike':
        ## for Spike, use hardcoded baseline mutation
        base_mut = mutant.Mutation.from_mstring(covid.get_baseline_mstring(lin_baseline))
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
