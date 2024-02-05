'''Compute entropy over the affected range to evaluate which tweaks are best'''

from collections import Counter
import itertools as it
from scipy import stats
import argparse
import verbose as v

from readseq import xopen
import covid
import mutant
import tweakalign as twa

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--mstrings","-m",nargs=2,action="append",
        help="pair of from/to mstrings")
    paa("--mfile","-M",action="append",
        help="read from/to mstring pairs from a file")
    paa("--output","-o",default="-",
        help="output table of entropies for various tweaks")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def sitewise_entropy(cntr,ndx):
    '''return the sum of sitewise entropies'''
    ## also return the entropy of the subsequnce as a whole
    ## The difference could be some kind of mutual information
    e_site = list()
    e_range = stats.entropy(list(cntr[k] for k in cntr))
    for n in ndx:
        cntr_n = Counter() ## count characters just at ndx=n
        for oseq,cnt in cntr.items():
            cntr_n[ oseq[n] ] += cnt
        e_site.append( stats.entropy(list(cntr_n.values())) )
    return e_range,sum(e_site)

def thesign(z):
    if z > 0: return '+'
    if z < 0: return '-'
    return ""

def _main(args):
    '''main'''
    v.vprint(args)
    
    ## Get all the mstring pair
    mstringpairs = twa.get_mstring_pairs(args)
    ## Identify where extra columns are needed (theoretical)
    xtras = twa.get_extra_chars(mstringpairs)
    v.vprint("xtras:",xtras)

    ## read the sequence file
    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,keepfirst=True)
    orig_seqlen = len(first.seq)
    mut_mgr = mutant.MutationManager(first.seq)

    ## Identify actual columns that need to be added
    xxtras = dict() ## xx is ndx based instead of site based
    for site in sorted(xtras):
        ## determine how many extra columns there already are
        mm_xtras = len(mut_mgr.indices_from_site(site))-1
        need_xtras = xtras[site] - mm_xtras
        if need_xtras > 0:
            ## we'll need to add some columns
            v.vprint("xtras:",site,mm_xtras,"->",xtras[site])
            ## ndx is the index of the column to whch new columns should be added
            ## should we use ndx+xtras[site] instead so we are adding columns
            ## "after" the ones that are alrady there??? just thinking out loud
            ndx = mut_mgr.index_from_site(site) ## hmm,adding from left side??
            xxtras[ndx] = need_xtras

    if len(xxtras):
        seqs = twa.add_xtra_dashes(seqs,xxtras)
    
    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    mut_mgr = mutant.MutationManager(first.seq)
    v.vprint("Sequence length:",orig_seqlen,"->",len(first.seq))
    ## Maybe only need to make a list if there is more than one mstringpair?
    if len(mstringpairs)>1:
        v.vvprint("reading sequences...",end="")
        seqs = list(seqs)
        v.vvprint("ok. Read",len(seqs),"sequences")
    
    ## Now we have our base sequences; lets go make some mutations!

    maxwidth = max(len(f'{ma}=>{mb}') for ma,mb in mstringpairs)
    fmt = "%%%ds" % maxwidth
    output=list()
    for ma,mb in mstringpairs:
        ndxlo,ndxhi,seq_a,seq_b = twa.mstrings_to_ndx_seqs(mut_mgr,ma,mb)
        v.vvprint(f"{ma} -> {mb}:   {seq_a} -> {seq_b}")

        ## only need the indices for which seq_a[n] != seq_b[n]
        ndiff = [n for n in range(len(seq_a))
                 if seq_a[n] != seq_b[n]]
        seq_an = "".join(seq_a[n] for n in ndiff)
        seq_bn = "".join(seq_b[n] for n in ndiff)

        v.vprint(f'ndiff={ndiff}: {seq_an} -> {seq_bn}')        
        
        ## oseqs are original seqs, confinded to the relevant ndx's
        oseqs = [s.seq[ndxlo:ndxhi] for s in seqs]
        ca = sum(1 for seq in oseqs if seq == seq_a)
        cb = sum(1 for seq in oseqs if seq == seq_b)
        #oseqs = ["".join(seq[n] for n in ndiff)
        #         for seq in oseqs]
        c_oseqs = Counter(oseqs)
        v.vprint('Distinct subseqs:',len(c_oseqs))

        ## count appearances of seq_a and seq_b
        cax = c_oseqs[seq_a]
        cbx = c_oseqs[seq_b]
        v.vvprint(f'old: ca={cax}, cb={cbx}; vs New: ca={ca}, cb={cb}')

        ## Set up the c_aseqs and c_bseqs Counters
        ## to characterize (using entropy, in particular)
        ## how the sequenecs differ under the a-version
        ## vs the b-version of the alignment.
        c_aseqs = Counter()
        c_bseqs = Counter()
        for seq,cnt in c_oseqs.items():
            aseq = seq_a if seq==seq_b else seq
            c_aseqs[aseq] += cnt
            bseq = seq_b if seq==seq_a else seq
            c_bseqs[bseq] += cnt

        v.vprint('o:',len(oseqs), 'c:',len(c_oseqs),
                 'a:',len(c_aseqs),'b:',len(c_bseqs))

        Eor,Eos = sitewise_entropy(c_oseqs,ndiff)
        Ear,Eas = sitewise_entropy(c_aseqs,ndiff)
        Ebr,Ebs = sitewise_entropy(c_bseqs,ndiff)

        v.vprint('O:',Eor,Eos)
        v.vprint('A:',Ear,Eas)
        v.vprint('B:',Ebr,Ebs)

        if args.output.endswith('.tsv'):
            output.append("\t".join(map(str,(ma,mb,ca,cb,Ebs-Eas,
                                             thesign(Ebs-Eas)))))
        else:
            output.append(" ".join([fmt % f'{ma}->{mb}',
                                    "%7d %7d %+12.9f" % (ca,cb,Ebs-Eas)]))

                        
    with xopen(args.output,"w") as fout:
        if args.output.endswith('.tsv'):
            tsv_header = '\t'.join(['mstring-a',
                                    'mstring-b',
                                    'seqcount-a',
                                    'seqcount-b',
                                    'entropy diffrence (b-a)',
                                    'sign of difference'])
            print(tsv_header,file=fout)
        for line in output:
            print(line,file=fout)
        

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
