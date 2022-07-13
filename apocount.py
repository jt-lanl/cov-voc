'''Find APOBEC patterns and mutation in DNA sequences'''

import argparse
from collections import defaultdict
import scipy.stats as stats

import verbose as v
import mutant
import sequtil
import apobec as apo

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input sequence file")
    paa("--nseq",type=int,default=0,
        help="if specified, then only analyze this many sequences")
    paa("--reversecomplement","-r",action="store_true",
        help="reverse complement the strings as they are read in")
    paa("--fwd",default="B",
        help="F = forward, R = reverse complement, B = both")
    paa("--strict",action="store_true",
        help="Use stricter APOBEC pattern definition")
    paa("--truesite",action="store_true",
        help="Use actual site number, not position in string")
    paa("--loose",action="store_false",dest='strict',
        help="Use looser APOBEC pattern definition")
    paa("--table",action="store_true",
        help="Show contingency table")
    paa("--summary",
        help="Write summary to file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args.loose = not args.strict
    return args

def get_sites(rseq,fwd='F'):
    '''
    identify sites (actually, 0-base indices) of interest in ref sequence
    output is 2-tuple of two lists:
    G_sites is set of sites with "G"
    apo_sites is subset of G_sites that have full apobec pattern
    '''
    if fwd == 'F':
        G_sites = set( apo.sites_match_char(rseq,'G') )
        G_sites -= set([len(rseq)-2,len(rseq)-1]) ## don't include last two sites
        apo_sites = set( n for n in G_sites if apo.is_apobec(rseq[n:n+3],fwd=fwd))
        return G_sites,apo_sites
    if fwd == 'R':
        C_sites = set( apo.sites_match_char(rseq,'C') )
        C_sites -= set([0,1]) ## don't include first two sites
        apo_sites = set( n for n in C_sites if apo.is_apobec(rseq[n-2:n+1],fwd=fwd))
        return C_sites,apo_sites
    if fwd == 'B':
        fwd_sites,fwd_apo = get_sites(rseq,'F')
        rvc_sites,rvc_apo = get_sites(rseq,'R')
        GC_sites = fwd_sites | rvc_sites
        apo_sites = fwd_apo | rvc_apo
        return GC_sites,apo_sites
    raise RuntimeError('fwd argument should be F, R, or B')

def count_patterns(seq,GC_sites,apo_sites,fwd='F'):
    '''count various patterns in mutations from ref rseq to seq'''
    if fwd == 'F':
        A_in_seq = set( apo.sites_match_char( seq,'A') )
        GA_sites = GC_sites & A_in_seq
        GA_apo = apo_sites & A_in_seq
        return len(GA_sites),len(GA_apo)
    if fwd == 'R':
        T_in_seq = set( apo.sites_match_char( seq, 'T') )
        CT_sites = GC_sites & T_in_seq
        CT_apo = apo_sites & T_in_seq
        return len(CT_sites),len(CT_apo)
    if fwd == 'B':
        GA_count,GA_apocount = count_patterns(seq,GC_sites,apo_sites,'F')
        CT_count,CT_apocount = count_patterns(seq,GC_sites,apo_sites,'R')
        return GA_count+CT_count, GA_apocount+CT_apocount
    raise RuntimeError('fwd argument should be F, R, or B')

def count_mutations(rseq,seq):
    '''
    count the mutations in seq relative to rseq
    discard all mutations involving dashes or Ns
    return 3-tuple: X->Y count, G->A count, C->T count
    where "X->Y" refers to any mutation except G->A or C->T 
    '''
    msites = apo.get_all_mutsites(rseq,seq)
    ctxlist = [apo.get_mut_context(n,rseq,seq) for n in msites]
    G_to_A = sum(1 for ctx in ctxlist if ctx[0]=='G')
    C_to_T = sum(1 for ctx in ctxlist if ctx[0]=='C')
    X_to_Y = len(ctxlist) - G_to_A - C_to_T
    return X_to_Y,G_to_A,C_to_T

def print_table(ctable,/,toplabels=('',''),sidelabels=('','')):
    print(f"         {toplabels[0]:>6s} {toplabels[1]:<6s}")
    print("       -----------------")
    print(f"{sidelabels[0]:>6s} | {ctable[0][0]:6d} {ctable[0][1]:6d} |"
          f" {sum(ctable[0][i] for i in [0,1]):6d}")
    print(f"{sidelabels[1]:>6s} | {ctable[1][0]:6d} {ctable[1][1]:6d} |"
          f" {sum(ctable[1][i] for i in [0,1]):6d}")
    print("       -----------------")
    print(f"        "
          f" {sum(ctable[i][0] for i in [0,1]):6d}"
          f" {sum(ctable[i][1] for i in [0,1]):6d} |"
          f" {sum(ctable[i][j] for i in [0,1] for j in [0,1]):6d}")
    
def _main(args):
    '''main'''
    v.vprint(args)

    if args.summary:
        fsummary = open(args.summary,'w')

    fwd = args.fwd.upper()[0]
    if args.loose:
        apo.loosen_apobec_rules()

    v.vprint("APOBEC_FD:",apo.APOBEC_FD)

    seqs = sequtil.read_seqfile(args.input)
    seqs = list(seqs)

    if args.reversecomplement:
        for s in seqs:
            s.seq = apo.reverse_complement(s.seq)
            
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    site_xlate = mutant.SiteIndexTranslator(first.seq)
    if args.nseq:
        seqs = seqs[:args.nseq]
    
    GC_sites, apo_sites = get_sites(first.seq,fwd=fwd)
    total_apo = len(apo_sites)
    total_GC = len(GC_sites)
    sitename = 'G' if fwd=='F' else 'C'
    v.vprint(f"{sitename} sites:",total_GC)
    v.vprint("apobecs:",total_apo)
    print("sequence-name other-mutations [contingency-table] p=p-value OR=odds-ratio")
    maxnamelen = max(len(s.name) for s in seqs)
    namefmt = "%%%ds" % maxnamelen
    for s in seqs:
        ga_count,apo_count = count_patterns(s.seq,GC_sites,apo_sites,fwd=fwd)
        xycnt,gacnt,ctcnt = count_mutations(first.seq,s.seq)
        if fwd == 'F':
            xycnt += ctcnt
        if fwd == 'R':
            xycnt += gacnt
        ctable = [[apo_count, ga_count-apo_count],
                  [total_apo - apo_count, total_GC - ga_count - total_apo + apo_count]]
        oddsratio,pvalue = stats.fisher_exact(ctable)
        ctable_line = (f'[{ctable[0][0]:3d} {ctable[1][0]:6d}'
                       f' {ctable[0][1]:3d} {ctable[1][1]:6d}]')
        print(f"{s.name:75s} {xycnt:4d} {ctable_line} p={pvalue:8.6f}, OR={oddsratio:.4g}")
        #ber = stats.barnard_exact(ctable)
        #print(f"p={ber.pvalue:8.6f} (Barnard's exact test)")
        if args.table:
            if fwd == 'F':
                sidelabels=("G->A","G->B")
            if fwd == 'R':
                sidelabels=("C->T","C->V")
            if fwd == 'B':
                sidelabels=('mutate','nonmut')
            print_table(ctable,
                        toplabels=("apobec","no-apo"),
                        sidelabels=sidelabels)
            print()
        if args.summary:
            fsummary.write(namefmt % s.name)
            msites = apo.get_all_mutsites(first.seq,s.seq)
            maposites = sorted(set(msites) & set(apo_sites))
            mxxxsites = sorted(set(msites) - set(apo_sites))
            #ctx = {n: apo.get_mut_context(n,first.seq,s.seq)
            #       for n in msites}

            muts = defaultdict(list)
            for n in maposites:
                muts['APO-'+first.seq[n]+s.seq[n]].append(n)
            for n in mxxxsites:
                muts[first.seq[n]+s.seq[n]].append(n)

            for mtype,mlist in muts.items():
                fsummary.write(f' {mtype}:')
                for n in mlist:
                    nsite = n
                    if args.truesite:
                        nsite = site_xlate.site_from_index(n)
                    fsummary.write(f' {nsite}')
            fsummary.write("\n")

    if args.summary:
        fsummary.close()
    
if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
