'''Compute entropy over the affected range to evaluate which tweaks are best'''

from collections import Counter
import argparse
from scipy import stats
import verbose as v

from readseq import xopen
import covid
import mutant
import tweak as tku
from tweak import IndexTweak

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--mstrings","-m",nargs=2,action="append",
        help="pair of from/to mstrings")
    paa("--mfile","-M",action="append",
        help="read from/to mstring pairs from a file")
    paa("--tfile","-T",action="append",
        help="read index-based tweaks from a file")
    paa("--output","-o",default="-",
        help="output table of entropies for various tweaks")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def sitewise_entropy(cntr,ndxlist):
    '''return the sum of sitewise entropies'''
    ## also return the entropy of the subsequnce as a whole
    ## The difference could be some kind of mutual information
    e_site = list()
    e_range = stats.entropy(list(cntr[k] for k in cntr))
    for ndx in ndxlist:
        cntr_n = Counter() ## count characters just at ndx
        for oseq,cnt in cntr.items():
            cntr_n[ oseq[ndx] ] += cnt
        e_site.append( stats.entropy(list(cntr_n.values())) )
    return e_range,sum(e_site)

def thesign(val):
    '''return plus or minus sign as string, depending on sign of val'''
    return '+' if val>0 else ('-' if val<0 else '')

def tweak_counter(cntr,sa,sb):
    '''modify counter by replacing sb->sa'''
    tcntr = cntr.copy()
    tcntr[sa] += cntr[sb]
    tcntr[sb] = 0
    return tcntr


def _main(args):
    '''main'''
    v.vprint(args)

    ## Read sequence file
    seqs = covid.read_filter_seqfile(args)

    ## Get all the mstring pairs
    mstringpairs = tku.get_mstring_pairs(args.mfile,args.mstrings)

    ## Adjust seqs if needed (not usually needed)
    seqs,xtra = tku.add_needed_dashes(seqs,mstringpairs)

    if xtra:
        v.print('Sequences expanded to accomondate mstrings')
    if xtra and args.tfile:
        raise RuntimeError(f'Since seqs expanded, cannot trust {args.tfile}')

    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    mut_mgr = mutant.MutationManager(first.seq)

    tweaklist = [tku.tweak_from_mstringpair(mut_mgr,ma,mb)
                 for ma,mb in mstringpairs]

    v.print("\n".join(str(tweak) for tweak in tweaklist))

    for tfile in args.tfile or []:
        tweaklist.extend( IndexTweak.tweaks_from_file(tfile) )

    ## Only need to make a list if there is more than one tweak
    if len(tweaklist)>1:
        v.vvprint("reading sequences...",end="")
        seqs = list(seqs)
        v.vvprint("ok. Read",len(seqs),"sequences")

    output=[]
    for tweak in tweaklist:

        ## oseqs are original seqs, confinded to the relevant ndx's
        oseqs = [s.seq[tweak.ndxlo:tweak.ndxhi] for s in seqs]
        c_oseqs = Counter(oseqs)
        v.vprint('Distinct subseqs:',len(c_oseqs))

        ## count appearances of tweak.sa and tweak.sb
        ca = c_oseqs[tweak.sa]
        cb = c_oseqs[tweak.sb]

        ## Set up the c_aseqs and c_bseqs Counters
        ## to characterize (using entropy, in particular)
        ## how the sequenecs differ under the a-version
        ## vs the b-version of the alignment.
        c_aseqs = tweak_counter(c_oseqs,tweak.sa,tweak.sb)
        c_bseqs = tweak_counter(c_oseqs,tweak.sb,tweak.sa)

        ## only need the indices for which tweak.sa[n] != tweak.sb[n]
        ndiff = [n for n in range(len(tweak.sa))
                 if tweak.sa[n] != tweak.sb[n]]

        Eor,Eos = sitewise_entropy(c_oseqs,ndiff)
        Ear,Eas = sitewise_entropy(c_aseqs,ndiff)
        Ebr,Ebs = sitewise_entropy(c_bseqs,ndiff)

        v.vprint('O:',Eor,Eos)
        v.vprint('A:',Ear,Eas)
        v.vprint('B:',Ebr,Ebs)

        #ma,mb = (tweak.ma,tweak.mb) if hasattr(tweak,'ma') else "",""

        if args.output.endswith('.tsv'):
            output.append("\t".join(map(str,(tweak,ca,cb,Ebs-Eas,
                                             thesign(Ebs-Eas)))))
        else:
            output.append(f'{tweak} {ca:7d} {cb:7d} {Ebs-Eas:+12.9f}')


    with xopen(args.output,"w") as fout:
        if args.output.endswith('.tsv'):
            tsv_header = '\t'.join(['tweak',
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
