'''
takes one or more mutant strings (either from a file or from the
command line), each of which looks something like
"[A222V,A262S,S494P,D614G]", and produces a fasta file, each sequence
of which corresponds to a mutant specified by the string, relative to
the reference sequence, which is the first sequence in the specified
reference fasta file.
'''

import re
import argparse

import verbose as v
import xopen
import sequtil
import wildtypes
import mutant
import tweak as tku

def _getargs():
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    paa = ap.add_argument_group("Mutation to Sequence Options").add_argument
    paa("--mstrings","-m",nargs='+',
        help="one or more strings of the form '[L5F,...,Q957R]'")
    paa("--mfile","-M",
        help="input file with list of mutant m-strings")
    paa("--nodash",action="store_true",
        help="remove dashes from the output sequences")
    paa("--skipref",action="store_true",
        help="Do not include reference sequence in output")
    paa("--output","-o",
        help="ouptut fasta file of mutated seqeunces")
    args = ap.parse_args()
    return args

def rd_mfile(filename):
    '''read file with mstrings'''
    with xopen.xopen(filename) as fin:
        yield from xopen.nonempty_lines(fin)

def _main(args):

    def de_dash(s):
        '''remove dashes from string'''
        if args.nodash:
            s.seq = re.sub("-","",s.seq)
        return s

    mstringlist = []
    if args.mfile:
        mstringlist.extend( rd_mfile(args.mfile) )
    if args.mstrings:
        mstringlist.extend( args.mstrings )

    v.vprint(f'Read {len(mstringlist)} m-strings')
    for mstring in mstringlist:
        v.vvprint("  ",mstring)

    if not mstringlist:
        v.print("Must specify mstrings on command line "
                "with --mfile/-M and/or --mstrings/-m")
        return

    ## Grab the Wuhan base sequence
    wuhan = sequtil.SequenceSample( wildtypes.SEQUENCE_NAME,
                                    wildtypes.SEQUENCE )

    xpander = tku.ExpandSeq(mstringlist,wuhan.seq)
    wuhan.seq = xpander.expand_seq(wuhan.seq)
    mut_mgr = mutant.MutationManager(wuhan.seq)

    mutseqs = []
    if not args.skipref:
        mutseqs.append(wuhan)
    for mstring in mstringlist:
        mutseq = mut_mgr.regex_from_mstring(mstring,exact=True)
        mutseqs.append( sequtil.SequenceSample( mstring, mutseq ) )

    mutseqs = [de_dash(s) for s in mutseqs]

    if args.output:
        sequtil.write_seqfile(args.output,mutseqs)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
