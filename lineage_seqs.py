'''
Produce a fasta file (of DNA sequences, ulitmately) which contains
pareent and child sequenecs for every lineage transition in a
newmuts.tsv file.  Want the earliest parent consistent with most
common form of lineage, and is itself of that lineage.  Also, the
earliest child (ie, seq with child lineage) consistent with the given
mutation (needn't be "exact" mutation; ie if it had other mutations,
that's okay)
for example: L5F: BA.2->BA.2.24 is line in newmuts.tsv file
Want earliest instance of most common form of BA.2; and
Earliest instance of BA.2.24 with an L5F mutation.
'''

import re
import argparse
from collections import Counter,defaultdict,namedtuple
import pandas as pd

import verbose as v
import sequtil
import covid
import mutant
import numu
from lineagenotes import LineageNotes
from commonforms import LineagePartition

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    covid.corona_args(argparser)
    paa = argparser.add_argument_group("Lineage Sequence Options").add_argument
    paa("--notesfile",
        help="file with lineage_notes")
    paa("--mutsfile","-M",required=True,
        help="File (typically tsv) with new mutations, by lineage")
    paa("--mincount","-m",type=int,default=0,
        help="only include mutations with minimum count")
    paa("--clade",
        help="Restrict data to this clade")
    paa("--output","-o",
        help="Output file is aa fasta file")
    paa("--dnainput","-D",
        help="reference DNA sequences")
    paa("--dnaoutput",
        help="output DNA sequences")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    if args.notesfile is None:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)

    return args

def earliest_seq(seqlist):
    '''given a list of SequenceSample's (ie: s.name, s.seq),
    return the single SequenceSample that is earliest
    '''
    earlydate=None
    earlyseq = seqlist[0] if len(seqlist) else None
    for s in seqlist:
        date = covid.get_date(s)
        if not date:
            continue
        if not earlydate or date < earlydate:
            earlydate = date
            earlyseq = s
    return earlyseq

class LineageCatalog(LineagePartition):
    '''Catalog of lineages'''
    def __init__(self,fullseqs):
        #super(LineageCatalog,self)
        #LineagePartition.__init__(self,fullseqs,bylineage=True)
        super().__init__(fullseqs,bylineage=True)
        v.vprint('LineagePartition finished...')
        ## for each lineage, list of earliest instance of each sequence pattern
        self.earlyseqs=defaultdict(list)
        self.mcf=dict() ## sequence with most common form for each lineage

        for lin in self.lineages:
            v.vprint_only(5,'lin',f'lin={lin}...')
            ## tmp structure: dict of lists, with
            ## list of sequences for each sequence pattern
            seqlist_bypattern = defaultdict(list)
            for s in self.sequences[lin]:
                seqlist_bypattern[s.seq].append(s)
            counts = Counter({seqpattern: len(seqlist_bypattern[seqpattern])
                              for seqpattern in seqlist_bypattern})
            [(mcf_seqpattern,_)] = counts.most_common(1)
            for seqpattern in counts:
                earlyseq = earliest_seq(seqlist_bypattern[seqpattern])
                if earlyseq is None:
                    v.print('No date for pattern with',
                            counts[seqpattern],'sequences')
                    continue
                self.earlyseqs[lin].append(earlyseq)
                if seqpattern == mcf_seqpattern:
                    self.mcf[lin] = earlyseq

            ## at this point could delete self.sequences dict
            ## since we won't use it

    def report_size(self,lin=None):
        '''return size of lineage catalog'''
        if lin:
            return len(self.earlyseqs[lin])
        return sum(len(self.earlyseqs[lin])
                   for lin in self.earlyseqs)

def get_transitions(mutsfile,cutoff=0):
    '''For all lineage transitions, parse out all the
    the parents and chldren; for children, include the mutation
    associated with it
    '''
    MutLin = namedtuple('MutLin','mut lin')

    df = pd.read_table(mutsfile)
    lintran = numu.match_column_name(df.columns,"lineage_trans")
    linseqs = numu.match_column_name(df.columns,"lineage_seq")
    v.vprint(f'{df.columns}')
    v.vprint(f'lintran={lintran}')
    parents = set()
    children = set()
    if not lintran:
        return parents,children
    for _,row in df.iterrows():
        mut = row["mutation"]
        if int(row[linseqs]) < cutoff:
            v.vprint_only(5,'skiplowcount',mut,row[linseqs])
            continue
        transitions = row[lintran]
        if not transitions or pd.isna(transitions):
            continue
        v.vvvprint_only(5,'trans:',transitions)
        transitions = transitions.split(",")
        transitions = [t.strip() for t in transitions]
        for trans in transitions:
            parent,child = trans.split("->")
            v.vvprint_only(5,'Trans:',f'mut={mut}: parent={parent} child={child}')
            parents.add(parent)
            children.add(MutLin(mut,child))

    v.vprint_only_summary('skiplowcount')
    v.vprint('Parents:',len(parents))
    v.vprint('Children:',len(children))

    return parents,children

def _main(args):
    '''main'''
    v.vprint(args)

    parents,children = get_transitions(args.mutsfile,cutoff=args.mincount)

    for child in children:
        v.vvprint_only(5,'Child:',f'child={child} mut={child.mut}, lin={child.lin}')

    lin_notes = LineageNotes.from_file(args.notesfile)
    lineage_set = lin_notes.get_lineage_set(args.clade)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    mmgr = mutant.MutationManager(first.seq)

    seqs = (s for s in seqs
            if covid.get_lineage(s) in lineage_set)

    lincat = LineageCatalog(seqs)
    v.vprint('lincat:',lincat.report_size(),'sequences')

    newnames = dict()
    outputseqs = set()
    for parent in parents:
        parentseq = lincat.mcf.get(parent,None)
        if not parentseq:
            continue
        newnames[parentseq.name] = "_".join([parent,
                                             lin_notes.parent_of(parent),
                                             covid.get_isl(parentseq.name),
                                             "MCF"])
        outputseqs.add( parentseq )

    v.vprint('Output will have',len(outputseqs),"mcf sequences")

    childseqs = set()
    for ssm,child in children:
        ssm = mutant.SingleSiteMutation.from_string(ssm)
        candidates=[s for s in lincat.earlyseqs[child]
                    if ssm in mmgr.seq_to_mutation(s.seq)]
        earlyseq = earliest_seq(candidates)
        if earlyseq is None:
            v.vprint('No seq for:',ssm,child)
            continue

        if earlyseq.name not in newnames:
            newnames[earlyseq.name] = "_".join([child,
                                                lin_notes.parent_of(child),
                                                covid.get_isl(parentseq.name)])
        newnames[earlyseq.name] += "_" + str(ssm)
        childseqs.add(earlyseq)

    v.vprint('Output will add',len(childseqs),"child sequences")
    v.vprint('with:',len(set(s.name for s in childseqs)),'distinct names')

    for s in childseqs:
        outputseqs.add(s)
    v.vprint('Output will have',len(outputseqs),"total sequences")
    v.vprint('with:',len(set(s.name for s in outputseqs)),'distinct names')
    isl_setofall = set(covid.get_isl(s.name) for s in outputseqs)
    v.vprint('isl_setofall:',len(isl_setofall))

    ## now that newnames is set, write out a table
    outputstem = re.sub(r"\.[^\.]*$","",args.output)
    with open(outputstem+"_table.tsv","w") as fout:
        print("Standard Name","New Name",sep='\t',file=fout)
        for name,newname in newnames.items():
            print(name,newname,sep='\t',file=fout)

    ## now write the aa sequence file
    outputseqlist=[first]
    for s in outputseqs:
        if s.name in newnames:
            s.name = newnames[s.name]
        else:
            v.print(f'Sequence [{s.name}] does not have a new name!')
        outputseqlist.append(s)

    sequtil.write_seqfile(args.output,outputseqlist)
    v.vprint(f'Wrote {len(outputseqlist)} sequences to {args.output}')

    ## Ok, next step is to write out DNA sequences
    if not args.dnainput:
        return


    v.vprint(f"Reading input dna file: {args.dnainput}")
    refseqs = sequtil.read_seqfile(args.dnainput)
    firstref,refseqs = sequtil.get_first_item(refseqs,keepfirst=False)
    outseqs = [firstref]
    dna_isl_set=set()
    for s in refseqs:
        isl_name = covid.get_isl(s.name)
        dna_isl_set.add(isl_name)
        if isl_name in isl_setofall:
            s.name = newnames.get(s.name,s.name)
            outseqs.append(s)
    v.vprint("Read",len(outseqs),"reference sequences")
    for isl in isl_setofall:
        if isl not in dna_isl_set:
            v.vprint(f'Uh-oh! ISL {isl} not in DNA set!')

    if outseqs and args.dnaoutput:
        sequtil.write_seqfile(args.dnaoutput,outseqs)
    v.vprint("Wrote",len(outseqs),
             "reference sequences to file:",args.dnaoutput)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
