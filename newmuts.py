'''
Find parent/child lineage pairs for which new mutations appear

1. From sequence file, obtain all the pango lineages
2. For each lineage, get most common sequence string
3.     Get mutation list (vis-a-vis Wuhan) for that string
4. For each lineage, get parent, and parent's mutation list
5.     Identify new mutations vis-a-vis the parent
6. While doing all this, be collecting a list of mutations
7. For each mutation, identify child/parent pairs, including counts
'''
## note, consensus is the most expensive part of the computation
## use --consensusnever to avoid that computation

import re
from collections import Counter,defaultdict
from functools import lru_cache
import argparse

import verbose as v

import covid
import mutant
import commonforms as cf
from xopen import xopen
from lineagenotes import LineageNotes

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    covid.corona_args(ap)
    paa = ap.add_argument_group("New Mutations Options").add_argument
    paa("--notesfile",
        help="lineage_notes.txt")
    paa("--cutoff","-c",type=int,default=0,
        help="Do not count any lineage mutations that occur fewer then c times")
    paa("--clade",
        help="Restrict attention to sequences in clade's lineage")
    paa("--mincount","-m",type=int,default=0,
        help="Show only lineage mutations with at least this many appearances")
    paa("--skipwuhan",action="store_true",
        help="skip mutations where Wuhan is the parent")
    paa("--mutationsfile","-M",#default="-",
        help="Write mutations to a tsv file")
    paa("--reversionsfile","-R",
        help="Write reversions to a tsv file")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    if args.notesfile is None:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)
    return args

def most_common_forms(seqlist,mut_manager):
    '''return a dict with most common form for each lineage'''
    ## mcf expressed as a set of mutations

    mcf_dict = dict() ## most common form for each lineage

    ## Partition seqlist by lineages, separate list for each lineage
    lin_partition = cf.LineagePartition(seqlist)

    for lin,seqlin in lin_partition.sequences.items():
        ## Get most common form
        cntr = Counter(s.seq for s in seqlin)
        mcf,_ = cntr.most_common(1)[0]
        ## Express as mutations
        mcf_muts = mut_manager.get_mutation(mcf)
        mcf_dict[lin] = set(mcf_muts)

    return mcf_dict

@lru_cache(maxsize=None)
def mutations_fromseq(mut_manager,seq):
    '''convert sequence into a set of mutations'''
    return set( mut_manager.get_mutation(seq) )

def count_appearances(lin_notes,mutsdict):
    '''parse mutsdict dict to count apppearances'''
    count_seq_child = Counter() ## how many sequences (of child lineage) have mutation
    count_seq_total = Counter() ## how many of ALL sequences have mutation
    count_lin_child = Counter() ## how many disctinct lineages have mutation
    count_lin_parent = Counter() ## how many distinct parents
    for mut in mutsdict:
        pset = set()
        for lin,cnt,tcnt in mutsdict[mut]:
            pset.add(lin_notes.parent_of(lin))
            count_seq_child[mut] += cnt
            count_lin_child[mut] += 1
            count_seq_total[mut] = tcnt
        count_lin_parent[mut] = len(pset)

    return count_seq_child, count_seq_total, count_lin_child, count_lin_parent

def write_mutations_summary(filename,lin_notes,mutsdict,denominator,mincount=0):
    '''write tsv file summarizing the mutations'''
    if not filename:
        return
    (count_seq_child,
     count_seq_total,
     count_lin_child,
     count_lin_parent) = count_appearances(lin_notes,mutsdict)

    header = "\t".join(["site",
                        "mutation",
                        "parent_lineages",
                        "child_lineages",
                        "lineage_sequences",
                        "total_sequences",
                        "denominator",
                        "lineage_transitions"
                        ])
    fmt = "\t".join("%d %s %d %d %d %d %d %s".split())
    with xopen(filename,"w") as fout:
        print(header,file=fout)
        v.vprint(f'Writing {len(mutsdict)} mutations to file={filename}')
        for mut in sorted(mutsdict.keys()):
            if count_seq_child[mut] < mincount:
                continue
            try:
                site = int(re.sub(r'\D*(\d+)\D*',r'\1',str(mut)))
            except TypeError as err:
                v.print(f'mut={mut}')
                raise TypeError from err
            lintrans = ", ".join(f'{lin_notes.parent_of(lin)}->{lin}'
                                for lin,_,_ in mutsdict[mut])
            print(fmt % (site,mut,count_lin_parent[mut],
                         count_lin_child[mut],count_seq_child[mut],
                         count_seq_total[mut],denominator,lintrans),
                  file=fout)


def main(args):
    '''newmuts main'''
    v.vprint(args)

    lin_notes = LineageNotes.from_file(args.notesfile)
    lineage_set = lin_notes.get_lineage_set(args.clade)

    mut_appearances = defaultdict(list)
    mut_reversions = defaultdict(list)

    firstseq,seqlist = cf.get_input_sequences(args)
    mut_manager = mutant.MutationManager(firstseq)
    mcf_dict = most_common_forms(seqlist,mut_manager)
    lin_partition = cf.LineagePartition(seqlist)

    denominator = sum(1 for s in seqlist
                      if covid.get_lineage(s) in lineage_set)
    v.vprint(f'denominator={denominator}/{len(seqlist)}')

    mut_total  = Counter() ## seqs with mut, regardless of lin
    for lin in lineage_set: ## only consider lineages in clade
        parent = lin_notes.parent_of(lin)
        if args.skipwuhan and parent=="Wuhan":
            continue
        parmutset = mcf_dict.get(parent,set([]))
        seqlin = lin_partition.sequences[lin]
        cntr = Counter(s.seq for s in seqlin)
        mut_appear = Counter() ## seqs where child has mat, parent does not
        mut_revert = Counter() ## seqs where parent has mut, child does not
        for seq,cnt in cntr.items():
            mutset = set( mut_manager.get_mutation(seq) )
            ## cull out the X's
            mutset = set(mut for mut in mutset if "X" not in str(mut))
            for mut in mutset:
                mut_total[mut] += cnt
            if cnt < args.cutoff:
                continue
            new_mutations = mutset - parmutset
            for mut in new_mutations:
                mut_appear[mut] += cnt
            rev_mutations = parmutset - mutset
            for mut in rev_mutations:
                mut_revert[mut] += cnt

        for mut,mcnt in mut_appear.items():
            mut_appearances[mut].append((lin,mcnt,mut_total[mut]))
        for mut,mcnt in mut_revert.items():
            mut_reversions[mut].append((lin,mcnt,mut_total[mut]))

    write_mutations_summary(args.mutationsfile,
                            lin_notes,mut_appearances,
                            denominator,
                            mincount=args.mincount)
    write_mutations_summary(args.reversionsfile,
                            lin_notes,mut_reversions,
                            denominator,
                            mincount=args.mincount)



if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
