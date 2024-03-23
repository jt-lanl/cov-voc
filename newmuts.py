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

from collections import Counter,defaultdict
import argparse

import verbose as v
import breakpipe

import covid
import mutant
from xopen import xopen
from lineagenotes import LineageNotes
from commonforms import LineagePartition

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
        help="Restrict attention to sequences in clade's lineages")
    paa("--xclade",
        help="Avoid sequences in xclade's lineages")
    paa("--mincount","-m",type=int,default=0,
        help="Show only lineage mutations with at least this many appearances")
    paa("--skipwuhan",action="store_true",
        help="skip mutations where Wuhan is the parent")
    paa("--mutationsfile","-M",#default="-",
        help="Write mutations to a tsv file")
    paa("--reversionsfile","-R",
        help="Write reversions to a tsv file")
    paa("--output","-o",
        help="Write two tsv files: equivealent to -M new-OUTPUT -R rev-OUTPUT ")
    args = ap.parse_args()
    args = covid.corona_fixargs(args)
    if args.notesfile is None:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)
    if args.output:
        if not args.mutationsfile:
            args.mutationsfile = covid.filename_prepend("new-",args.output)
        if not args.reversionsfile:
            args.reversionsfile = covid.filename_prepend("rev-",args.output)
    if args.clade and "-" in args.clade:
        if args.xclade:
            raise RuntimeError('Cannot use minus (-) in clade and also have an xclade')
        args.clade,args.xclade = args.clade.split('-')

    return args

def most_common_forms(mut_manager,lin_partition):
    '''return a dict with most common form for each lineage'''
    ## mcf expressed as a set of mutations

    mcf_dict = dict() ## most common form for each lineage

    for lin,seqlin in lin_partition.sequences.items():
        ## Get most common form
        cntr = Counter(s.seq for s in seqlin)
        mcf,_ = cntr.most_common(1)[0]
        ## Express as mutations
        mcf_muts = mut_manager.seq_to_mutation(mcf)
        mcf_dict[lin] = set(mcf_muts)

    return mcf_dict

def count_appearances(lin_notes,mutsdict):
    '''parse mutsdict dict to count apppearances'''
    count_seq_child = Counter() ## how many sequences (of child lineage) have mutation
    count_lin_child = Counter() ## how many disctinct lineages have mutation
    count_lin_parent = Counter() ## how many distinct parents
    for ssm in mutsdict:
        pset = set()
        for lin,cnt in mutsdict[ssm]:
            pset.add(lin_notes.parent_of(lin))
            count_seq_child[ssm] += cnt
            count_lin_child[ssm] += 1
        count_lin_parent[ssm] = len(pset)

    return count_seq_child, count_lin_child, count_lin_parent

def write_mutations_summary(filename,lin_notes,mutsdict,count_seq_total,denominator,mincount=0):
    '''write tsv file summarizing the mutations'''
    if not filename:
        return
    (count_seq_child,
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
        for ssm in sorted(mutsdict):
            if count_seq_child[ssm] < mincount:
                continue
            lineage_transitions = ", ".join(f'{lin_notes.parent_of(lin)}->{lin}'
                                            for lin,_ in mutsdict[ssm])
            print(fmt % (ssm.site,ssm,count_lin_parent[ssm],
                         count_lin_child[ssm],count_seq_child[ssm],
                         count_seq_total[ssm],denominator,lineage_transitions),
                  file=fout)

@breakpipe.no_broken_pipe
def main(args):
    '''newmuts main'''
    v.vprint(args)

    lin_notes = LineageNotes.from_file(args.notesfile)
    lineage_set = lin_notes.get_lineage_set(args.clade)
    if args.xclade:
        xlineage_set = lin_notes.get_lineage_set(args.xclade)
        lineage_set = lineage_set - xlineage_set
    v.vprint('Total lineages:',len(lineage_set))

    if args.clade and args.clade in lineage_set:
        # Except for no clade (or clade=ALL),
        # We don't want to include parent of clade to clade transitions
        args.skipwuhan=True

    ssm_appearances = defaultdict(list)
    ssm_reversions = defaultdict(list)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs)
    mut_manager = mutant.MutationManager(first.seq)

    lin_partition = LineagePartition(seqs,restrict_to=lineage_set)
    denominator = sum(lin_partition.counts.values())
    v.vprint(f'denominator={denominator}')
    v.vprint(f'lineages in partition: {len(lin_partition.lineages)}')

    mcf_dict = most_common_forms(mut_manager,lin_partition)


    ssm_total  = Counter() ## seqs with ssm, regardless of lin
    for lin in lineage_set: ## only consider lineages in clade
        parent = lin_notes.parent_of(lin)
        if parent not in lineage_set:
            ## this happens when lin is the clade lineage
            parent="Wuhan"
        if args.skipwuhan and parent=="Wuhan":
            continue
        parmutset = mcf_dict.get(parent,set())
        seqlin = lin_partition.sequences[lin]
        cntr = Counter(s.seq for s in seqlin)
        ssm_appear = Counter() ## seqs where child has mat, parent does not
        ssm_revert = Counter() ## seqs where parent has ssm, child does not
        for seq,cnt in cntr.items():
            mutset = set( mut_manager.seq_to_mutation(seq) )
            ## cull out the X's
            mutset = set(ssm for ssm in mutset if "X" not in ssm.mut)
            for ssm in mutset:
                ssm_total[ssm] += cnt
            if cnt < args.cutoff:
                continue
            new_mutations = mutset - parmutset
            for ssm in new_mutations:
                ssm_appear[ssm] += cnt
            rev_mutations = parmutset - mutset
            for ssm in rev_mutations:
                ssm_revert[ssm] += cnt

        for ssm,mcnt in ssm_appear.items():
            ssm_appearances[ssm].append((lin,mcnt))
        for ssm,mcnt in ssm_revert.items():
            ssm_reversions[ssm].append((lin,mcnt))

    write_mutations_summary(args.mutationsfile,
                            lin_notes,ssm_appearances,
                            ssm_total,denominator,
                            mincount=args.mincount)
    write_mutations_summary(args.reversionsfile,
                            lin_notes,ssm_reversions,
                            ssm_total,denominator,
                            mincount=args.mincount)



if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
