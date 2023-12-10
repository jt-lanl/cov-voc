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
from functools import cache
import argparse

import verbose as v

import covid
import mutant
import commonforms as cf
from readseq import xopen
from lineagenotes import LineageNotes

DEFAULT_LINEAGE_NOTES_FILE="/home/jt/src/corona/data/lineage_notes.txt"

def getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--notesfile",default=DEFAULT_LINEAGE_NOTES_FILE,
        help="lineage_notes.txt")
    paa("--cutoff","-c",type=int,default=0,
        help="Do not count any lineage mutations that occur fewer then c times")
    paa("--clade",
        help="Restrict attention to sequences in clade's lineage")
    paa("--mincount","-m",type=int,default=0,
        help="Show only lineage mutations with at least this many appearances")
    paa("--mutationsfile","-M",default="-",
        help="Write mutations to file")
    paa("--reversionsfile","-R",
        help="Write reversions to file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
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
    
@cache
def mutations_fromseq(mut_manager,seq):
    return set( mut_manager.get_mutation(seq) )

def count_appearances(lin_notes,muts):
    '''parse muts dict to count apppearances'''
    mcount = Counter() ## how many sequences (of daughter lineage) have mutation
    mtotal = Counter() ## how many of ALL sequences have mutation
    mtimes = Counter() ## how many disctinct lineages have mutation
    ptimes = Counter() ## how many distinct parents
    for m in sorted(muts.keys()):
        pset = set()
        for lin,cnt,tcnt in muts[m]:
            pset.add(lin_notes.parent_of(lin))
            mcount[m] += cnt
            mtimes[m] += 1
            mtotal[m] = tcnt
        ptimes[m] = len(pset)

    return mcount, mtotal, mtimes, ptimes
    
def write_mutations(filename,lin_notes,muts,mincount=0):
    mcount,mtotal,mtimes,ptimes = count_appearances(lin_notes,muts)
    with xopen(filename,"w") as fout:
        for m in sorted(muts.keys()):
            for lin,mcnt,tcnt in muts[m]:
                if mcnt < mincount:
                    continue
                print("%10s %3d %3d %7d %7d %6d %12s->%s" %
                      (m,ptimes[m],mtimes[m],mtotal[m],mcount[m],
                       mcnt,
                       lin_notes.parent_of(lin),lin),
                      file=fout)

    

def main(args):
    '''newmuts main'''

    lin_notes = LineageNotes(args.notesfile)
    for bad in lin_notes.inconsistencies(remove=True):
        v.print(bad)
    v.vprint(lin_notes.report_size())

    lineage_set = set(lin_notes.lineages)
    if args.clade and args.clade in lineage_set:
        lin_notes_x = lin_notes.restrict_to_clade(args.clade,
                                                  excludeparent=True)
        lineage_set = set(lin_notes_x.lineages)
    
    mut_appearances = defaultdict(list)
    mut_reversions = defaultdict(list)
    
    firstseq,seqlist = cf.get_input_sequences(args)
    mut_manager = mutant.MutationManager(firstseq)
    mcf_dict = most_common_forms(seqlist,mut_manager)
    #lin_info = get_lineage_info(firstseq,seqlist)
    lin_partition = cf.LineagePartition(seqlist)
            
    mut_total  = Counter() ## seqs with mut, regardless of lin
    for lin in lin_partition.lineages:
        if lin not in lineage_set:
            continue
        parent = lin_notes.parent_of(lin)
        parmut = mcf_dict.get(parent,set([]))
        seqlin = lin_partition.sequences[lin]
        cntr = Counter(s.seq for s in seqlin)
        mut_appear = Counter() ## seqs where child has mat, parent does not
        mut_revert = Counter() ## seqs where parent has mut, child does not
        for seq,cnt in cntr.items():
            mut = set( mut_manager.get_mutation(seq) )
            for m in mut:
                mut_total[m] += cnt
            if cnt < args.cutoff:
                continue
            new_mutations = mut - parmut
            for m in new_mutations:
                mut_appear[m] += cnt
            rev_mutations = parmut - mut
            for m in rev_mutations:
                mut_revert[m] += cnt

        for m,mcnt in mut_appear.items():
            mut_appearances[m].append((lin,mcnt,mut_total[m]))
        for m,mcnt in mut_revert.items():
            mut_reversions[m].append((lin,mcnt,mut_total[m]))

    write_mutations(args.mutationsfile,
                    lin_notes,mut_appearances,
                    mincount=args.mincount)
    write_mutations(args.reversionsfile,
                    lin_notes,mut_reversions,
                    mincount=args.mincount)
                  

if __name__ == "__main__":

    _args = getargs()
    v.verbosity(_args.verbose)

    main(_args)
