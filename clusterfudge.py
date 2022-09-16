'''Ad-hoc approach to clustering sequences, 
   Groups are split recursively,
   based on the site with the most mutations.
'''
import re
import argparse
from collections import Counter
from verbose import verbose as v
import sequtil

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="Input fasta file, with first seq the reference")
    paa("--output","-o",
        help="Output fasta file")
    paa("--minclustersize",type=int,default=2,
        help="Stop clustering if clusters are smaller than this")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

class Cluster(list):
    '''a cluster of sequences is just a list with a few flags'''
    def __init__(self,seqs=None):
        if seqs:
            self.extend(seqs)
        self.splittable=True
        self.definedby=[]
        self.name=""
        
    def __str__(self):
        return (f'Cluster {self.name}: split={self.splittable}, '
                f'len={len(self)}, definedby: {self.definedby}')
    
def split_seqs(first,seqs,sitelist=None):
    sitelist = sitelist or range(len(first.seq))
    muts_by_site = Counter()
    sitelist = [n for n in sitelist
                if first.seq[n] in 'ACGT']
    for s in seqs:
        muts_by_site.update((n,s.seq[n]) for n in sitelist
                            if first.seq[n] != s.seq[n]
                            and s.seq[n] in 'ACGT')
    try:
        (sitesplit,mut),count = muts_by_site.most_common(1)[0]
        v.vprint(f'Site {sitesplit} has {count}/{len(seqs)} mutations '
                 f'of {first.seq[sitesplit]}->{mut}')
        sitesplit_ties = [site for (site,base),cnt in muts_by_site.items()
                          if cnt==count]
        v.vprint('Ties:',sitesplit,'=',sitesplit_ties)
    except IndexError:
        v.print('muts:',muts_by_site)
        v.vprint(f'No new mutations in cluster of {len(seqs)} sequences')
        sitesplit=None
    
    aseqs = Cluster()
    bseqs = Cluster()
    for s in seqs:
        if ((sitesplit is not None) and
            (first.seq[sitesplit] != s.seq[sitesplit]) and
            (s.seq[sitesplit] in 'ACGT')):
            aseqs.append(s)
        else:
            bseqs.append(s)

    return aseqs,bseqs,sitesplit

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = sequtil.read_seqfile(args.input)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    seqs = list(seqs)

    sitelist = range(len(first.seq))
    clusters = [Cluster(seqs)]  ## initialize list with one cluster
    while True:
        newclusters = []
        v.print([len(cluster) for cluster in clusters])
        for cluster in clusters:
            v.vprint(cluster)
            if not cluster.splittable:
                newclusters.append(cluster)
                continue
            aseqs,bseqs,sitesplit = split_seqs(first,cluster,sitelist)
            v.vprint(f'{sitesplit=}')
            if len(aseqs) < args.minclustersize:
                ## then don't split
                cluster.splittable=False
                newclusters.append(cluster)
            elif len(bseqs) == 0:
                ## don't split this time, but try again later
                sitelist = sorted(set(sitelist)-set([sitesplit]))
                newclusters.append(cluster)                
            else:
                ## split
                sitelist = sorted(set(sitelist)-set([sitesplit]))
                aseqs.name = cluster.name + "A"
                if cluster.definedby:
                    aseqs.definedby.extend(cluster.definedby)
                    bseqs.definedby.extend(cluster.definedby)
                    
                aseqs.definedby.append("%s%d%s" % (first.seq[sitesplit],
                                                sitesplit,
                                                aseqs[0].seq[sitesplit]))
                v.vprint(aseqs.name,sitesplit,aseqs.definedby)
                bseqs.name = cluster.name + "B"
                newclusters.extend([aseqs,bseqs])
        if len(newclusters) == len(clusters):
            break
        clusters = newclusters

    maxlenclustername = max(len(cluster.name) for cluster in clusters)
    fmtclustername = "%%%ds" % maxlenclustername


    ## get class counts
    if len(clusters) == 1:
        raise RuntimeError(f'No clustering: no clusters > {args.minclustersize}')

    
    classcounts = Counter()
    class_number = 1
    for cluster in clusters:
        v.vprint('get class counts:',cluster)
        if cluster.name[class_number-1]=='B':
            class_number += 1
        classcounts[class_number] += len(cluster)
    
    seqs = [first]
    class_number = 1
    print("Class %2d: (%d sequences)" % (class_number,classcounts[class_number]))
    print("%s %6s %s" % ((fmtclustername % 'name'),'count','defined_by'))
    for n,cluster in enumerate(clusters):
        if cluster.name[class_number-1]=='B':
            class_number += 1
            print("Class %2d: (%d sequences)" % (class_number,classcounts[class_number]))
            print("%s %6s %s" % ((fmtclustername % 'name'),'count','defined_by'))
        print("%s %6d" % ((fmtclustername % cluster.name),len(cluster)),end=" ")
        print("[%s]" % ",".join(cluster.definedby))
        for s in cluster:
            print(fmtclustername % "",s.name)
        seqs.extend(cluster)
    
    if args.output:
        sequtil.write_seqfile(args.output,seqs)    

    

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
