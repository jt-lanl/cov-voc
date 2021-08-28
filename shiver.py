'''
SHIVER: SARS CoV-2 Historically Identified Variants in Epitope Regions
'''

import sys
from collections import Counter
import argparse
import numpy as np


import readseq
import sequtil
import intlist
import mutant
from spikevariants import SpikeVariants
import wrapgen
import covid

DESCRIPTION=__doc__ + '''

SHIVER identifies sets of variant forms of the SARS CoV-2 virus with a
focus on just the NTD and RBD neutralizing antibody epitope regions of
the Spike protein, chosen to maximize coverage globally and/or on
separate continents[*], depending on which of several strategies is
employed.

The first variant in the input alignment is taken as the reference
sequence, and should be the ancestral Wuhan variant to ensure epitope
regions are chosen appropriately. The epitope regions in Spike that
are featured as are defined as: The NTD supersite includes Spike
positions 13-20, 140-158, and 242-264 (note, however, that site 18 is
not included in the analysis because it is so variable that both the
ancestral L18 form and the common variant L18F are very often both
found in significant numbers among Variants of Interest).

The NTD supersite sites selected are for inclusion are based on:

Sites 14-20, 140-158, and 245-264:
N-terminal domain antigenic mapping reveals a site of vulnerability for SARS-CoV-2
McCallum, M. et al. bioRxiv
doi: 10.1101/2021.01.14.426475

Site 13:
SARS-CoV-2 immune evasion by variant B.1.427/B.1.429
McCallum, M. et al. bioRxiv, 2021/04/07
doi: 10.1101/2021.03.31.437925 PMC8020983

Sites 242-244:
SARS-CoV-2 501Y.V2 escapes neutralization by South African COVID-19 donor plasma
Wibmer, C. et al. bioRxiv,
doi: 10.1101/2021.01.18.427166

Sites 330-521:
The RBD region includes positions 330-521, based on a synthesis of
the literature from early 2020.

All distinct variants found within these boundaries are identified and
tallied, and the most common variants are selected.  Windows in time
can be selected to reflect more recently emerging patterns in
variation in key epitope regions.

[*] Note that the UK is treated as a separate continent because so much
of the sequencing has been from the UK.
'''

T_STRATEGY='''This run uses the T=taketurns strategy for identifying further
variants.  Each continent, in turn, chooses the next variant, based on
which is the most common variant in that continent that has not
already been chosen.  The order of the continents
is based on number of samples available in those continents.  '''

G_STRATEGY='''This run uses the G=globalonly strategy for identifying further
variants.  Coverage is optimzed globally, without consideration of
continents.'''

M_STRATEGY='''This run uses the M=mostimproved strategy for identifying further
variants.  At each iteration, we determine for each continent how much
gain in fraction coverage would be obtained within that continent if
we chose the variant that maximized that fraction.  We choose the
variant associated with the continent that sees the largest
improvement in fractional coverage.'''

TGM_STRATEGY = {
    'T': T_STRATEGY,
    'G': G_STRATEGY,
    'M': M_STRATEGY,
}

## pattern, motif, design, stencil, prototype? not signature!

TABLE_VARIANTS='''

Table of Variants

In the table below, the first column is the pattern at sites where differences occur,
relative to initial (Wuhan) sequence, with site numbers read down vertically).

LPM = Local Pattern Matches = # of seqs in continent that match over RBD+NTD
GPM = Global Pattern Matches = # of seqs in world that match over RBD+NTD
GSM = Global Sequence Matches = # of seqs that match over whole Spike protein

'''



def _getargs():
    '''get arguments from command line'''
    ap = argparse.ArgumentParser(description=__doc__)
    paa = ap.add_argument
    covid.corona_args(ap)
    paa("--output","-o",
        help="output fasta file with variant sequences")
    paa("-n",type=int,default=29,
        help="number of components in cocktail")
    paa("--strategy","-s",default="taketurns",
        help="how to pick variants: (T)aketurns, (M)ostimproved, (G)lobalonly")
    paa("--region",default="NTD-18+RBD",
        help="region of spike sequence over which patterns are defined")
    paa("--colormut",
        help="name of color mutation file (mutation_string,lineage_name) are 2nd,3rd columns")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = ap.parse_args()
    return args

def print_sequence_counts_by_continent2(Continents,counts):
    '''as a nicely formatted table'''
    #maxlen = max(len(cx) for cx,_,_ in Continents)
    #fmt = "%%-%ds" % maxlen
    print("  #Seqs Continent") #fmt % ("Continent",),"#Seqs")
    for cx,_,_ in [("Global","Global","")] + Continents:
        print("%7d %s" % (sum(counts[cx].values()),cx))

def print_sequence_counts_by_continent(Continents,counts):
    ''' as just a list, that can be part of a paragraph '''
    print("Total: %d" % sum(counts["Global"].values()),end="")
    for cx,_,_ in Continents:
        print(", %s: %d" % (cx,sum(counts[cx].values())),end="")
    print(".")

def main(args):
    '''shiverx main()'''

    strategy_name = { "T": "Taketurns", "G": "Globalonly", "M": "Mostimproved" }
    TGM = args.strategy.upper()[0] ## = T, G, or M
    if TGM not in "TGM":
        raise RuntimeError(f"Strategy [{args.strategy}] should be one of T,G,M")

    allfullseqs = covid.read_seqfile(args)
    allfullseqs = vcount(allfullseqs,"Total sequences:")
    allfullseqs = covid.filter_seqs_by_date(allfullseqs,args)
    allfullseqs = covid.fix_seqs(allfullseqs,args)
    allfullseqs = list(allfullseqs)
    vprint("Read",len(allfullseqs),"sequences")
    fullseqs = covid.filter_seqs_by_pattern(allfullseqs,args)
    fullseqs = list(fullseqs)
    vprint("Read",len(fullseqs),"sequences after filtering")
    firstseq = fullseqs[0].seq

    ## patt seqs is a copy of full seqs,
    ## but with just the sites in restricted region (NTD+RBD)
    sitenums = covid.spike_sites(args.region)

    MM = mutant.MutationManager(firstseq)
    site_indices = []
    for site in sitenums:
        site_indices.extend( MM.indices_from_site(site) )


    pattseqs = sequtil.copy_seqlist(fullseqs)
    for s in pattseqs:
        s.seq = "".join(s.seq[n] for n in site_indices)
    ## Filter out sequences with X's
    pattseqs = [s for s in pattseqs if "X" not in s.seq]
    vprint("Sequences:",len(pattseqs),
           "w/o X in pattern region:",args.region)

    allpattseqs = sequtil.copy_seqlist(allfullseqs)
    for s in allpattseqs:
        s.seq = "".join(s.seq[n] for n in site_indices)
    allpattseqs = [s for s in allpattseqs if "X" not in s.seq]

    firstpatt = pattseqs[0].seq
    pattseqs = pattseqs[1:]
    allpattseqs = allpattseqs[1:]

    global_cnt = Counter(s.seq for s in pattseqs)
    all_global_cnt = Counter(s.seq for s in allpattseqs)

    continent_cnt = dict()
    all_cont_cnt = dict()
    continent_cnt["Global"] = global_cnt
    all_cont_cnt["Global"] = all_global_cnt

    ConExclude = covid.parse_continents()
    for cx,c,x in ConExclude:
        if x:
            def keepseq(s):
                return c in s.name and x not in s.name
        else:
            def keepseq(s):
                return c in s.name

        cseqs = filter(keepseq,pattseqs)
        cnt = Counter(s.seq for s in cseqs)
        continent_cnt[cx] = cnt
        vprint(f"{cx:25s} {sum(cnt.values()):7d}")
        cseqs = filter(keepseq,allpattseqs)
        cnt = Counter(s.seq for s in cseqs)
        all_cont_cnt[cx] = cnt

    ## Sort ConExclude by total counts
    ## so continent w/ most sequences gets to go first
    def conex_sortkey(cxcx):
        cx,_,_ = cxcx
        return sum(continent_cnt[cx].values())
    ConExclude = sorted(ConExclude, key=conex_sortkey, reverse=True)

    vprint("                Continent   #Seqs [Coverage in first five patterns  ]")
    for cx,cnt in continent_cnt.items():
        total = sum(cnt.values()) ## same as len(cseqs)
        vlist = sorted(list(cnt),key=cnt.get,reverse=True)
        vprint("%25s %7d [%s]"%(cx,total,
                               ", ".join(["%.3f"%(cnt[v]/total,)
                                          for v in vlist[:5]])))
    vprint()

    cocktail=[]
    cockname=[]
    cockcont=[]

    includeFirst=True
    if includeFirst:
        cocktail.append(firstpatt)
        cockname.append("1-Initial")
        cockcont.append("Global")


    def check_next_variant(cnt):
        vlist = sorted(list(cnt),key=cnt.get,reverse=True)
        for v in vlist:
            if v in cocktail:
                continue
            return v,cnt[v]/sum(cnt.values())
        return None,None


    def make_name(c):
        n = len(cocktail)
        nc = sum(1 for cname in cockname if c in cname)
        cc = "Local" if c=="Global" and args.filterbyname else c
        return "%d-%s-%d" % (n,cc,nc+1)

    def alt_get_next_variant():
        candidates = []
        for cx,c,x in ConExclude:
            v,val = check_next_variant(continent_cnt[cx])
            if v:
                candidates.append( (val,v,cx,c,x) )

        if candidates:
            val,v,cx,c,x = max(candidates)
            vprint("Appending v from",cx,"with value=",val)
            cocktail.append(v)
            cockcont.append(cx)
            cockname.append(make_name(c))
            return v
        return None


    def get_next_variant(cnt,cx,c):
        vlist = sorted(list(cnt),key=cnt.get,reverse=True)
        for v in vlist:
            if v in cocktail:
                ndx = cocktail.index(v)
                vvprint("Skipping v with cnt=",cnt[v],"#",ndx,cockname[ndx])
                continue
            vvprint("Appending v with cnt=",cnt[v],global_cnt[v])
            cocktail.append(v)
            cockcont.append(cx)
            cockname.append(make_name(c))
            return v
        return None

    while len(cocktail) < args.n:
        if TGM == "G":
            v = get_next_variant(global_cnt,"Global","Global")
        elif TGM == "M":
            v = alt_get_next_variant()
        elif TGM == "T":
            v = None  # only if all get_next's are None will v stay equal to None
            for cx,c,x in ConExclude:
                vv = get_next_variant(continent_cnt[cx],cx,c)
                v = v or vv
                if len(cocktail) >= args.n:
                    break
        else:
            raise RuntimeError(f"Invalid strategy: {args.strategy}")
        if v is None:
            break

    vvprint(sitenums)

    for line in intlist.write_numbers_vertically(sitenums,plusone=0):
        vvprint("       %s" % (line,))

    vo = cocktail[0]
    for v in cocktail:
        rseq = sequtil.relativename(vo,v) if v!=vo else vo
        vvprint("%6d %s" % (global_cnt[v],rseq))

    altered_sites= np.zeros(shape=(len(vo),),dtype=bool)
    for v in cocktail[1:]:
        altered = np.array([v[n] != vo[n] for n in range(len(v))])
        altered_sites |= altered

    vprint("altered sites:",altered_sites.shape)
    vprint("sitenums:",len(sitenums))
    vprint("vo:",len(vo))

    altered_site_indices = [site_indices[n] for n in range(len(vo)) if altered_sites[n]]
    altered_sitenums = [MM.site_from_index(n) for n in altered_site_indices]
    vprint("Altered sites:",altered_sitenums)

    if not altered_sitenums and args.n > 1:
        ## Should only happen if all sequences are identical in RBD/NTD regions
        raise RuntimeError("No altered sites: must be an error!")


    rows = []
    srseq = dict()
    if altered_sitenums:
        vertlines = intlist.write_numbers_vertically(altered_sitenums,plusone=0)
        for line in vertlines:
            title = "Global  Local" if line == vertlines[-1] else ""
            vprint("%13s %s" % (title,line))
        for v,c,name in zip(cocktail,cockcont,cockname):
            rseq = sequtil.relativename(vo,v) if v!=vo else vo
            srseq[v] = "".join(rseq[n] for n in range(len(rseq)) if altered_sites[n] )
            rows.append((global_cnt[v],continent_cnt[c][v],srseq[v],name))
        for r in rows:
            vprint("%6d %6d %s %s" % r)
        for line in vertlines:
            title = "Global  Local" if line == vertlines[-1] else ""
            vprint("%13s %s" % (title,line))
        for r in sorted(rows,reverse=True):
            vprint("%6d %6d %s %s" % r)
    else:
        for v,c,name in zip(cocktail,cockcont,cockname):
            srseq[v] = ""
        vertlines=[""]

    ## OK, now that we have our cocktail, can we expand to what the full
    ## spike sequences would be

    ## But first lets build a way to map mutation patterns to lineage names
    lineages = []

    svar = SpikeVariants.from_colormut(args.colormut) if args.colormut \
        else SpikeVariants.default(refseq=firstseq)
    lineages = list(svar.vocs)
    vprint("lineages:",lineages)

    variant_table = dict()
    cocktail_fasta = []
    vo = cocktail[0]
    for v,c,name in zip(cocktail,cockcont,cockname):
        ## list of all seqeunces that match the viral pattern
        vnames = set(s.name for s in pattseqs if v in s.seq)
        vseqs = (s for s in fullseqs if s.name in vnames)
        vcnt = Counter(s.seq for s in vseqs)
        vcntlist = sorted(list(vcnt),key=vcnt.get,reverse=True)
        vvprint(name,":",[vcnt[vc] for vc in vcntlist[:5]])
        if v==vo:
            v_fullseq = firstseq
        elif len(vcntlist)>0:
            ## Take the most common global form
            v_fullseq = vcntlist[0]
        else:
            ## There is no global form: just use X's then
            v_fullseq = "".join("X" for _ in firstseq)

        ## make full list of mutations
        mutliststr = str( MM.get_mutation(v_fullseq,exact=False) )
        vvprint(v,mutliststr)

        ## find best lineage name:
        matched_voc_names = [voc.name for voc in lineages
                             if MM.seq_fits_pattern(voc,v_fullseq)]
        mut_lineage = "".join(f"({name})" for name in matched_voc_names)

        v_fullseq_name = name

        variant_table[v] = "%s %-20s %6d %6d %6d   %5.1f%% %s %s" % (
            srseq[v],
            name,
            continent_cnt[c][v],
            global_cnt[v],
            vcnt[v_fullseq],
            100*vcnt[v_fullseq]/global_cnt[v] if vcnt[v_fullseq]>0 else 0,
            mutliststr,
            mut_lineage,
        )
        cocktail_fasta.append(
            readseq.SequenceSample(v_fullseq_name,
                                   v_fullseq))

    print(DESCRIPTION)
    print()
    print(TGM_STRATEGY[TGM])
    print()
    print("This run uses sequences sampled from %s to %s."
          % sequtil.range_of_dates(pattseqs))
    if args.filterbyname:
        print("Filtered by geographic region(s):","+".join(args.filterbyname))
    if args.xfilterbyname:
        print("Excluding geographic region(s):","+".join(args.xfilterbyname))
    print("The number of sequences, broken out by continent is:")
    print_sequence_counts_by_continent(ConExclude,continent_cnt)

    print("Note: the focus here is specifically on the epitope region:",args.region)
    print("Sites:",intlist.intlist_to_string(sitenums,sort=True))

    print(TABLE_VARIANTS)
    TabVar_Heading="Name                    LPM    GPM    GSM  GSM/GPM [Mutations] %s" \
        % ("(Lineage)" if args.colormut else "")
    for line in vertlines[:-1]:
        print(line)
    print(vertlines[-1],TabVar_Heading)
    for v in variant_table.values():
        print(v)
    #if not args.colormut:
    #    print("\n* Note: lineages not available with this run")

    if args.output:
        fname = args.output
        fname = fname if fname.endswith(".fasta") else fname + ".fasta"
        readseq.write_fasta(fname,cocktail_fasta)

    dropFirst=False
    if dropFirst:
        cocktail = cocktail[1:]

    def coverage_table(cocktail_array):
        ### Cocktail coverage per continent
        ### (fraction of sequences witht exact match in RBD/NTD regions)
        CocktailName = f"{TGM}-{len(cocktail_array)}"
        for cx,_,_ in [("Global","Global","")] + ConExclude:
            cnt = all_cont_cnt[cx]  ## not continent_cnt[cx]
            tot = sum(cnt.values())
            cov = sum(cnt[v] for v in cocktail_array)
            if tot>0:
                f = cov/tot
                yield "%25s %-4s   %.4f" % (cx,CocktailName,f)

    if len(allpattseqs) > len(pattseqs):
        print("\nNote: coverage plot below is based on possibly larger sequence set:")
        print_sequence_counts_by_continent(ConExclude,all_cont_cnt)

    COVERAGETABLE = f'''
Table of Coverages

In table below, {TGM}-n refers to a batch of the first n variants.
Coverage is defined as fraction of sequences in the continent with an
exact match (over the region {args.region}) to one of the first n variants.
(Here, '{TGM}' corresponds to the {strategy_name[TGM]} strategy.)
The coverage table is based on {len(allpattseqs)} sequences.
'''

    print(COVERAGETABLE)

    print("%25s %-4s %s" % ("Continent","Name","Coverage"))
    for n in range(1,len(cocktail)+1,7):
        print()
        for line in coverage_table(cocktail[:n]):
            print(line)


if __name__ == "__main__":

    _args = _getargs()
    def vprint(*p,**kw):
        '''verbose print'''
        if _args.verbose:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vvprint(*p,**kw):
        '''very verbose print'''
        if _args.verbose>1:
            print(*p,file=sys.stderr,flush=True,**kw)
    def vcount(seqs,*p,**kw):
        '''if verbose, count seqs as they go by'''
        return wrapgen.keepcount(seqs,*p,**kw) if _args.verbose else seqs

    main(_args)
