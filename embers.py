'''
Stacked barplots (also, optionally, line-plots) of variant counts over time
'''
from collections import Counter
import argparse
import warnings

import sequtil
import intlist
from spikevariants import SpikeVariants
import covid

from verbose import verbose as v
import embersutil as emu

def _getargs():
    ''' read command line arguments using argparse package'''
    ap = argparse.ArgumentParser(description=__doc__,
                                 conflict_handler='resolve')
    paa = ap.add_argument
    covid.corona_args(ap)
    emu.embers_args(ap)
    paa("--colormut","-c",
        help="read SpikeVariants structure from color-mut file")
    paa("--ctable","-t",
        help="write a count table to this file")
    paa("--verbose","-v",action="count",default=0,
        help="verbosity")
    ap.set_defaults(legend=0)
    args = ap.parse_args()
    return args

def relativepattern(master,mutant,dittochar='_',noditto='-'):
    '''AB-C,AB-D -> __-D'''
    return "".join((dittochar if (a==b and a not in noditto) else b)
                   for a,b in zip(master,mutant))

def reltoabspattern(master,mutant,dittochar='_'):
    '''ABC,__D -> ABD'''
    return "".join((a if b==dittochar else b)
                   for a,b in zip(master,mutant))


def check_dups(xlist):
    '''see if there are any duplicates in the list xlist'''
    duplist=[]
    xset=set()
    for x in xlist:
        if x in xset:
            duplist.append(x)
        xset.add(x)
    return duplist

def lineage_counts(sitelist,master,voclist,cpatt,n_sequences):
    ''' yield print strings for lineage counts '''
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent Lineage"
    for voc in voclist[::-1]:
        patt = voc.flat_pattern
        fpatt = cpatt[patt]/n_sequences if n_sequences else 0
        rpatt = relativepattern(master,patt)
        yield "%s %6d %6.2f%% %s" % (rpatt,cpatt[patt],100*fpatt,voc.name)

def missing_patterns_with_nearby(sitelist,master,voclist,xpatt,n_sequences):
    ''' for each of the seq patterns in xpatt, find nearby patterns in voclist '''
    ## routine yields lines that are meant to be printed
    yield f"Missing: {len(xpatt)} patterns, for a total of {sum(xpatt.values())} sequences"
    yield from intlist.write_numbers_vertically(sitelist)
    yield f"{master} Counts Percent NearbyLineage(s)"
    for seq in sorted(xpatt,key=xpatt.get,reverse=True)[:50]: #hardcoded 15
        rseq = relativepattern(master,seq)
        dv = {str(v): sum(bool(x != y and y != ".")
                          for x,y in zip(seq,v.flat_pattern))
              for v in voclist}
        vnearby_names = [v.name.strip() for v in voclist if dv[str(v)]<2]
        fpatt = xpatt[seq]/n_sequences if n_sequences else 0
        yield "%s %6d %6.2f%% %s" % (rseq,xpatt[seq],100*fpatt,
                             ", ".join(vnearby_names))

def write_counts_file(filename,sitelist,master,voclist,cpatt,xpatt,n_sequences):
    '''write both the counts.out and x-count.out summary tables'''
    if not filename:
        return
    with open(filename,'w') as fout:
        for line in lineage_counts(sitelist,master,voclist,cpatt,n_sequences):
            print(line,file=fout)

    xfilename = covid.filename_prepend('x-',filename)
    v.vprint("Unmatched sequence patterns in file:",xfilename)
    with open(xfilename,'w') as fout:
        for line in missing_patterns_with_nearby(sitelist,master,voclist,xpatt,n_sequences):
            print(line,file=fout)


def main(args):
    ''' embers main '''
    v.vprint(args)

    ## read sequences, filter by pattern and by padded dates
    seqs = covid.read_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)
    seqs = covid.filter_seqs_by_pattern(seqs,args,keepfirst=False)
    seqs = emu.filter_seqs_by_padded_dates(seqs,args)
    if not args.keepx:
        seqs = (s for s in seqs if "X" not in s.seq)
    seqs = sequtil.checkseqlengths(seqs)

    if args.colormut:
        svar = SpikeVariants.from_colormut(args.colormut,refseq=first.seq)
    else:
        warnings.warn("Default color-mutation list may be out of date")
        svar = SpikeVariants.default(refseq=first.seq)

    OTHERNAME  = svar.OTHERNAME
    OTHERCOLOR = svar.OTHERCOLOR

    master = svar.master
    sitelist = svar.ssites()
    voclist = svar.vocs

    for voc in voclist:
        voc.flat_pattern = svar.flatpattern(voc)
        v.vvprint(f"{voc} {voc.name}: {voc.exact} {voc.flat_pattern}")

    mutants = [v.flat_pattern for v in voclist]
    patterns = mutants + [OTHERNAME]
    dups = check_dups(patterns)
    if dups:
        for dup in dups:
            for voc in voclist:
                if voc.flat_pattern == dup:
                    print("dup=",dup,"voc=",voc,voc.name)
        raise RuntimeError(f"Duplicated patterns {dups}")

    colors = [voc.color for voc in voclist] + [OTHERCOLOR]
    mcolors = dict(zip(patterns,colors))
    dups = check_dups(colors)
    if dups:
        v.vprint("Duplicated colors:",dups)

    namelist = [voc.name for voc in voclist] + [OTHERNAME]
    maxnamelen = max(len(n) for n in namelist)
    namefmt = "%%-%ds" % maxnamelen
    for voc in voclist:
        voc.name = namefmt % voc.name
    mnames =  {m: namefmt%n for m,n in zip(patterns,namelist)}
    dups = check_dups(namelist)
    if dups:
        v.vprint("Duplicated names:",dups)

    def relpattern(mut):
        if mut == OTHERNAME:
            return "." * len(master)
        return relativepattern(master,mut,dittochar='_')

    mrelpatt = {p: relpattern(p) for p in patterns}

    for p in patterns:
        v.vprint(mnames[p],mrelpatt[p],mcolors[p])

    if args.legend == 0:
        fullnames = None
    if args.legend == 1:
        fullnames = mnames
    if args.legend == 2:
        fullnames = {p: " ".join([mnames[p],mrelpatt[p]]) for p in patterns}

    svar.checkmaster(first.seq) ## ensure master agrees with first seqlist

    seqlist = list(seqs)
    n_sequences = len(seqlist)-1  ## -1 not to count the reference sequence

    if not args.keepx:
        seqlist = [s for s in seqlist if "X" not in s.seq]
        v.vprint("Removed",n_sequences+1-len(seqlist),"sequences with X")
        n_sequences = len(seqlist)-1

    ## How many of each sequence
    c = Counter(s.seq for s in seqlist[1:])

    ## How many of each mutant
    cpatt = Counter()
    xpatt = Counter()
    for seq in c:

        vocmatch = svar.vocmatch(seq)

        for voc in vocmatch[:1]:
            cpatt[voc.flat_pattern] += c[seq]

        ## Ideally just one match, if zero or more than one, then...
        if len(vocmatch)==0:
            sseq = svar.shorten(seq)
            xpatt[sseq] += c[seq]
        elif len(vocmatch)>1:
            if any(voc.name != vocmatch[0].name for voc in vocmatch):
                warn_msg = f"\n{svar.shorten(seq)} matches\n"
                warn_msg += " and\n".join(f"{relpattern(v.flat_pattern)} {v.name} {v}"
                                          for v in vocmatch)
                warn_msg += f"\n{svar.master} Master"
                v.vprint_only(10,'overlap:',warn_msg)

    v.vprint_only_summary('overlap:')
    v.vprint("Unmatched sequences:",sum(xpatt.values()))

    write_counts_file(args.ctable,sitelist,master,voclist,cpatt,xpatt,n_sequences)

    date_counter = {m: Counter() for m in patterns}
    for s in seqlist[1:]:
        if not args.keepx and "X" in s.seq:
            raise RuntimeError("X's should have already been filtered out")

        seqdate = emu.date_from_seqname(s.name)
        if not seqdate:
            continue

        vocmatch = svar.vocmatch(s.seq)
        for voc in vocmatch[:1]:
            ## if multiple matches, only count the first
            ## if that's an error, then you'll see a warning from earlier in the computation
            date_counter[voc.flat_pattern][seqdate] += 1
        if not vocmatch:
            date_counter[OTHERNAME][seqdate] += 1

    nmatches = sum(sum(date_counter[p].values()) for p in patterns)
    v.vprint("matched sequences:",nmatches)
    if nmatches==0:
        raise RuntimeError("No sequences for: " + " ".join(args.filterbyname))

    for p in date_counter:
        v.vprint(p,sum(date_counter[p].values()))

    ## Onset times for each mutant
    onsets=dict()
    if args.onsets:
        onsets.update( {m: min(date_counter[m]) for m in mutants if date_counter[m]} )

    ord_range, ord_plot_range = emu.get_ord_daterange(date_counter,args.dates)
    cum_counts = emu.get_cumulative_counts(date_counter,ord_range,
                                           daysperweek=args.daily)

    emu.make_emberstyle_plots(args,None,cum_counts,fullnames,mcolors,ord_range[0],
                              ordplotrange = ord_plot_range,
                              title=": ".join([covid.get_title(args),
                                               f"{n_sequences} sequences"]),
                              onsets=onsets)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    main(_args)
