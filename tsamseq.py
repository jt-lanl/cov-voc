'''Time-Sampled Sequences
   Produce a fasta file with specified number of sequences per (day/week/month) sampled from the larger database of sequences
'''

## y'know, this could probably be an option to fixfasta rather than
## a whole separate proram

from pathlib import Path
import itertools as it
import argparse
import random

import verbose as v
import sequtil
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    covid.corona_args(argparser)
    paa = argparser.add_argument_group('Time-Sampling Options').add_argument
    paa("-K",type=int,default=30,
        help="number of sequences per D days")
    paa("-D",type=int,default=30,
        help="number of days per sampling period")
    paa("--output","-o",type=Path,
        help="file for writing sampled sequences")
    paa("--seed",type=int,default=17,
        help="random number seed")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    args = covid.corona_fixargs(args)
    return args

def by_date(s):
    '''given sequence s, return date'''
    date = covid.get_date(s.name)
    return date.toordinal() if date else None

def add_date_attribute(seqs):
    for s in seqs:
        s.date = by_date(s)
        if s.date:
            yield s

def _main(args):
    '''main'''
    v.vprint(args)

    random.seed(args.seed)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = covid.get_first_item(seqs,keepfirst=False)
    seqs = add_date_attribute(seqs)
    seqs = list(seqs)
    
    v.vprint('seqs:',len(seqs))
    
    seqs = sorted(seqs, key=lambda s: s.date)

    firstdate = min(s.date for s in seqs)
    lastdate = max(s.date for s in seqs)

    ndx = 0
    ndxlist = [0]
    for curdate in range(firstdate,lastdate+1):
        while ndx < len(seqs) and seqs[ndx].date <= curdate:
            ndx += 1
        v.vvprint(curdate-firstdate,ndx-1)
        ndxlist.append(ndx-1)

    v.vvprint('ndxlist:',ndxlist)

    def ndx_bydate(date):
        if date > lastdate:
            return max(ndxlist)
        try:
            ndx = ndxlist[date-firstdate]
        except IndexError:
            v.print('Index:',date-firstdate,'vs',lastdate-firstdate,'vs',len(ndxlist))
        return ndx

    tsamseqs = [first]
    for curdate in range(firstdate,lastdate+1,args.D):
        ndxlo = ndx_bydate(curdate)
        ndxhi = ndx_bydate(curdate+args.D)
        if ndxhi+1-ndxlo < args.K:
            ## If not enough samples, then double-dip!
            ndxsamples = random.choices(range(ndxlo,ndxhi+1),k=args.K)
        else:
            ndxsamples = random.sample(range(ndxlo,ndxhi+1),k=args.K)
        v.vvprint(curdate-firstdate,
                  'ndxrange:',ndxlo,ndxhi,ndxhi-ndxlo+1,
                  ndxsamples)

        for ndx in sorted(ndxsamples):
            tsamseqs.append(seqs[ndx])

    if args.output:
        v.vprint(f'Writing {len(tsamseqs)} sequences to {args.output}')
        sequtil.write_seqfile(args.output,tsamseqs)
        

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
