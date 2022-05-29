'''Template doc string'''

import argparse
from collections import Counter
from verbose import verbose as v

import sequtil
import covid

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    covid.corona_args(argparser)
    paa("--output","-o",
        help="write output file with fixed names")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def fixname(seqname):
    '''fix common errors in sequence names
    returns tuple: fixed name, and bad location'''

    ## Most common error is too many components to geographical location
    tokens = seqname.split('.')
    try:
        isl_token = tokens[5]
        if isl_token[:8] == "EPI_ISL_":
            return seqname,None
    except IndexError:
        pass
    for ndx,token in enumerate(tokens):
        if token[:8] == "EPI_ISL_":
            location = ".".join(tokens[1:ndx-1])
            if ndx > 5:
                new_location = ".".join(tokens[1:4]) + "_" + "_".join(tokens[4:ndx-1])
            else:
                new_location = ".".join(tokens[1:ndx-1] + [""]*(5-ndx))
            v.vvprint("location:",location,"->",new_location)
            v.vprint_only(1,f'bad location: {location} -> {new_location}')
            seqname = ".".join(tokens[:1] + [new_location] + tokens[ndx-1:])
            return seqname,location
    return None,None

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)
    if args.output:
        seqs = list(seqs)

    bad_locations = Counter()
    for s in seqs:
        new_name,old_location = fixname(s.name)
        if new_name:
            if new_name != s.name and old_location not in bad_locations:
                v.vprint('new name:',new_name)
            if old_location is not None:
                bad_locations[old_location] += 1
            s.name = new_name
        else:
            v.vprint("Bad name:",s.name)

    for bad_location,cnt in bad_locations.items():
        v.vprint(f"{cnt:6d} {bad_location}")

    v.print("Fixed names:",sum(bad_locations.values()))

    if args.output:
        sequtil.write_seqfile(args.output,seqs)
        v.print(f"Fixed names: {100*sum(bad_locations.values())/len(seqs):.2f}%")


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
