'''Template doc string'''

import re
import argparse
from collections import Counter
import verbose as v

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

DATE_REGEX = re.compile(r'\d\d\d\d-\d\d-\d\d')
BETA_REGEX = re.compile(r'B\.1\.351(\..*)?')
BETA_DATE = covid.date_fromiso('2020-06-01')

def fixlocation(seqname):
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

def fixlocations(seqs):
    '''apply fixlocation to full seqlist'''
    bad_locations = Counter()
    seqcount=0
    for s in seqs:
        seqcount += 1

        new_name,old_location = fixlocation(s.name)
        if new_name:
            if new_name != s.name and old_location not in bad_locations:
                v.vprint('new name:',new_name)
            if old_location is not None:
                bad_locations[old_location] += 1
            s.name = new_name
        else:
            v.vprint("Bad name:",s.name)

        yield s

    for bad_location,cnt in bad_locations.items():
        v.vprint(f"{cnt:6d} {bad_location}")

    v.print("Fixed names:",sum(bad_locations.values()))
    v.print(f"Fixed names: {100*sum(bad_locations.values())/seqcount:.2f}%")


def fixdate(seqname):
    '''fix common errors in date'''
    tokens = seqname.split('.',6)
    date_token = tokens[4]
    date = covid.date_fromiso(date_token)
    if date is None:
        v.vprint_only(5,'bad date:',date_token,seqname)
        return None

    lineage_token = tokens[6]
    if BETA_REGEX.match(lineage_token) and date < BETA_DATE:
        tokens[4] = "2021" + date_token[4:]
        return ".".join(tokens)

    return seqname

def fixdates(seqs):
    '''apply fixdate to each sequence in seqlist -- return new iterator'''
    seqcount=0
    badcount=0
    fixcount=0
    for s in seqs:
        seqcount += 1
        new_name = fixdate(s.name)
        if new_name is None:
            badcount += 1
            continue
        if new_name != s.name:
            fixcount += 1
            s.name = new_name
        yield s

    v.vprint(f'Fix Dates: {seqcount} total, '
             f'{fixcount} fixed, {badcount} removed')

def _main(args):
    '''main'''
    v.vprint(args)

    seqs = covid.read_filter_seqfile(args)
    first,seqs = sequtil.get_first_item(seqs,keepfirst=False)

    seqs = fixlocations(seqs)
    seqs = fixdates(seqs)

    if args.output:
        sequtil.write_seqfile(args.output,
                              [first] + list(seqs))


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
