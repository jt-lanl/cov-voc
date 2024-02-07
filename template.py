'''Template doc string'''

import argparse
import verbose as v

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    generic_paa = argparser.add_argument
    generic_paa("--verbose","-v",action="count",default=0,
                help="verbose")
    paa = argparser.add_argument_group('Program Options').add_argument
    args = argparser.parse_args()
    return args

def _main(args):
    '''main'''
    v.vprint(args)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
