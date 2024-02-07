'''Parallel python module'''

import subprocess
import argparse
import verbose as v

def _getargs_xtras():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument_group('Program Options').add_argument
    paa("pmodule",
        help="python module to be run through GNU parallel")
    paa("--jobs","-j",
        help="Number of jobs (if not specified, then as many as can fit)")
    paa("--input","-i",
        help="input file")
    paa("--output","-o",
        help="output file")
    paa("--ppmverbose","-V",action="count",default=0,
        help="verbose for ppm (not for proram module)")
    paa("--ppmexecute","-X",action="store_true",
        help="Actually run the command")

    args,xtras = argparser.parse_known_args()
    return args,xtras

def _main(args,xtras):
    '''main'''
    v.vvprint(args)
    v.vvprint("XTRA:",xtras)

    jobsopt="--jobs={args.jobs}" if args.jobs else ""

    ## To-Do: if input endsiwth ".tbl" then xtras.append("-i /dev/stdin.tbl")
    ## To-DO: if input endswith ".xz" then we should prepend "unxz input |"
    input_cmd = ""
    if args.input:
        xtras.extend(["-i","/dev/stdin"])
        cat = "cat"
        if args.input.endswith('.xz'):
            cat = "xz --decompress --stdout"
        if args.input.endswith('.gz'):
            cat = "gunzip --stdout"
        input_cmd = f"{cat} {args.input} |"

    output_cmd = ""
    if args.output:
        xtras.extend(["-o","/dev/stdout"])
        if args.output not in ["-", "/dev/stdout"]:
            output_cmd = f"> {args.output}"

    xtras = [f"'{xtra}'" for xtra in xtras]
    xtra_cmd = " ".join(xtras)
    cmdline = (f"{input_cmd} "
    f"parallel {jobsopt} -k --header '>[^>]*' --recstart '>' --blocksize 400M --pipe "
               f"python -m {args.pmodule} {xtra_cmd} "
               "--jobno={#} "  ## not an f-string here!
               f"{output_cmd}")

    v.print(cmdline)

    if args.ppmexecute:
        v.vprint("Running...")
        retval = subprocess.call(cmdline, shell=True)
        v.vprint('returned value:', retval)
    else:
        v.print('Use -X option to actually execute this command')

if __name__ == "__main__":

    _args,_xtras = _getargs_xtras()
    v.verbosity(_args.ppmverbose)
    _main(_args,_xtras)
