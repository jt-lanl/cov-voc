'''xopen is a replacement for open that handles compressed files'''

import sys
import gzip
import lzma
import zstdx

def xopen(filepath,mode):
    '''like builtins.open() but works with compressed files too'''
    filestr=str(filepath)

    if filestr == '-':
        filepath = sys.stdout if 'w' in mode else sys.stdin

    if filestr.endswith('.zst'):
        return zstdx.open(filepath,mode)

    if 'b' not in mode:
        ## if not binary, then text (gzip,lzma want you to be explicit)
        mode += 't'
    if filestr.endswith('.gz'):
        return gzip.open(filepath,mode)
    if filestr.endswith('.xz'):
        return lzma.open(filepath,mode)
    ## default is just the usual
    return open(filepath,mode)

if __name__ == "__main__":
    import argparse
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",help="input file")
    paa("--output","-o",help="output file")
    args = argparser.parse_args()

    with xopen(args.input,'r') as fin:
        with xopen(args.output,'w') as fout:
            for line in fin:
                print(line,file=fout)
