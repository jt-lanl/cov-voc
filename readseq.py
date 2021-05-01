'''
Library of routines for reading sequence files (mase, fasta, tbl, seq)
And for writing them
And for dealing with them when the files are gzip'd
Type "seq" is raw sequences (no names)

Recommended usage to read a seqeuence file in fasta format:
   seqlist = read_seqfile(filename,type="fasta")
but if the filename ends in a "standard" extension, then the "type=" is not needed
   seqlist = read_seqfile("file.fasta")
'''

## read_seqfile() encapsulates all of these functions
## write_seqfile() for writing

import sys
import os.path
import os
import re
from gzip import open as gzip_open

from seqsample import SequenceSample
#class SequenceSample:
#    def __init__(self,name,seq):
#        self.name = name
#        self.seq = seq

def open_fcn(gzip=False):
    if gzip:
        return lambda filename: gzip_open(filename,"rt")
    else:
        return lambda filename: open(filename,"r")

def rd_rawseq(fp):
    for line in fp:
        line = line.strip()
        if line:
            yield SequenceSample('',line.strip())
        
def rd_mase(fp):
    re_comment = re.compile('^;.*')
    cur_name=''
    cur_seq=''

    for line in fp:
        line = line.strip()
        if cur_name:
            if re_comment.match(line):
                yield SequenceSample(cur_name,cur_seq)
                cur_name=''
                cur_seq=''
            else:
                cur_seq += line
        else:
            line = re_comment.sub('',line)
            if not line:
                continue
            cur_name=line
    if cur_name:
        yield SequenceSample(cur_name,cur_seq)

def rd_fasta(fp):
    re_seqname = re.compile('^>')

    cur_name=''
    cur_seq=''
    
    for line in fp:
        line = line.strip()
        if cur_name:
            if re_seqname.match(line):
                yield SequenceSample(cur_name,cur_seq)
                cur_name = re_seqname.sub('',line)
                cur_seq=''
            else:
                cur_seq += line
        else:
            if re_seqname.match(line):
                cur_name = re_seqname.sub('',line)
                cur_seq = ''
            else:
                raise RuntimeError("Seq name not yet defined")

    if cur_name:
        ## last sequence in the file
        yield SequenceSample(cur_name,cur_seq)

def read_fasta(filename,gzip=False):
    ## kept for backward compatibaility
    ## better to use: read_seqfile(filename,filetype="fasta")
    opfcn = open_fcn(gzip)
    with opfcn(filename) as fp:
        seqlist = list( rd_fasta(fp) )
    return seqlist
    
def rd_tbl(fp):
    re_comment    = re.compile('^\#.*')

    for line in fp:
        line = re_comment.sub('',line)
        line = line.strip()
        if not line:
            continue
        try:
            name,seq = line.split(None,1)
            yield SequenceSample(name,seq)
        except ValueError:
            print("line=[",line,"]")
            raise RuntimeError("Invalid line in tbl file")


FILEFUNCS = {
    "mase"  : rd_mase,
    "fasta" : rd_fasta,
    "fa"    : rd_fasta,
    "fst"   : rd_fasta,
    "fna"   : rd_fasta,
    "faa"   : rd_fasta,
    "tbl"   : rd_tbl,
    "table" : rd_tbl,
    "seq"   : rd_rawseq,
}
FILETYPES = list(FILEFUNCS)

def auto_filetype(filename,filetypes=FILETYPES):
    if filename.endswith(".gz"):
        gzip = True
        filename = re.sub(".gz$","",filename)
    else:
        gzip = False
        
    for t in filetypes:
        if filename.lower().endswith("."+t):
            return t,gzip

    _,ext = os.path.splitext(filename)
    raise RuntimeWarning(f"Filename extension [{ext}] "
                         f"not among supported filetypes: "
                         f"{filetypes}")
    return "",gzip

def rd_seqfile(filename,filetype="auto"):
    '''
    returns a GENERATOR of SequenceSample's
    filename=string (or libpath.Path): specifies sequence filename
    filetype=string: to specify file format: fasta,tbl, etc.
    '''
    ## potential advantage of this generator approach is that we could
    ## read arbitrarily large files, but only save into memory those
    ## sequences we need  (ie, filter as we read)
    
    filename = os.fspath(filename)

    filetype = filetype.lower()
    if filetype not in ["auto"] + FILETYPES:
        raise RuntimeError("filetype="+filetype+" not supported")

    atype,gzip = auto_filetype(filename)

    if filetype != "auto" and atype and FILEFUNCS[atype] != FILEFUNCS[filetype]:
        raise RuntimeWarning("filetype="+filetype+
                             " but file appears of filetype "+atype)

    if filetype == "auto":
        filetype = atype
    if filetype not in FILETYPES:
        raise RuntimeError("Unknown filetype ["+filetype+
                           "] of sequence file ["+filename+"]")

    ## having gone through all that to determine what kind
    ## of file this is, now start reading it
    o_fcn = open_fcn(gzip)
    read_fcn = FILEFUNCS[filetype]
    with o_fcn(filename) as fp:
        yield from read_fcn(fp)
            
def filter_seqs(seqs,
                rmdash=False,toupper=True,badchar=None,
                pattern=None,xpattern=None,maxseqs=None,
                rmdup=False,revseq=False):

    '''
    filters a GENERATOR of SequenceSample's, and returns a GENERATOR
    seqgen = generator (eg, output of rd_seqfile(filename) )
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    pattern=str: only keep seqeunces whose name matches this pattern
    xpattern=str: only keep seqeunces whose name does NOT match this pattern
    revseq=True/False: reverse the s.seq's
    maxseqs=int: only keep this many sequences
    '''

    ## all these different filters could be separate functions
    ## eg, instead of maxseqs keyword, could call
    ## seqs = itertools.islice(seqs,maxseqs) if maxseqs else seqs

    re_dash = re.compile('-')
    re_badchar = re.compile('[\#\$\*Xx]')
    if pattern:
        re_pattern = re.compile(pattern)
    if xpattern:
        re_xpattern = re.compile(xpattern)
    if rmdup:
        ## nb, if there are not a lot of repeats then the uniq set can become
        ## almost as large as the data file being read.  one way to save some
        ## space, if the sequences are long (longer than 32 bytes), is to use
        ## md5(s.seq).hexdigest() instead of s.seq as the comparative string
        ## [but then again, maybe this hashing already happens, and besides
        ## identical copies of strings may point to the same string??? or not.]
        uniq = set() 

    nseqs=0
    for s in seqs:
        if pattern and not re_pattern.search(s.name):
            continue
        if xpattern and re_xpattern.search(s.name):
            continue
        if rmdash:
            s.seq = re_dash.sub('',s.seq)
        if toupper:
            s.seq = s.seq.upper()
        if badchar:
            s.seq = re_badchar.sub(badchar,s.seq)
        if revseq:
            s.seq = s.seq[::-1]
        if rmdup:
            if s.seq in uniq:
                continue
            else:
                uniq.add(s.seq)
        yield s
        nseqs += 1
        if maxseqs and nseqs >= maxseqs:
            break

def read_seqfile(filename,filetype="auto",**kw):
    '''
    returns a list of SequenceSample's, with names and sequences for each one
    filename=string or pathlib.Path
    filetype=string: to specify file format: fasta,tbl, etc.
    **kw for filter_seqs:
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    pattern=str: only keep seqeunces whose name matches this pattern
    xpattern=str: only keep seqeunces whose name does NOT match this pattern
    '''
    seqs = rd_seqfile(filename,filetype=filetype)
    seqs = filter_seqs(seqs,**kw)
    return list(seqs)
    
def read_seqfile_old(filename,filetype="auto",rmdash=False,toupper=True,badchar=None):
    '''
    returns a list of SequenceSample's, with names and sequences for each one
    filetype=string: to specify file format: fasta,tbl, etc.
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    '''

    filename = os.fspath(filename)
    
    atype,gzip = auto_filetype(filename)

    if filetype != "auto" and filetype not in FILETYPES:
        raise RuntimeError("filetype="+filetype+" not supported")
    if filetype != "auto" and atype and FILEFUNCS[atype] != FILEFUNCS[filetype]:
        raise RuntimeWarning("filetype="+filetype+" but file appears of filetype "+afiletype)

    if filetype == "auto":
        filetype = atype
    if filetype in FILETYPES:
        o_fcn = open_fcn(gzip)
        read_fcn = FILEFUNCS[filetype]
        with o_fcn(filename) as fp:
            seqlist = list( read_fcn(fp) )
    else:
        raise RuntimeError("Unknown filetype ["+filetype+
                           "] of sequence file ["+filename+"]")

    if rmdash:
        re_dash = re.compile('-')
        for s in seqlist:
            s.seq = re_dash.sub('',s.seq)

    if toupper:
        for s in seqlist:
            s.seq = s.seq.upper()

    if badchar:
        re_badchar = re.compile('[\#\$\*Xx]')
        for s in seqlist:
            s.seq = re_badchar.sub(badchar,s.seq)

    return seqlist


################### WRITE SEQUENCE FILES



def columnize(str,col=70):
    strarr = []
    while len(str)>=col:
        strarr.append( str[:col] )
        str = str[col:]
    if str:
        strarr.append( str )
    return strarr

def write_fasta(filename,seq_samples):
    with open(filename,"w") as fout:
        for s in seq_samples:
            fout.write(">" + s.name + "\n")
            for str in columnize(s.seq):
                ## better to break these into 70 column units
                fout.write(str)
                fout.write("\n")

def write_mase(filename,seq_samples):
    with open(filename,"w") as fout:
        for s in seq_samples:
            fout.write(";\n%s\n" % s.name)
            for str in columnize(s.seq):
                fout.write(str + "\n")

def write_tbl(filename,seq_samples):
    with open(filename,"w") as fout:
        for s in seq_samples:
            sname = s.name.strip()
            sname = re.sub(" ","_",s.name) ## no spaces in name
            fout.write("%s\t%s\n" % (sname,s.seq))

def write_rawseq(filename,seq_samples):
    with open(filename,"w") as fout:
        for s in seq_samples:
            fout.write("%s\n" % s.seq)

W_FILEFUNCS = {
    #"mase"  : write_mase,
    "fasta" : write_fasta,
    "fst"   : write_fasta,
    "fa"    : write_fasta,
    "fna"   : write_fasta,
    "faa"   : write_fasta,
    "mase"  : write_mase,
    "tbl"   : write_tbl,
    "table" : write_tbl,
    "seq"   : write_rawseq,
}
            
            
W_FILETYPES = list(W_FILEFUNCS)
    
def write_seqfile(filename,seq_samples,filetype="auto"):

    filename = os.fspath(filename) ## enables pathlib input
    
    atype,gzip = auto_filetype(filename,filetypes=W_FILETYPES)

    if filetype != "auto" and filetype not in W_FILETYPES:
        raise RuntimeError("filetype="+filetype+" not supported")
    if filetype != "auto" and atype and atype != filetype:
        raise RuntimeError("filetype="+filetype+" but file appears of filetype "+atype)
    if gzip:
        raise RuntimeError("gzip not currently supported for writing files")

    if filetype == "auto":
        filetype = atype
    if filetype in W_FILETYPES:
        write_fcn = W_FILEFUNCS[filetype]
        write_fcn(filename,seq_samples)
    else:
        raise RuntimeError("Unknown filetype ["+filetype+
                           "] of sequence file ["+filename+"]")


if __name__ == "__main__":
    
    import argparse
    argparser = argparse.ArgumentParser()
    paa = argparser.add_argument
    paa("filename",
        help="Name of input file")
    paa("--filetype","-t",default="auto",
        help="Type of input file (mase, fasta, tbl, etc)")
    paa("--rmdash",action="store_true",
        help="Remove dashes from sequences")
    paa("--output","-o",
        help="Write seqeunce to file")
    paa("--ofiletype",default="auto",
        help="Type of output file")
    
    args = argparser.parse_args()
    
    seqlist = read_seqfile(args.filename,filetype=args.filetype,rmdash=args.rmdash)
    for s in seqlist[:5]:
        print(s.name,s.seq[:10],"...")

    if args.output:
        write_seqfile(args.output,seqlist,filetype=args.ofiletype)
        
