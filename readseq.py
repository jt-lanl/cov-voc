'''
Library of routines for reading sequence files (mase, fasta, tbl, seq; also pkl)
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
import pickle
from gzip import open as gzip_open
from lzma import open as lzma_open

import warnings

from seqsample import SequenceSample

DEV_NAMES = "- /dev/null /dev/stdin /dev/stderr /dev/stdout".split()

def xopen(filename,rw,gz=False,xz=False,binaryfile=False):
    ''' equivalent of open that works w/ and w/o gzip '''

    if filename == "-":
        return sys.stdout if rw=="w" else sys.stdin

    rwstr = rw + "b" if binaryfile else rw
    if (gz or xz) and not binaryfile:
        rwstr += "t"

    if gz: return gzip_open(filename,rwstr)
    if xz: return lzma_open(filename,rwstr)

    return open(filename,rwstr)

def rd_names(fp):
    '''only names, no sequences'''
    re_seqname_prefix = re.compile(r'^>\s*')
    for line in fp:
        line = re_seqname_prefix.sub('',line)
        yield SequenceSample(line.strip(),"")

def rd_pickle(fp):
    '''pickled file contains sequence list as a single object'''
    return pickle.load(fp)

def rd_incr_pickle(fp):
    while True:
        try:
            yield pickle.load(fp)
        except EOFError:
            return

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

def read_fasta(filename,gz=False,xz=False):
    ## kept for backward compatibaility
    ## better to use: read_seqfile(filename,filetype="fasta")
    with xopen(filename,'r',gz=gz,xz=xz) as fp:
        seqlist = list( rd_fasta(fp) )
    return seqlist

def rd_tbl(fp):
    ## table format is two whitespace-separated strings on each line (Hyejin says: tab separated)
    ## first string it name, second string is sequence
    ## Spaces in the name string are not allowed ... but sometimes creep in
    ## To deal with that, take name as first token and sequence as last token

    re_comment    = re.compile(r'^\#.*')

    for line in fp:
        line = re_comment.sub('',line)
        line = line.strip()
        if not line:
            continue
        try:
            tokens = line.split('\t') if '\t' in line else line.split()
            if len(tokens) != 2:
                warnings.warn(f"Invalid tbl line: {' '.join(tokens[0:-1])} {tokens[-1][:5]}...")
            name,seq = tokens[0],tokens[-1]
            yield SequenceSample(name.strip(),seq.strip())
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
    "pkl"   : rd_pickle,
    "ipkl"  : rd_incr_pickle,
    "nm"    : rd_names,
}
FILETYPES = list(FILEFUNCS)

def auto_filetype(filename,filetypes=None):
    '''determind file type from filename extension(s)'''
    filetypes = filetypes or FILETYPES
    gz = xz = False
    if filename.endswith(".gz"):
        gz = True
        filename = re.sub(".gz$","",filename)
    elif filename.endswith(".xz"):
        xz = True
        filename = re.sub(".xz$","",filename)

    for t in filetypes:
        if filename.lower().endswith("."+t):
            return t,gz,xz

    if filename in DEV_NAMES:
        ## default for stdin or stdout is fasta
        return 'fasta',gz,xz

    _,ext = os.path.splitext(filename)
    raise RuntimeWarning(f"Filename extension [{ext}] "
                         f"not among supported filetypes: "
                         f"{filetypes}")

def get_fileinfo(filename,filetype,filetypelist):
    '''return: filename,filetype,gz,xz,binaryfile'''

    ## check that keyword argument filetype is valid
    filetype = filetype.lower()
    if filetype not in ["auto"] + list(filetypelist):
        raise RuntimeError(f"filetype={filetype} not supported")

    filename = os.fspath(filename) ## if pathlib input, make a string
    atype,gz,xz = auto_filetype(filename,filetypes=filetypelist)

    if filetype != "auto" and atype and atype != filetype:
        raise RuntimeError(f'Specified filetype={filetype} but filename={filename} '
                           f'suggests file type is {atype}')

    if filetype == "auto":
        filetype = atype

    if filetype not in filetypelist:
        raise RuntimeError(f"Unknown filetype={filetype} for sequence filename={filename}")

    binaryfile = bool(filetype in ["ipkl","pkl"])

    for dev in DEV_NAMES:
        if filename.startswith(dev):
            filename=dev

    return filename,filetype,gz,xz,binaryfile


def rd_seqfile(filename,filetype="auto"):
    '''
    returns a GENERATOR of SequenceSample's
    filename=string (or libpath.Path): specifies sequence filename
    filetype=string: to specify file format: fasta,tbl, etc.
    '''
    ## potential advantage of this generator approach is that we could
    ## read arbitrarily large files, but only save into memory those
    ## sequences we need  (ie, filter as we read)

    (filename,filetype,
     gz,xz,binaryfile) = get_fileinfo(filename,
                                      filetype,
                                      FILETYPES)

    read_fcn = FILEFUNCS[filetype]
    with xopen(filename,'r',gz=gz,xz=xz,binaryfile=binaryfile) as fp:
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

    re_dash = re.compile(r'-')

    ## old-style badchar as regexp
    re_badchar = re.compile(r'[\$\#\*Xx]')
    re_badlastchar = re.compile(r'[\#Xx]')

    ## new badchar based on translation (faster)
    ## Always keep '*'
    ## Translate '$' into '*'
    ## Convert '#' or 'x' into badchar (badchar='X' is typical)
    badchar_dict = dict((c,badchar) for c in r'#Xx')
    badchar_dict['$'] = '*'
    tr_badchar = "".maketrans( badchar_dict )

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
        if len(s.seq):
            ## if s.seq is empty, don't bother with all this stuff
            if rmdash:
                s.seq = re_dash.sub('',s.seq)
            if toupper:
                s.seq = s.seq.upper()
            if badchar:
                #translate bad character into badchar
                s.seq = s.seq.translate(tr_badchar)
            if revseq:
                s.seq = s.seq[::-1]
            if rmdup:
                if s.seq in uniq:
                    continue
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
    return seqs
    #return list(seqs)

def read_seqfile_old(filename,filetype="auto",rmdash=False,toupper=True,badchar=None):
    '''
    returns a list of SequenceSample's, with names and sequences for each one
    filetype=string: to specify file format: fasta,tbl, etc.
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    '''

    filename = os.fspath(filename)

    atype,_,_ = auto_filetype(filename)

    if filetype != "auto" and filetype not in FILETYPES:
        raise RuntimeError("filetype="+filetype+" not supported")
    if filetype != "auto" and atype and FILEFUNCS[atype] != FILEFUNCS[filetype]:
        raise RuntimeWarning("filetype="+filetype+" but file appears of filetype "+afiletype)

    if filetype == "auto":
        filetype = atype
    if filetype in FILETYPES:
        read_fcn = FILEFUNCS[filetype]
        with xopen(filename,'r') as fp:
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
        re_badchar = re.compile(r'[\#\$\*Xx]')
        for s in seqlist:
            s.seq = re_badchar.sub(badchar,s.seq)

    return seqlist


################### WRITE SEQUENCE FILES


def columnize(str,col=0):
    '''break up a long string into a sequence of strings with col columns or less'''
    if col == 0:
        ## default is not to columnize
        yield str
    else:
        while len(str)>=col:
            yield str[:col]
            str = str[col:]
        if str:
            yield str

def write_fasta(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz) as fout:
        for s in seq_samples:
            fout.write(">" + str(s.name) + "\n")
            for seq in columnize(s.seq):
                fout.write(seq)
                fout.write("\n")

def write_mase(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz) as fout:
        for s in seq_samples:
            fout.write(";\n%s\n" % s.name)
            for str in columnize(s.seq):
                fout.write(str + "\n")

def write_tbl(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz) as fout:
        for s in seq_samples:
            sname = s.name.strip()
            sname = re.sub(" ","_",s.name) ## no spaces in name
            fout.write("%s\t%s\n" % (sname,s.seq))

def write_rawseq(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz) as fout:
        for s in seq_samples:
            fout.write("%s\n" % s.seq)

def write_pickle(filename,seq_samples,gz=False,xz=False):
    ## use list(seq_samples) in case it's a generator
    with xopen(filename,'w',gz=gz,xz=xz,binaryfile=True) as fout:
        pickle.dump(list(seq_samples),fout)

def write_incr_pickle(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz,binaryfile=True) as fout:
        for s in seq_samples:
            pickle.dump(s,fout)

def write_names(filename,seq_samples,gz=False,xz=False):
    with xopen(filename,'w',gz=gz,xz=xz,binaryfile=False) as fout:
        for s in seq_samples:
            fout.write("%s\n" % s.name)


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
    "pkl"   : write_pickle,
    "ipkl"  : write_incr_pickle,
    "nm"    : write_names,
}


W_FILETYPES = list(W_FILEFUNCS)

def write_seqfile(filename,seq_samples,filetype="auto"):

    (filename,filetype,
     gz,xz,binaryfile) = get_fileinfo(filename,
                                      filetype,
                                      W_FILETYPES)
    write_fcn = W_FILEFUNCS[filetype]
    write_fcn(filename,seq_samples,gz=gz,xz=xz)

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

    seqs = read_seqfile(args.filename,
                        filetype=args.filetype,
                        rmdash=args.rmdash)

    if args.output:
        write_seqfile(args.output,seqs,
                      filetype=args.ofiletype)
