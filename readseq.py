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

import re
import pickle
import warnings

from xopen import xopen
from seqsample import SequenceSample

DEV_NAMES = "- /dev/null /dev/stdin /dev/stderr /dev/stdout".split()


################### READ SEQUENCE FILES


def rd_names(fp):
    '''read names file; only names, no sequences'''
    re_seqname_prefix = re.compile(r'^>\s*')
    for line in fp:
        line = re_seqname_prefix.sub('',line)
        yield SequenceSample(line.strip(),"")

def rd_rawseq(fp):
    '''read raw sequences file; only sequences, no names'''
    for line in fp:
        line = line.strip()
        if line:
            yield SequenceSample('',line.strip())

def rd_mase(fp):
    '''read mase file'''
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
    '''read fasta file'''
    re_seqname = re.compile('^>')

    cur_name=''
    cur_seq=''

    for line in fp:
        line = line.strip()
        if cur_name:
            if re_seqname.match(line):
                yield SequenceSample(cur_name,cur_seq)
                cur_name = re_seqname.sub('',line).strip()
                cur_seq=''
            else:
                cur_seq += line
        else:
            if re_seqname.match(line):
                cur_name = re_seqname.sub('',line).strip()
                cur_seq = ''
            else:
                raise RuntimeError("Seq name not yet defined")

    if cur_name:
        ## last sequence in the file
        yield SequenceSample(cur_name,cur_seq)

def rd_tbl(fp):
    '''read tbl file; name and sequence on same line'''
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
            warnings.warn(f'Invalid line=[{line}] in tbl file')

## Note: pickle files are experimental

def rd_pickle(fp):
    '''read pickle file (contains sequence list as a single object)'''
    return pickle.load(fp)

def rd_incr_pickle(fp):
    '''read incremental pickle file'''
    while True:
        try:
            yield pickle.load(fp)
        except EOFError:
            return


R_FILEFUNCS = {
    "fasta" : rd_fasta,
    "fa"    : rd_fasta,
    "fst"   : rd_fasta,
    "fna"   : rd_fasta,
    "faa"   : rd_fasta,
    "tbl"   : rd_tbl,
    "table" : rd_tbl,
    "mase"  : rd_mase,
    "seq"   : rd_rawseq,
    "pkl"   : rd_pickle,
    "ipkl"  : rd_incr_pickle,
    "nm"    : rd_names,
}

def auto_filetype(filename,filetypelist):
    '''determind file type from filename extension(s)'''

    for ftype in filetypelist:
        if re.search(r'\.'+ftype+r'\b',filename):
            return ftype

    return 'fasta'

def get_fileinfo(filepath,filetype,filetypelist):
    '''return: filename,filetype,is_pklfile'''

    ## check that keyword argument filetype is valid
    filetype=filetype.lower()
    if filetype not in ["auto"] + list(filetypelist):
        raise RuntimeError(f"filetype={filetype} not supported")

    filename = str(filepath)
    atype = auto_filetype(filename.lower(),filetypelist)
    if filetype not in ["auto",atype]:
        raise RuntimeError(f'Specified filetype={filetype} but filename={filename} '
                           f'suggests file type is {atype}')

    if filetype == "auto":
        filetype = atype

    is_pklfile = bool(filetype in ["ipkl","pkl"])

    for dev in DEV_NAMES:
        if filename.startswith(dev):
            filename=dev

    return filename,filetype,is_pklfile


def rd_seqfile(filename,filetype="auto"):
    '''
    returns a GENERATOR of SequenceSample's
    filename=string (or libpath.Path): specifies sequence filename
    filetype=string: to specify file format: fasta,tbl, etc.
    '''
    ## potential advantage of this generator approach is that we could
    ## read arbitrarily large files, but only save into memory those
    ## sequences we need  (ie, filter as we read)

    filename,filetype,is_pklfile = get_fileinfo(filename,
                                                filetype,
                                                list(R_FILEFUNCS))

    read_fcn = R_FILEFUNCS[filetype]
    rwstr = 'rb' if is_pklfile else 'r'
    with xopen(filename,rwstr) as fp:
        yield from read_fcn(fp)

def filter_seqs(seqs,
                rmdash=False,toupper=True,badchar=None,
                maxseqs=None,rmdup=False,revseq=False):

    '''
    filters a GENERATOR of SequenceSample's, and returns a GENERATOR
    seqgen = generator (eg, output of rd_seqfile(filename) )
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    revseq=True/False: reverse the s.seq's
    maxseqs=int: only keep this many sequences
    '''

    ## all these different filters could be separate functions
    ## eg, instead of maxseqs keyword, could call
    ## seqs = itertools.islice(seqs,maxseqs) if maxseqs else seqs

    re_dash = re.compile(r'-')

    ## old-style badchar as regexp
    #re_badchar = re.compile(r'[\$\#\*Xx]')
    #re_badlastchar = re.compile(r'[\#Xx]')

    ## new badchar based on translation (faster)
    ## Always keep '*'
    ## Translate '$' into '*'
    ## Convert '#' or 'x' into badchar (badchar='X' is typical)
    badchar_dict = dict((c,badchar) for c in r'#Xx')
    badchar_dict['$'] = '*'
    tr_badchar = "".maketrans( badchar_dict )

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

def read_seqfile(filename,filetype="auto",nofilter=False,**kw):
    '''
    returns a list of SequenceSample's, with names and sequences for each one
    filename=string or pathlib.Path
    filetype=string: to specify file format: fasta,tbl, etc.
    **kw for filter_seqs:
    rmdash=True/False: to remove dashes from every sequence (loses alignment)
    toupper=True/False: convert all sequence characters to upper case
    badchar=char: replace "bad" characters (currently: #$*X) with specified character
    '''
    seqs = rd_seqfile(filename,filetype=filetype)
    if not nofilter:
        seqs = filter_seqs(seqs,**kw)
    return seqs


################### WRITE SEQUENCE FILES


def columnize(longstring,col=0):
    '''break up a long string into a sequence of strings with col columns or less'''
    if col == 0:
        ## default is not to columnize
        yield longstring
    else:
        while longstring:
            yield longstring[:col]
            longstring = longstring[col:]

def wt_fasta(fp,s):
    '''write single SequenceSample(name,seq) in fasta format'''
    print(f'>{s.name}',file=fp)
    for line in columnize(s.seq):
        print(line,file=fp)

def wt_mase(fp,s):
    '''write single SequenceSample(name,seq) in mase format'''
    print(f';\n{s.name}',file=fp)
    for line in columnize(s.seq):
        print(line,file=fp)

def wt_tbl(fp,s):
    '''write single SequenceSample(name,seq) in tbl (table) format'''
    sname = re.sub(" ","_",s.name.strip())
    print(f'{sname}\t{s.seq}',file=fp)

def wt_rawseq(fp,s):
    '''write single sequence in raw sequence format'''
    print(s.seq,file=fp)

def wt_incr_pickle(fp,s):
    '''write single SequenceSample(name,seq) in incremental pickle format'''
    pickle.dump(s,fp)

def wt_names(fp,s):
    '''write single name in .nm format'''
    print(s.name,file=fp)

W_FILEFUNCS = {
    "fasta" : wt_fasta,
    "fst"   : wt_fasta,
    "fa"    : wt_fasta,
    "fna"   : wt_fasta,
    "faa"   : wt_fasta,
    "mase"  : wt_mase,
    "tbl"   : wt_tbl,
    "table" : wt_tbl,
    "seq"   : wt_rawseq,
    "ipkl"  : wt_incr_pickle,
    "nm"    : wt_names,
}


W_FILETYPES = list(W_FILEFUNCS)

def write_seqfile(filepath,seq_samples,filetype="auto"):
    '''write sequence file in format determined by filename and/or filetype keyword'''

    filename,filetype,is_pkl = get_fileinfo(filepath,
                                            filetype,
                                            list(W_FILEFUNCS))
    rwstr = 'wb' if is_pkl else 'w'
    with xopen(filename,rwstr) as fout:
        if filetype == 'pkl':
            pickle.dump(list(seq_samples),fout)
        elif filetype in W_FILEFUNCS:
            wt_fcn = W_FILEFUNCS[filetype]
            for s in seq_samples:
                wt_fcn(fout,s)
        else:
            raise RuntimeError(f'Invalid filetype={filetype}')

if __name__ == "__main__":

    import argparse
    argparser = argparse.ArgumentParser()
    paa = argparser.add_argument
    paa("--input","-i",
        help="Name of input file")
    paa("--filetype","-t",default="auto",
        help="Type of input file (mase, fasta, tbl, etc)")
    paa("--rmdash",action="store_true",
        help="Remove dashes from sequences")
    paa("--output","-o",
        help="Name of output file")
    paa("--otype",default="auto",
        help="Type of output file")

    _args = argparser.parse_args()

    _seqs = read_seqfile(_args.input,
                         filetype=_args.filetype,
                         rmdash=_args.rmdash)

    if _args.output:
        write_seqfile(_args.output,_seqs,
                      filetype=_args.otype)
