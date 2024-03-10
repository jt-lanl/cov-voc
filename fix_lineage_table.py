'''automatically write regex's for (some of) the entries in the lineage table'''

import re
import argparse
import verbose as v
import covid

from lineagenotes import LineageNotes

DEFAULT_LINEAGE_NOTES_FILE="/home/jt/src/corona/data/lineage_notes.txt"

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--tablefile",
        help="lineage table file")
    paa("--notesfile",
        help="lineage_notes.txt")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    if not args.notesfile:
        args.notesfile = covid.find_seqfile(LineageNotes.default_file)
        v.vprint(f'Using notes file: {args.notesfile}')
    return args

def regexp_from_kidlist(kidlist):
    '''eg, EK EL EM EU FD DG FH => (E[KLMU])|(F[DGH])'''
    kidlist2 = [kid for kid in kidlist if len(kid)==2]
    kidxtras = [kid for kid in kidlist if len(kid)!=2]
    firstchars = sorted(set(kid[0] for kid in kidlist2))
    regexp=[]
    for kid in kidxtras:
        regexp.append(f'({kid})')
    for firstchar in firstchars:
        seconds = [kid[1] for kid in kidlist2 if kid[0]==firstchar]
        seconds = "".join(sorted(seconds))
        if len(seconds) > 1:
            regexp.append(f'({firstchar}[{seconds}])')
        else:
            regexp.append(f'({firstchar}{seconds})')
    return "|".join(regexp)

def get_regexp(parent,kidlist):
    '''return a regexp associated with parent and kidlist'''
    ## replace "." with "\."
    regexp = re.sub(r'\.',r'\.',parent)
    ## get regexp for kids
    regexp_kids = regexp_from_kidlist(kidlist)
    ## add kids to regexp
    if regexp_kids:
        regexp = f'(({regexp})|{regexp_kids})'
    ## append the magic "(\..*)?" to enable children of children
    regexp = f'{regexp}(\\..*)?'
    return regexp

Pango_WHO = {
    'Alpha':   'B.1.1.7',
    'Beta':    'B.1.351',
    'Gamma':   'B.1.1.28',
    'Mu':      'B.1.621',
    'Epsilon': 'B.1.42[97]',
    'Iota':    'B.1.526',
    'Delta':   'B.1.617.2',
}

def _main(args):
    '''main'''
    v.vprint(args)

    lin = LineageNotes.from_file(args.notesfile)
    v.vprint(lin.report_size())
    for bad in lin.inconsistencies(remove=True):
        v.print(bad)
    v.vprint(lin.report_size())

    with open(args.tablefile,'r') as fin:
        for line in fin:
            line = line.strip()
            if not line or line[0]=='#':
                print(line)
                continue
            tokens = line.split()
            if len(tokens) < 3:
                v.vprint('Invalid line:',line)
                continue
            color,name,_ = tokens[:3]
            parent = name.split('/')[0]
            parent = Pango_WHO.get(parent,parent)
            if not re.match(r'[A-Z]+(\.[0-9\.]*)?$',parent):
                print(line)
                continue

            xparent = lin.fullname.get(parent,parent)
            v.vprint(f'{parent}={xparent}')
            children=list()
            for pango in lin.lineages:
                if pango.startswith(xparent+"."):
                    children.append(pango)
                elif lin.fullname.get(pango,'').startswith(xparent+"."):
                    children.append(pango)
            kids = set()
            for pango in children:
                kid,_ = pango.split('.',1)
                if parent.startswith(kid):
                    continue
                kids.add(kid)
            print(color,name,get_regexp(parent,kids),sep='\t')


if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
