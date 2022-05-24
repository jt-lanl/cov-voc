'''From a file of mstrings, produce output of the same file,
   but formatted so that the components of the mstrings are aligned
'''

import re
import itertools as it
import argparse

from verbose import verbose as v
import intlist
import mutant

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("mstringfile",
        help="file with mstrings")
    paa("--nolabels",action="store_true",
        help="only print mstrings, no text to left")
    paa("--compact",action="store_true",
        help="print a compact version of the aligned mutations")
    paa("--output","-o",
        help="output file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

def read_label_mstring(filename):
    '''
    read line of input of form "LABEL [MSTRING]
    and yield label,mstring tuple
    '''
    with open(filename) as fin:
        for fullline in fin:
            line = fullline.strip()
            line = re.sub(r'#.*','',line)
            if not line:
                continue
            label_mstring = re.match(r'(.*)\s*(\[.*\])',line)
            if not label_mstring:
                v.vprint(f'no match for: {fullline}')
            label = label_mstring[1]
            mstring = label_mstring[2]
            label= label.strip()
            yield label,mstring

def get_allsites(mutlist):
    '''
    return a list of all sites that appear in all mutations
    list is sorted, and duplicates are removed
    '''
    msites = [ssm.site for mut in mutlist for ssm in mut]
    msites = sorted(set(msites))
    return msites

def get_sub_site_ref(mutlist):
    '''
    get a dict of ref values for all the sites
    associated with substitutions and deletions
    '''
    return {ssm.site: ssm.ref
            for mut in mutlist
            for ssm in mut
            if ssm.ref != '+'}

def get_ins_site_len(mutlist):
    '''
    get a dict of maximum string length associated with each insertion site
    '''
    ins_site_len = dict()
    for mut in mutlist:
        for ssm in mut:
            if ssm.ref != '+':
                continue
            ins_site_len[ssm.site] = max( [ins_site_len.get(ssm.site,0), len(ssm.mut)] )
    return ins_site_len

def get_sssm_bysite(mut):
    '''
    return dict of ssm's keyed by site: only substitution and deletion ssm's
    '''
    return {ssm.site: ssm for ssm in mut if ssm.ref != "+"}

def get_issm_bysite(mut):
    '''
    return dict of ssm's keyed by site: only insertion ssm's
    '''
    return {ssm.site: ssm for ssm in mut if ssm.ref == "+"}

def emptysite(site):
    '''return empty string associated iwth AnnnB for site=nnn'''
    return " "*(2+len(str(site)))

class MStringFramework:
    '''
    framework for dealing with mstring's
    '''
    def __init__(self,mutlist):
        self.allsites = get_allsites(mutlist)
        self.ins_site_len = get_ins_site_len(mutlist)
        self.sub_site_ref = get_sub_site_ref(mutlist)

    def full_sitelist(self):
        '''
        create list of ints for sites, with extras for insertions
        yield list items one at a time
        '''
        for site in self.allsites:
            if site in self.sub_site_ref:
                yield site
            for _ in range(self.ins_site_len.get(site,0)):
                yield site


    def format_mstring(self,mut):
        '''
        format the mstring associated with mutation
        by adding extra spaces so that the all the mstrings will line up
        '''

        sssm_bysite = get_sssm_bysite(mut)
        issm_bysite = get_issm_bysite(mut)

        ## begin with a list of ssm strings
        ## later, will convert list to a single mstring using join
        mstring = []

        ## build up mstring for every potential site
        for site in self.allsites:
            if site in self.sub_site_ref:
                ssm = sssm_bysite.get(site,None)
                ssm_string = str(ssm) if ssm else emptysite(site)
                mstring.append( ssm_string )
            if site in self.ins_site_len:
                ## insertion site formatting is a little convoluted
                if site in issm_bysite:
                    sstr = str(issm_bysite[site])
                    if not issm_bysite[site].mut:
                        sstr += "x"
                    ## pad with tildes so comma is not removed
                    sstr += "~" * (self.ins_site_len[site] + 1 + len(str(site)) - len(sstr))
                else:
                    sstr = " " * (self.ins_site_len[site] + 1 + len(str(site)))
                mstring.append(sstr)

        ## joint ssm's in mstring list into a single mstring
        mstring = "[" + ",".join(mstring) + "]"
        ## remove extraneous commas
        mstring = re.sub(" ,","  ",mstring)
        ## remove tilde characters
        mstring = re.sub(r'([~]+),',r',\1',mstring)
        mstring = re.sub("~"," ",mstring)
        return mstring

    def compact_format_mstring(self,mut=None):
        '''
        return compactly-formatted mstring for input mutation
        '''

        mut = [] if mut is None else mut
        sssm_bysite = get_sssm_bysite(mut)
        issm_bysite = get_issm_bysite(mut)

        ## begin with a list of ssm strings
        ## later, will convert list to a single mstring using join
        mstring = []

        for site in self.allsites:
            if site in self.sub_site_ref:
                if mut:
                    ssm = sssm_bysite.get(site,None)
                    ssm_char = ssm.mut if ssm else "."
                else:
                    ssm_char = self.sub_site_ref[site]
                mstring.append(ssm_char)
            if site in self.ins_site_len:
                if mut:
                    ssm = issm_bysite.get(site,None)
                    ssm_str = ssm.mut if ssm else ""
                else:
                    ssm_str = ""
                ssm_str += "~"*(self.ins_site_len[site] - len(ssm_str))
                mstring.append(ssm_str)

        return "".join(mstring)

def _main(args):
    '''main'''
    v.vprint(args)

    label_mstring_list = list(read_label_mstring(args.mstringfile))

    labels = [label for label,_ in label_mstring_list]
    mstrings = [mstring for _,mstring in label_mstring_list]
    mutants = [mutant.Mutation(mstring) for _,mstring in label_mstring_list]

    maxlabellen = max(len(label) for label in labels)
    fmtlabel = "%%%ds " % maxlabellen

    for label,mstring in zip(labels,mstrings):
        v.print((fmtlabel + "%s") % (label,mstring))

    frame = MStringFramework(mutants)

    for label,mut in zip(labels,mutants):
        labeltext = fmtlabel % label if not args.nolabels else ""
        mstring = frame.format_mstring(mut)
        print(labeltext + mstring)


    allsites = list( frame.full_sitelist() )
    for headerstring in intlist.write_numbers_vertically(allsites):
        labeltext = fmtlabel % "" if not args.nolabels else ""
        print(labeltext + headerstring)
    for label,mut in it.chain( [("",None)], zip(labels,mutants) ):
        labeltext = fmtlabel % label if not args.nolabels else ""
        mstring = frame.compact_format_mstring(mut)
        print(labeltext + mstring)

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
