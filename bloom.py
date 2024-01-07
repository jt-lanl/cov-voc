'''Read summary tsv file, based on experiments by Jesse Bloom'''

import re
import argparse
import numpy as np
import pandas as pd

import verbose as v

def _getargs():
    '''parse options from command line'''
    argparser = argparse.ArgumentParser(description=__doc__)
    paa = argparser.add_argument
    paa("--input","-i",
        help="input tsv file")
    paa("--verbose","-v",action="count",default=0,
        help="verbose")
    args = argparser.parse_args()
    return args

class BloomSummary:
    def __init__(self, filepath):
        self.df = pd.read_table(filepath)

    def match_column_name(self,input_name):
        '''find the true column name that column_name aspires to'''
        the_true_names = [] ## list of true names that match
        if input_name in self.df.columns:
            return input_name
        for true_name in self.df.columns:
            if input_name.strip().lower() in true_name.lower():
                the_true_names.append(true_name)
        if len(the_true_names) == 1:
            return the_true_names[0]
        if len(the_true_names) == 0:
            raise RuntimeError(f'No match found for name={input_name}')
        if len(the_true_names) > 1:
            raise RuntimeError(f'Too many matches for name={input_name}:\n'
                               f'{the_true_names}')        

    def get_ssm(self,row):
        return (row.wildtype + str(row.site) + row.mutant)

    def get_value(self,row,column):
        true_column = self.match_column_name(column)
        return row[true_column]
        
    
    def info_by_mutation(self,ssm):
        '''ssm is single-site mutation string; eg "D614G"'''
        match = re.match(r'([A-Z])(\d+)([A-Z-])',ssm.strip())
        if not match:
            return None
        wildtype,site,mutant = match.group(1,2,3)
        rows = self.df[ self.df["site"] == int(site) ]
        rows = rows[ rows["mutant"] == mutant ]
        return rows

    def toprows(self,column,n=100,ascending=False):
        ''' get the top 100 items, based on column '''
        true_column = self.match_column_name(column)
        sorted_df = self.df.sort_values(by=true_column,
                                        ascending=ascending,
                                        na_position='last')
        sorted_df = sorted_df[:n]
        sorted_df = sorted_df.sort_values(by="site",ascending=True)
        
        return sorted_df

    def plotable(self,df,column):
        ''' return data in a dataframe as a plot-able pair of lists'''
        true_column = self.match_column_name(column)
        xvals = []
        yvals = []
        for _,row in df.iterrows():
            xvals.append(self.get_ssm(row))
            yvals.append(row[true_column])
        return xvals,yvals

def _main(args):
    '''main'''
    v.vprint(args)
    bloomer = BloomSummary(args.input)
    v.print('F4S:\n',
            bloomer.info_by_mutation('F4S'))

    ## get the biggest "human sera" type:
    top = bloomer.toprows("escape",n=10)
    v.print(top)

    xvals,yvals = bloomer.plotable(top,"escape")
    for x,y in zip(xvals,yvals):
        print(x,y)
    
    

if __name__ == "__main__":

    _args = _getargs()
    v.verbosity(_args.verbose)
    _main(_args)
