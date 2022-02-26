'''utiitiy routines for embers-like programs'''
import re
import argparse
import datetime

import lineagetable
from embersplot import embersplot
import verbose as v

import warnings

OTHER='other'

def embers_args(ap):
    '''
    call this in the getargs() function
    so these options will be added in
    '''    
    paa = ap.add_argument
    
    paa("--daily",type=int,default=7,
        help="daily=1 for daily, daily=7 (default) for weekly, daily=0 for cumulative")
    paa("--weekly",action="store_const",const=7,dest='daily',
        help="Make weekly average plots instead of cumulative (set daily=7)")
    paa("--cumulative",action="store_const",const=0,dest='daily',
        help="Make cumulative plots, instead of daily or weekly averaged")
    ap.set_defaults(daily=7)
    
    g = ap.add_argument_group('Graphing options')
    gaa = g.add_argument
    gaa("--title",
        help="Use this TITLE in plots")
    gaa("--lineplot",action="store_true",
        help="Make log-linear line plot instead of linear stacked bar plot")
    gaa("--onsets",action="store_true",
        help="plot onset dates for each mutant")
    gaa("--legend",type=int,default=1,choices=(0,1,2),
        help="0: no legend, 1: legend, 2: big legend (with seq patterns)")
    gaa("--output","-o",help="write plot to file")

    return

def days_in_month(yyyy,mm):
    # https://stackoverflow.com/questions/42950/how-to-get-the-last-day-of-the-month
    next_month = datetime.date(int(yyyy),int(mm),28) + datetime.timedelta(days=4)
    return (next_month - datetime.timedelta(days=next_month.day)).day

def date_fromiso(s):
    ## over-rides covid.date_fromiso by giving yyyy-mm a random day
    if type(s) == datetime.date:
        return s
    try:
        yyyy,mm,dd = s.split("-")
        dt = datetime.date(int(yyyy),int(mm),int(dd))
        return dt
    except ValueError: 
        if s == ".":
            return None
        try:
            yyyy,mm = s.split("-")
            day = random.randint(1,days_in_month(yyyy,mm))
            dt = datetime.date(int(yyyy),int(mm),day)
            return dt
        except ValueError as e:
            return None
    warnings.warn(f'Invalid Date {s}')
    return None


def date_from_seqname(sname):
    ## faster, but less robust, than covid.date_from_seqname
    tokens = sname.split('.')
    try:
        date = date_fromiso(tokens[4])
    except IndexError:
        warnings.warn(f"bad date: {sname}")
        date = None
    return date


def get_daterange(datecounter,argsdates):
    ''' find range of dates, in ordinal numbers, based on:
    datecounter[mutant][date] = count of mutants in date, and
    argsdates is basically args.dates '''
    
    all_dateset=set()
    for p in datecounter:
        all_dateset.update(datecounter[p])

    v.vprint("Range of dates:",min(all_dateset),max(all_dateset))
    ordmin = min(all_dateset).toordinal()
    ordmax = max(all_dateset).toordinal()

    ordplotmin = ordmin
    ordplotmax = ordmax
    if argsdates and argsdates[0] and argsdates[0] != ".":
        ordplotmin = date_fromiso(argsdates[0]).toordinal()
    if argsdates and argsdates[1] and argsdates[1] != ".":
        ordplotmax = date_fromiso(argsdates[1]).toordinal()

    ordmin = min([ordmin,ordplotmin])
    ordmax = max([ordmax,ordplotmax])
    return ordmin, ordmax, ordplotmin, ordplotmax


def epoch_average(DG_cum,Ndays,daysperweek=7):
    DG_weekly=dict()
    for m in DG_cum:
        DG_weekly[m] = DG_cum[m][:]
        for n in range(daysperweek,Ndays):
            DG_weekly[m][n] = DG_cum[m][n] - DG_cum[m][n-daysperweek]
        DG_cum[m] = DG_weekly[m]
    return DG_cum

def get_cumulative_counts(DG_datecounter,ord_tuple,daysperweek=7):
    
    ordmin, ordmax, ordplotmin, ordplotmax = ord_tuple
    Ndays = ordmax+1-ordmin
    v.vprint("Days:",Ndays,ordmin,ordmax)
   
    DG_cum={m: list() for m in DG_datecounter}
    for ord in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord)
        for m in DG_datecounter:
            DG_cum[m].append( sum(DG_datecounter[m][dt]
                                  for dt in DG_datecounter[m] if dt <= day ) )

    DG_cum = epoch_average(DG_cum,Ndays,daysperweek=daysperweek)
    return DG_cum

def rd_lineage_table(filename):
    ## this code is currently broken
    lineage_table=[]
    with open(filename) as ftable:
        for line in ftable:
            line = re.sub(r'#.*','',line)
            line = line.strip()
            if not line:
                continue
            try:
                color,name,pattern = line.split()
                lineage_table.append( (color,name,pattern) )
            except:
                v.vprint(f'Bad line in lineage file: {line}')
    return lineage_table

def get_lineage_table(lineagetable_file=None):
    
    table = lineagetable.LineageTable
    if lineagetable_file:
        ## if file exists, over-ride default
        table = rd_lineage_table(lineagetable_file)
    table.insert(0,('Gainsboro', 'other', OTHER))
    return table
