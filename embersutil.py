'''utiitiy routines for embers-like programs'''
import datetime
import random
import warnings
import matplotlib.pyplot as plt

import verbose as v
from embersplot import embersplot
import covid

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

def days_in_month(yyyy,mm):
    '''return number of days in month yyyy-mm'''
    # https://stackoverflow.com/questions/42950/how-to-get-the-last-day-of-the-month
    next_month = datetime.date(int(yyyy),int(mm),28) + datetime.timedelta(days=4)
    return (next_month - datetime.timedelta(days=next_month.day)).day

def date_fromiso(s):
    '''return datetime.date object from iso string yyyy-mm-dd'''
    ## over-rides covid.date_fromiso by giving yyyy-mm a random day
    if isinstance(s,datetime.date):
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
        except ValueError:
            return None
    warnings.warn(f'Invalid Date {s}')
    return None

def date_from_seqname(sname):
    '''parses date (as a string) from sequence name'''
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


def get_running_weekly(cumulative_counts,num_days,daysperweek=7):
    '''convert cumulative counts into running weekly counts'''
    running_weekly=dict()
    for m in cumulative_counts:
        running_weekly[m] = cumulative_counts[m][:]
        for n in range(daysperweek,num_days):
            running_weekly[m][n] = cumulative_counts[m][n] - cumulative_counts[m][n-daysperweek]
    return running_weekly

def get_cumulative_counts(date_counter,ord_range,daysperweek=7):
    '''convert date_counter into list of cumulative counts'''
    ordmin, ordmax = ord_range
    num_days = ordmax+1-ordmin
    v.vprint("Days:",num_days,ordmin,ordmax)

    cumulative_counts={m: list() for m in date_counter}
    for ord_date in range(ordmin,ordmax+1):
        day = datetime.date.fromordinal(ord_date)
        for m in date_counter:
            cumulative_counts[m].append( sum(date_counter[m][dt]
                                  for dt in date_counter[m] if dt <= day ) )

    running_weekly = get_running_weekly(cumulative_counts,num_days,daysperweek=daysperweek)
    return running_weekly

def get_plot_filename(args,fraction,xtra=None):
    '''return name of file in which to save plot'''
    ## eg, something like wk-f-bar-xtra-out.pdf
    if not args.output:
        return ''

    linbar = "line" if args.lineplot else "bar"
    fc = "f" if fraction else "c"
    wk = "wk" if args.daily==7 \
        else "dy" if args.daily == 1 \
             else "cm"

    pplist = [fc,wk,linbar,""]
    if xtra:
        pplist.insert(3,xtra)
    prepend = "-".join(pplist)
    outfile = covid.filename_prepend(prepend,args.output)
    return outfile

def make_emberstyle_plots(args,id,cum_counts,names,colors,ordmin,ordplotrange,
                          title=None, onsets=None):
    '''plot both fraction and counts against time for variants of interest'''
    for fraction in (False,True):
        embersplot(cum_counts,names,colors,ordmin,
                   ordplotrange = ordplotrange,
                   legend=args.legend,lineplot=args.lineplot,
                   title = title, onsets=onsets,
                   fraction=fraction)

        if args.output:
            outfile = get_plot_filename(args,fraction,xtra=id)
            plt.savefig(outfile)

    if not args.output:
        plt.show()
