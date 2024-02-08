'''utiitiy routines for embers-like programs'''
import datetime
import random
import warnings
import itertools as it
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

import verbose as v
import wrapgen
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
        v.vprint_only(5,'try:',f's={s}')
        try:
            yyyy,mm = s.split("-")
            day = random.randint(1,days_in_month(yyyy,mm))
            dt = datetime.date(int(yyyy),int(mm),day)
            v.vprint_only(5,'Random day:',f'string {s} assigned {yyyy}-{mm}-{day}')
            return dt
        except ValueError:
            return None
    warnings.warn(f'Invalid Date {s}')
    return None

def date_from_seqname(sname):
    '''parses date (as a string) from sequence name'''
    ## faster, but less robust, than covid.get_date
    ## Q: is it really faster?
    ## Main difference is that it calls local date_from_iso, which chooses random days if days are not supplied
    ## But that only happens if those bad dates have not already been filtered.  They will have been if any --dates
    ## or --days have been specified, because covid.filter_by_date does that. Since we virtually ALWAYS do confine
    ## those dates, this random day thing is mostly just vestigal code; it's maybe a dubious idea anyway; not as though
    ## there are a lot of them.
    tokens = sname.split('.')
    try:
        date = date_fromiso(tokens[4])
    except IndexError:
        warnings.warn(f"bad date: {sname}")
        date = None
    return date

def get_ord_daterange(datecounter,argsdates):
    '''
    find range of dates, in ordinal numbers, based on:
    datecounter[mutant][date] = count of mutants in date, and
    argsdates is basically args.dates;
    returns a tuple of two two-element tuples: (ord_range,ord_plot_range)
    with each range containing a (ord_min,ord_max) pair
    '''
    all_dateset=set()
    for p in datecounter:
        all_dateset.update(datecounter[p])

    v.vprint("Range of dates:",min(all_dateset),max(all_dateset))
    ordmin = min(all_dateset).toordinal()
    ordmax = max(all_dateset).toordinal()

    ordplotmin = ordmin
    ordplotmax = ordmax
    v.vprint('argsdates:',argsdates)
    if argsdates and argsdates[0] and argsdates[0] != ".":
        ordplotmin = date_fromiso(argsdates[0]).toordinal()
    if argsdates and argsdates[1] and argsdates[1] != ".":
        ordplotmax = date_fromiso(argsdates[1]).toordinal()

    ordmin = min([ordmin,ordplotmin])
    ordmax = max([ordmax,ordplotmax])

    return (ordmin, ordmax), (ordplotmin, ordplotmax)

def filter_seqs_by_padded_dates(seqs,args,firstisref=False):
    '''
    filter input data to keep only data in date range specified...except
    pad the date range so that weekly or cumulative plots are still correct
    side effect: args.dates potentially gets modified
    '''
    start_date,stop_date = covid.date_range_from_args(args)
    args.dates = [start_date,stop_date]
    start_date,stop_date = covid.expand_date_range([start_date,stop_date],args.daily)
    if not (start_date is None and stop_date is None):
        seqs = covid.filter_by_date(seqs,start_date,stop_date,firstisref=firstisref)

    if args.nseq:
        seqs = it.islice(seqs,args.nseq+1)
    if args.verbose:
        seqs = wrapgen.keepcount(seqs,"Sequences filtered:")

    return seqs

def get_running_weekly(cumulative_counts,num_days,daysperweek=7):
    '''convert cumulative counts into running weekly counts'''
    if daysperweek==0:
        return cumulative_counts

    running_weekly=dict()
    for m in cumulative_counts:
        running_weekly[m] = cumulative_counts[m][:]
        for n in range(daysperweek,num_days):
            running_weekly[m][n] = cumulative_counts[m][n] \
                - cumulative_counts[m][n-daysperweek]
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

    running_weekly = get_running_weekly(cumulative_counts,num_days,
                                        daysperweek=daysperweek)
    return running_weekly

def mk_counts_table(date_counter,names):
    '''useful diagnostic; a table of names, counts, onsets, and patterns'''
    ## Uses 'yield' to return lines one at a time
    ## Usage:
    ## for line in mk_counts_table(...):
    ##     print(line)
    ##
    maxnamelen = max(len(names[p]) for p in date_counter)
    fmt = f'%{maxnamelen}s %7s %10s %s'
    yield fmt % ('Name','Count','Onset','Pattern')
    for p in date_counter:
        name = names[p]
        count = sum(date_counter[p].values())
        onset = min(date_counter[p]) if date_counter[p] else ''
        yield fmt % (name,str(count),onset,p)


def get_plot_filename(args,fraction,case_is_none=True,xtra=None):
    '''return name of file in which to save plot'''
    ## eg, something like wk-f-bar-xtra-out.pdf
    if not args.output:
        return ''

    linbar = "line" if args.lineplot else "bar"
    if fraction:
        fc = "f" if case_is_none else "t"
    else:
        fc = "c"
    wk = "wk" if args.daily==7 \
        else "dy" if args.daily == 1 \
             else "cm"

    pplist = [fc,wk,linbar,""]
    if xtra:
        pplist.insert(3,xtra)
    prepend = "-".join(pplist)
    outfile = covid.filename_prepend(prepend,args.output)
    return outfile

def date_friendly(dt):
    '''simple string format for datetime.date object used on x-axis of plots'''
    mmm = dt.strftime("%b")
    dd = dt.strftime("%-d")
    return "%3s %2d" % (mmm, int(dd))

def half_labels(labels,n=2):
    '''replace list of labels ["one","two","three","four","five"]
    with ["one", "", "three", "", "five"]
    '''
    if n<2:
        return labels
    return [label if i%n==0 else ""
            for i,label in enumerate(labels)]

def skip_elements(elist,n=1):
    '''only include every n'th element of a list of elements'''
    if n<1:
        return elist
    return [elt for i,elt in enumerate(elist) if i%n==0]

def embersplot(counts,
               fullnames,
               mcolors,
               ordmin,
               ordplotrange=None,
               num_cases=None,
               title=None,legendtitle=None,legend=0,onsets=False,
               fraction=False,lineplot=False,show=False,daily=7):
    '''
    embersplot is an embers-style plot of counts vs date for multiple variants
    counts: dictionary of arrays; each array is counts vs date;
            keys of dictionary correspond to different variants
    fullnames: dictionary of names associated with dictionary keys (=None if keys are names)
    mcolors: dictionary of colors
    ordmin: first day in array, as ordinal number
    ordplotrange: tuple (ordplotmin, ordplotmax) indicates range of days to plot
    '''

    patterns = list(counts)
    fullnames = fullnames if fullnames else {m:m for m in patterns}
    maxnamelen = max(len(n) for n in fullnames.values())
    namefmt = "%%-%ds" % maxnamelen
    mnames = {m : namefmt % fullnames[m] for m in patterns}

    num_days = len(counts[patterns[0]])
    for m in patterns: ## should make sure len is same for all patterns
        assert num_days == len(counts[m])

    if legend == 0:
        plt.figure(figsize=(6,3))
    elif legend == 1:
        F = 4.7 if lineplot else 5.0
        plt.figure(figsize=(12,1+len(patterns)/F))
    elif legend == 2:
        F = 4.7 if lineplot else 4.5
        plt.figure(figsize=(12,1+len(patterns)/F))
        #plt.figure(figsize=(12,0.5+0.5*len(patterns)/F))
    else:
        raise RuntimeError(f"legend should be 0, 1, or 2: not {legend}")

    if title:
        #title = title + ": %d sequences" % (Nsequences,)
        plt.title(title,fontsize='x-large')

    counts_bottom = dict()
    counts_total = [0] * len(counts[patterns[0]])
    for m in counts:
        counts_bottom[m] = counts_total
        counts_total = [x+y for x,y in zip(counts_total,counts[m])]

    ## at this point counts_total is array (vs date) of total counts over all patterns

    if legendtitle and legend>1 and not lineplot:
        plt.bar(range(num_days),[0]*num_days,width=1,
                label=legendtitle,color="white")

    Y = []
    Ylabels = []
    Ycolors = []

    name_color_sofar = set()
    for m in  patterns: #patterns[::-1]:

        name = mnames[m] ## mnames
        #if legend>1:   ### take this outside
        #    name += " " + mrelpatt[m]
        name = " " + name ## hack! leading underscore doesn't make it to legend??

        if lineplot:
            kwargs = dict(color=mcolors[m],label=name)
        else:
            ## For bar plot, repeated name,color tuples only appear once in legend
            name_color = (name,mcolors[m])
            kwargs = dict(color=mcolors[m])
            if name_color not in name_color_sofar:
                name_color_sofar.add(name_color)
                kwargs['label']=name

        if fraction:
            fm = [a/(b+0.001) for a,b in zip(counts[m],counts_total)]
            bm = [a/(b+0.001) for a,b in zip(counts_bottom[m],counts_total)]
            if num_cases is not None:
                fm = [c*f for c,f in zip(fm,num_cases)]
                bm = [c*f for c,f in zip(bm,num_cases)]

        else:
            fm = counts[m]
            bm = counts_bottom[m]

        if lineplot:
            dy = np.array(range(num_days))
            fm = np.array(fm)
            dy = dy[fm>0]
            fm = fm[fm>0]
            plt.semilogy(dy,fm,lw=2, **kwargs)
        else:
            #plt.bar(range(num_days),fm,bottom=bm,width=1,**kwargs)
            Y.append(fm)
            Ylabels.append(kwargs.get('label',None))
            Ycolors.append(kwargs['color'])

    ## HACK!
    if 0: #num_cases is not None:
        def datemax(values,*xtravalues):
            for xtra in xtravalues:
                values = [v+x for v,x in zip(values,xtra)]
            maxval = max(values)
            ordval = ordmin + values.index(maxval)
            dt = datetime.date.fromordinal(ordval).isoformat()
            return round(maxval),dt
            
            
        valdict = dict()
        datedict = dict()
        for lab,vals in zip(Ylabels,Y):
            lab = lab[1:9]
            valdict[lab] = vals

        valdict['BA.5'] =    [u+v+w for u,v,w in zip(valdict['BA.5/B[E'],
                                                     valdict['BA.4/CS '],
                                                     valdict['BA.4.6/D'])]
        valdict['XBB']  =    [v for v in valdict['XBB     ']]
        valdict['BQ'] =      [u+v for u,v in zip(valdict['BQ.1    '],
                                                 valdict['BQ.1.1/C'])]
        valdict['BA.2.75'] = [u+v for u,v in zip(valdict['BA.2.75.'],
                                                 valdict['BA.2.75/'])]
        print(' PEAKS ',end='')
        for ndx in ['BA.5','BA.2.75','XBB','BQ']:
            print('%9d %s' % datemax(valdict[ndx]),
                  end='')
        print()

        


            
    if not lineplot:
        ## instead of multiple calls to plt.bar, single call to plt.stackplot
        plt.stackplot(range(num_days),Y,
                      labels=Ylabels,colors=Ycolors)

    if fraction and not lineplot and num_cases is None:
        plt.ylim([0,1.05])

    #if fraction:
    #    plt.yticks([0.0001,0.001,0.01,0.1,1])
    #    plt.yticks([0.0001,0.01,1])

    if legend:
        ## reverse the order of the handles/labels in the legend
        ## deleting the empty labels (corresponding to repeated colors)
        handles, labels = plt.gca().get_legend_handles_labels()
        hl_reverse_list = list(zip(handles,labels))[::-1]
        handles = []
        labels = []
        for h,l in hl_reverse_list:
            if l is not None:
                handles.append(h)
                labels.append(l)

        plt.legend(handles,labels,
                   bbox_to_anchor=(1.02, 1),
                   #handlelength=3,
                   #markerfirst=False,
                   frameon=False,
                   handletextpad=0,
                   labelspacing=0.45,
                   loc='upper left', borderaxespad=0.,
                   prop={'family' : 'monospace'})

    if not ordplotrange:
        ordplotmin = ordmin+1
        ordplotmax = ordmin+num_days
    else:
        ordplotmin, ordplotmax = ordplotrange

    plt.xlim(ordplotmin-ordmin-1,ordplotmax-ordmin+1)
    xticks = list(range(ordplotmin-ordmin,ordplotmax-ordmin+1,7)) ## was n+6
    xlabels = [datetime.date.fromordinal(int(ord+ordmin)) for ord in xticks]
    min_date = datetime.date.fromordinal(ordplotmin)
    max_date = datetime.date.fromordinal(ordplotmax)
    xlabels = [date_friendly(dt) for dt in xlabels]
    nhalf = 1 + len(xlabels)//16
    if nhalf<4:
        xlabels = half_labels(xlabels,nhalf)
    else:
        xlabels = skip_elements(xlabels,nhalf)
        xticks  = skip_elements(xticks, nhalf)
    plt.xticks(xticks,xlabels,fontsize='medium', ## was small
               rotation=45,ha='right',position=(0,0.01))

    plt.xlabel(f'{min_date} to {max_date}')

    if onsets:
        ylo,yhi = plt.gca().get_ylim()
        for m in patterns:
            if m not in onsets:
                continue
            x = onsets[m].toordinal() - ordmin
            kwargs=dict(lw=1,color=mcolors[m])
            if not fraction:
                kwargs['zorder']=0
            plt.plot([x,x],[ylo,yhi],**kwargs)

    if lineplot:
        if legend==0:
            plt.subplots_adjust(bottom=0.15,right=0.95) ## hardcoded hack!
        else:
            #plt.subplots_adjust(bottom=0.15,right=0.8) ## still hardcoded!
            #plt.subplots_adjust(bottom=0.15,right=0.7) ## still hardcoded!
            pass

        def ytickformat(x,pos):
            if x>1:
                return "" if (fraction and num_cases is None) else "%g" % x
            else:
                return "%.1g" % x

        plt.gca().yaxis.set_major_formatter(mpl.ticker.FuncFormatter(ytickformat))

    plt.yticks(fontsize=12)

    if daily==7:
        ylabel_prefix="Weekly"
    elif daily==1:
        ylabel_prefix="Daily"
    elif daily==0:
        ylabel_prefix="Cumulative"
    else:
        ylabel_prefix=f'{daily}-day'



    plt.ylabel("Fraction" if fraction else ylabel_prefix+" Sequence Counts")
    if fraction and num_cases is not None:
        plt.ylabel(ylabel_prefix+" Cases")
    plt.tight_layout()
    if show:
        plt.show()


def make_emberstyle_plots(args,extra_id,cum_counts,
                          names,colors,ordmin,ordplotrange,
                          num_cases=None,
                          title=None,nmatches=0,ncases=0,
                          daily=7,
                          onsets=None):
    '''
    plot fraction, sequence counts, and cases
    against time for variants of interest
    '''

    ## restrict data to within the plot range
    ## then, auto scaling on the y-axis will be based on available data
    ordplotmin,ordplotmax = ordplotrange
    assert ordmin <= ordplotmin
    for m in cum_counts:
        cum_counts[m] = cum_counts[m][ordplotmin-ordmin: ordplotmax-ordmin+1]
    ordmin = ordplotmin

    ## Now make three plots, one w/ seq-counts, one w/ fractions,...
    ## and one w/ case counts (if num_counts is available)
    fractions = (False,True)
    cases = (None,None)
    totals = (f'{nmatches} sequences',
              f'{nmatches} sequences')
    if num_cases is not None:
        fractions = (False,True,True)
        cases = (None,None,num_cases)
        totals = (f'{nmatches} sequences',
                  f'{nmatches} sequences',
                  #f'{nmatches} sequences',
                  f'{ncases} cases',
        )

    for fraction,case,total in zip(fractions,cases,totals):
        the_title = title + f': {total}'
        embersplot(cum_counts,names,colors,ordmin,
                   ordplotrange = ordplotrange,
                   num_cases=case,
                   legend=args.legend,lineplot=args.lineplot,
                   title = the_title, onsets=onsets,
                   fraction=fraction, daily=daily)

        if args.output:
            outfile = get_plot_filename(args,fraction,
                                        case is None,xtra=extra_id)
            plt.savefig(outfile)

    if not args.output:
        plt.show()
