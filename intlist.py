def str_intgen(x,end=None):
    """
    Convert a string such as 2,5-8,10
    Into a generator of integers (2,5,6,7,8,10)
    Also works with :-notation, eg 2,5:9:2,10 becomes (2,5,7,10)
    Note: no sorting or removal of duplicates; sorted(set(...)) will do that
    """
    if x:
        for part in x.split(','):
            if ':' in part:
                abc = part.split(':')
                intabc = [int(_) for _ in abc]
                for n in range(*intabc):
                    yield n
            elif '-' in part:
                a, b = part.split('-')
                if not b and end is not None:
                    b=end
                a, b = int(a), int(b)
                for n in range(a, b + 1):
                    yield n
            else:
                n = int(part)
                yield n

def string_to_intlist(x,end=None):
    return list(str_intgen(x,end=end))

    
def string_to_intlist_old(x,end=None):
    """
    Convert a string such as 2,5-8,10
    Into a list of integers [2,5,6,7,8,10]
    Also works with :-notation, eg 2,5:9:2,10 becomes [2,5,7,10]
    Note: no sorting or removal of duplicates; sorted(set(...)) will do that
    """
    result = []
    for part in x.split(','):
        if ':' in part:
            abc = part.split(':')
            intabc = [int(_) for _ in abc]
            result.extend(range(*intabc))
            #print "part=",part,"intabc=",intabc,"result=",result #debug
        elif '-' in part:
            a, b = part.split('-')
            if not b and end is not None:
                b=end
            a, b = int(a), int(b)
            result.extend(range(a, b + 1))
        elif "" == part:
            pass
        else:
            a = int(part)
            result.append(a)
    return result

def format_intlist(nlist,width=79,intro=""):
    '''write out a list of integers neatly over several lines,
    separated by commas, and lined up nicely in columns. eg:
    nlist = range(25), width=50, intro="Sites:" produces
    Sites:  0,  1,  2,  3,  4,  5,  6,  7,  8,  9,
           10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
           20, 21, 22, 23, 24
    '''
    maxwid = max(len(str(n)) for n in nlist)
    fmt = f"%{maxwid}d"
    nperline = (width - len(intro)) // (maxwid+2)
    nlines = 1 + len(nlist) // nperline
    lines=[]
    for nl in range(nlines):
        introstr = intro if nl==0 else " "*len(intro)
        nbot = nl*nperline
        ntop = min([nl*nperline + nperline,len(nlist)])
        nslice = slice(nbot, ntop)
        if ntop>nbot:
            lines.append(introstr + ", ".join(fmt % n for n in nlist[nslice]))
    return lines

def intlist_to_string(nlist,sort=False):
    '''inverse of string_to_intlist, 
    eg, [1,2,3,6,7,9,10,11] --> "1-3,6,7,9-11"
    '''
    if len(nlist)==0:
        return ""
    if sort:
        ## sort and also remove duplicates
        nlist = sorted(set(nlist))
    slist=[]
    nlo = nhi = nlist[0]
    for n in list(nlist[1:]) + [None]:
        if n == nhi+1:
            nhi = n
        else:
            if nlo==nhi:
                slist.append( str(nlo) )
            elif nhi==nlo+1:
                slist.append( str(nlo) + "," + str(nhi) )
            elif nhi>nlo+1:
                slist.append( str(nlo) + "-" + str(nhi) )
            nlo = nhi = n
    
    return ",".join(slist)
        
    

def write_numbers_vertically(nlist,plusone=0,leadingzero=' '):
    '''
    write numbers vertically so one number per column;
    eg, input: nlist = [80, 96, 175, 515]
    output three lines:
    ..15
    8971
    0655
    '''

    if nlist is None or len(nlist)==0:
        return []
    
    lines = []
    rlist = []
    nlist = [n+plusone for n in nlist]
    while max(nlist)>0:
        rlist.append( [str(n%10) for n in nlist] )
        nlist = [n//10 for n in nlist]
    rlist = rlist[::-1]

    ## Replace leading 0's with space
    for j in range(len(rlist[0])):
        for k in range(len(rlist)):
            if rlist[k][j] == '0':
                rlist[k][j] = leadingzero
            else:
                break
            
    for r in rlist:
        lines.append("".join(r))
    return lines


if __name__ == "__main__":

    s="2,5-8,10"
    print(s,string_to_intlist(s))
    s="2,5-8,10-"
    print(s,string_to_intlist(s,end=20))

    print("Empty:",string_to_intlist_old(""))
    print("None:",string_to_intlist(None))
    print("Empty:",string_to_intlist(""))
    print("Zero:",string_to_intlist("0"))

    nlist = [80,90,145,2,3567,3456,124]
    lines = write_numbers_vertically(nlist)
    print(nlist,":")
    print("\n".join(lines))
    
          
