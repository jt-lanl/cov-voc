'''
Routine for wrapping a generator, and counting the items as they go by;
finally printing the number of items either when the generator is exhausted,
or some other error condition arises, or the program exits.
'''

import sys
from collections.abc import Iterable, Iterator

def _keepcount_iterator(seqs,msg_prefix="",file=sys.stderr):
    '''
    seqs is an iterator (eg, a generator)
    return value is an iterator that wraps seqs 
    and keeps count of how many items are consumed
    '''
    assert isinstance(seqs,Iterator)
    n=0
    try:
        while True:
            n += 1
            yield next(seqs)
    except StopIteration:
        n -= 1 ## undo the n+=1 before for next(seqs) above
    finally:
        print(msg_prefix,n,file=file)


def keepcount(seqs,msg_prefix="",file=sys.stderr):
    '''
    seqs is an iterable; could be a list or generator;
    if seqs is a list, print length of list, and return seqs (still as a list)
    '''

    if isinstance(seqs,Iterator):
        return _keepcount_iterator(seqs,msg_prefix=msg_prefix,file=file)
    elif isinstance(seqs,Iterable):
        print(msg_prefix,len(seqs),file=file)
        return seqs
    else:
        raise RuntimeError(f"[msg_prefix] input is not an iterable!")


        
if __name__ == '__main__':

    ### Here's some code to demonstrate use of wrapgen.keepcount

    gen = (n for n in range(10))
    gen = keepcount(gen,"gen items:")

    print([x for x in gen])

    gen = (n for n in range(10))
    gen = keepcount(gen,"xgen items:")

    #finish early
    for m,x in enumerate(gen):
        print(m,x)
        if m>5:
            break
    
    lst = [n for n in range(10)]
    lxst = keepcount(lst,"lst items:")
    print("lxst:",type(lst),type(lxst))
    #print("lxst:",lxst[:2])
    print("lst:",[x for x in lst])
    
    for m,x in enumerate(lxst):
        print(m,x)
        if m>5:
            break

    j = keepcount(7,"hello")  ## This SHOULD raise an error
        
    #lxst = list(lxst)
    
