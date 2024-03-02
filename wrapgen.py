'''
Routine for wrapping a generator, and counting the items as they go by;
finally printing the number of items either when the generator is exhausted,
or some other error condition arises, or the program exits.
'''

import sys
from collections.abc import Iterable, Iterator

def _keepcount_iterator_action(seqs,action):
    '''
    seqs is an iterator (eg, a generator)
    action is a function with argument n
    and that function will be called when the generator is done
    return value is an iterator that keeps count of n, the number of items consumed
   '''
    assert isinstance(seqs,Iterator)
    ncount=0
    for s in seqs:
        yield s
        ncount += 1
    action(ncount)

def _keepcount_iterator(seqs,msg_prefix=None,file=sys.stderr):
    '''
    seqs is an iterator (eg, a generator)
    return value is an iterator that keeps count of how many items are consumed
    when the iterator is done, a message is printed to file
    and
    '''
    def action(ncnt):
        if msg_prefix:
            print(msg_prefix,ncnt,file=file)
        else:
            print(ncnt,file=file)

    yield from _keepcount_iterator_action(seqs,action)

def keepcount(seqs,msg_prefix="",file=sys.stderr):
    '''
    seqs is an iterable; could be a list or generator;
    if seqs is a list, print length of list, and return seqs (still as a list)
    '''
    if isinstance(seqs,Iterator):
        return _keepcount_iterator(seqs,msg_prefix=msg_prefix,file=file)
    if isinstance(seqs,Iterable):
        print(msg_prefix,len(seqs),file=file)
        return seqs
    ## if we've made it this far, something is wrong!
    raise RuntimeError(f"[{msg_prefix}] input is not an iterable!")

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
