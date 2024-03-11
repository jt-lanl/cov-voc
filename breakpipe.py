'''trap the BrokenPipeError in python programs'''

## Usage
##   from breakpipe import no_broken_pipe
##   ...
##   @no_broken_pipe
##   def main(args):
##       ....
##

import os
import sys
import functools

def wrapper(fcn,*pargs,**kwargs):
    '''
    avoids the bulky BrokenPipeError that arises, eg,
    if output is piped through 'head'
    '''
    try:
        fcn(*pargs,**kwargs)
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        ## following two lines /had/ been advised, and /had/ worked
        ## but /now/ it seems to trigger "I/O operation on closed file"
        ## commenting out those two line, now it's quite again...
        #devnull = os.open(os.devnull, os.O_WRONLY)
        #os.dup2(devnull, sys.stdout.fileno())
        sys.exit(1)  # Python exits with error code 1 on EPIPE
    except KeyboardInterrupt:
        print(" Keyboard Interrupt",file=sys.stderr)
        sys.exit(1) ## Just exit

def no_broken_pipe(fcn):
    '''decorator version of wrapper to avoid BrokenPipeError'''
    @functools.wraps(fcn)
    def wrapped_fcn(*p,**kw):
        wrapper(fcn,*p,**kw)
    return wrapped_fcn
