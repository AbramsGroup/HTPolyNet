"""

.. module:: countTime
   :synopsis: defines a wrapper implementing method timing
   
.. moduleauthor: Ming Huang, <mh3429@dragons.drexel.edu>

"""
import os
import time
from functools import wraps

def countTime(fn):
    @wraps(fn)
    def measure_time(*args, **kwargs):
        t1 = time.time()
        result = fn(*args, **kwargs)
        t2 = time.time()
        if os.path.isfile('time.txt'):
            with open('time.txt', 'a') as f:
                f.write(f"@timefn: {fn.__name__} took: {t2 - t1: .5f} s\n")
        else:
            with open('time.txt', 'w') as f:
                f.write(f"@timefn: {fn.__name__} took: {t2 - t1: .5f} s\n")
        return result
    return measure_time
