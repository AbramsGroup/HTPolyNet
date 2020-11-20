# -*- coding: utf-8 -*-
"""
Created on Mon Nov  9 15:30:10 2020

@author: huang
"""
import time
from functools import wraps

def countTime(fn):
    @wraps(fn)
    def measure_time(*args, **kwargs):
        t1 = time.time()
        result = fn(*args, **kwargs)
        t2 = time.time()
        print(f"@timefn: {fn.__name__} took: {t2 - t1: .5f} s")
        return result
    return measure_time

#@countTime
#def readFile(name):
#    with open(name, 'r') as f:
#        f1 = f.readlines()
#        for i in f1:
#            pass
#
#readFile('tmp.gro')