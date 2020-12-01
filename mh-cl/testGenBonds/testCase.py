# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 13:54:55 2020

@author: huang
"""
import os
from countTime import *
import pandas as pd

class testCase(object):
    def __init__(self):
        self.gro = []
        self.top = []

    @countTime
    def testReadParam(self):
        import readParameters
        a = readParameters.parameters()
        name = 'basic/options.txt'
        a.setName(name)
        a.readParam()
        return a

    def testGenBonds(self, pairs):
        import genBonds
        rctMols = ['2', '14']
        def getChargeMaps():
            maps = {}
            with open('charges.txt', 'r') as f:
                idx = 0
                for i in f.readlines():
                    key, value = i.split(':')
                    if key in maps.keys():
                        pass
                    else:
                        maps[key] = value
            return maps
        
        b = getChargeMaps()

        a = genBonds.genBonds(self.gro, self.top, pairs, b, rctMols)
        a.main()
        return a
        
if __name__ == "__main__":
    a = testCase()
    names = ['acro', 'amon']; tmp = [['123', '1177']]
    df_pairs = pd.DataFrame(tmp, columns=names)
    b2 = a.testGenBonds(df_pairs)
