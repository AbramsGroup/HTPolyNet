# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 13:54:55 2020

@author: huang
"""
import os
from countTime import *
import pandas as pd
import readGro
import readTop
import groInfo
import topInfo

class testCase(object):
    def __init__(self):
        self.gro = []
        self.top = []

    @countTime
    def testReadParam(self):
        import readParameters
        a = readParameters.parameters()
        name = '../basic/options.txt'
        a.setName(name)
        a.readParam()
        return a

    def testGenBonds(self, pairs, cat='pd'):
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

        a1 = genBonds.genBonds(self.gro, self.top, pairs, b, rctMols, cat=cat)
        a1.main()
        return a1
        
if __name__ == "__main__":
    a = testCase()
    names = ['acro', 'amon']; tmp = [['121', '1175']]
    df_pairs = pd.DataFrame(tmp, columns=names)

    a1 = a.testReadParam()
    atomsDf = groInfo.gro()
    topDf = topInfo.top()

    a2 = readGro.initGro()
    a3 = readTop.initTop()

    a2.setName('tmp')
    df_init, sysName, atNum, boxSize = a2.readGRO()
    atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
    atomsDf.initRctInfo(a1)

    a3.setName('tmp.top', 'tmp.itp')
    a3.genTopSession()
    topDf.setInfo(a3.sumTop)
    topDf.checkCharge()

    a.gro = atomsDf
    a.top = topDf


    # b2 = a.testGenBonds(df_pairs, cat='pd')

    b2 = a.testGenBonds(df_pairs, cat='map')
