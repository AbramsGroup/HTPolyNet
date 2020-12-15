# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 13:54:55 2020

@author: huang
"""
import os
from countTime import *

class testCase(object):
    def __init__(self):
        self.gro = []
        self.top = []

    def testReadParam(self):
        import readParameters
        a = readParameters.parameters()
        name = '../basic/options.txt'
        a.setName(name)
        a.readParam()
        return a

    @countTime
    def testSearchBonds(self):
        import searchBonds
        import readGro
        import readTop2
        import topInfo
        import groInfo
        
        a = self.testReadParam()
        a1 = []
        a2 = readGro.initGro()
        a3 = readTop2.initTop()
        
        a2.setName('nvt-1')
        # a2.setName('init')
        df_init, sysName, atNum, boxSize = a2.readGRO()
        atomsDf = groInfo.gro()
        atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        atomsDf.initRctInfo(a)
        self.gro = atomsDf
        
        topDf = topInfo.top()
        a3.setName('tmp.top', 'tmp.itp')
        # a3.setName('init.top', 'init.itp')
        a3.genTopSession()
        topDf.setInfo(a3.sumTop)
        topDf.checkCharge()
        self.top = topDf

        print('start searching bonds!!!!')
        aa = searchBonds.searchBonds(a, a1, atomsDf, topDf)
        out = aa.main()
        return out

if __name__ == "__main__":
    a = testCase()
    b1 = a.testSearchBonds()