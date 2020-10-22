# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 13:54:55 2020

@author: huang
"""
import os

class testCase(object):
    def __init__(self):
        self.gro = []
        self.top = []
    
    def testReadParam(self):
        import readParameters
        a = readParameters.parameters()
        name = 'basic/options.txt'
        a.setName(name)
        a.readParam()
        return a
    
    def testreadGro(self, name='systems/VEA'):
        import readGro
        import groInfo
        pp = self.testReadParam()
        a = readGro.initGro()
        a.setName(name)
        df_init, sysName, atNum, boxSize = a.readGRO()
        a1 = groInfo.gro()
        a1.setGroInfo(df_init, sysName, atNum, boxSize)
        a1.initRctInfo(pp)
        return a1
    
    def testMain(self):
        import os
        import main
        
        path = os.getcwd()
        basicFolder = '{}/{}'.format(path, 'basic')
        systemsFolder = '{}/{}'.format(path, 'systems')
        mdpFolder = '{}/{}'.format(path, 'mdp')
        a = main.main(basicFolder, systemsFolder, mdpFolder)
        a.setParam('options.txt')
        a.initSys('1')
        return a.basicParameter
    
    def testSearchBonds(self):
        import searchBonds
        import readGro
        import readTop
        import topInfo
        import groInfo
        
        a = self.testReadParam()
        a1 = []
        a2 = readGro.initGro()
        a3 = readTop.initTop()
        
        a2.setName('nvt-1')
        df_init, sysName, atNum, boxSize = a2.readGRO()
        atomsDf = groInfo.gro()
        atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        atomsDf.initRctInfo(a)
        self.gro = atomsDf
        
        topDf = topInfo.top()
        a3.setName('init.top', 'init.itp')
#        a3.setName('systems/VEA.top', 'systems/VEA.itp')
        a3.genTopSession()
        topDf.setInfo(a3.sumTop)
        topDf.checkCharge()
        self.top = topDf
        aa = searchBonds.searchBonds(a, a1, atomsDf, topDf)
        out = aa.main()
        return out
    
    def testGenBonds(self, pairs):
        import genBonds
        a = genBonds.genBonds(self.gro, self.top, pairs)
        a.main()
        return a
    
    def testmolRctInfo(self):
        import molRctInfo
        
        path = os.getcwd() + '/' + 'systems/rctSystem'
        a = molRctInfo.molRctInfo(os.getcwd(), path)
        a.getNames(['VEA'])
        
if __name__ == "__main__":
    a = testCase()
#    b1 = a.testreadGro()
    
#    b2 = a.testreadGro('systems/STY')
#    a1 = b2.df_atoms
    
#    b = a.testReadParam()        
    
#    b = a.testMain()
    
    b1 = a.testSearchBonds()
#    b2 = a.testGenBonds(b1)
    
#    a.testmolRctInfo()