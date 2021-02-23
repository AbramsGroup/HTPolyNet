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

    @countTime
    def testReadParam(self):
        import readParameters
        a = readParameters.parameters()
        name = 'basic/options.txt'
        a.setName(name)
        a.readParam()
        return a
    
    @countTime
    def testreadGro(self, name='systems/unrctSystem/VEA'):
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
    
    @countTime
    def testMain(self): # cannot test locally, due to the parmed module cannot install locally
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
        
        # a2.setName('init')
        a2.setName('init')
        df_init, sysName, atNum, boxSize = a2.readGRO()
        atomsDf = groInfo.gro()
        boxSize = [5.00000, 5.00000, 5.00000] # TODO: remove latter
        atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        atomsDf.initRctInfo(a)
        self.gro = atomsDf
        
        topDf = topInfo.top()
        a3.setName('init.top', 'init.itp')
        # a3.setName('init-1.top', 'init-1.itp')
#        a3.setName('systems/VEA.top', 'systems/VEA.itp')
        a3.genTopSession()

        topDf.setInfo(a3.sumTop)
        topDf.checkCharge()
        self.top = topDf
        atomsDf.df_atoms.to_csv('atomsDf.csv')
        aa = searchBonds.searchBonds(a, a1, atomsDf, topDf, 0, 20)
        pairs, rctMols, cutoff = aa.main()
        return pairs, rctMols, topDf, aa
        # return aa

    @countTime
    def testGenBonds(self, pairs, rctMols):
        import genBonds
        # rctMols = ['2', '14']
        def getChargeMaps():
            maps = {}
            with open('basic/charges.txt', 'r') as f:
                idx = 0
                for i in f.readlines():
                    key, value = i.split(':')
                    if key in maps.keys():
                        pass
                    else:
                        maps[key] = value
            return maps
        
        b = getChargeMaps()

        a = genBonds.genBonds(self.gro, self.top, pairs, b, rctMols, cat='map')
        a.main()
        return a
    
    def testmolRctInfo(self):
        import molRctInfo
        
        path = os.getcwd() + '/' + 'systems/rctSystem'
        a = molRctInfo.molRctInfo(os.getcwd(), path)
        a.getNames(['VEA'])
        
if __name__ == "__main__":
    a = testCase()
    b1 = a.testreadGro()
    # b2 = a.testreadGro('systems/unrctSystem/STY')
    b = a.testReadParam()
    
    pairs, rctMols, topDf, aa = a.testSearchBonds()
    # df_pairs = pairs
    # print('{} bonds will be formed!'.format(len(df_pairs)))
    #
    # b2 = a.testGenBonds(df_pairs, rctMols)
    # #
    # a = b2.top.outDf('tmp11', k=0.9, stepRelax=True)