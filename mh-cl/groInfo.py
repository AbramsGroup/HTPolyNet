# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:01:01 2020

@author: huangming
"""

import pandas as pd
from countTime import *

class gro(object):
    def __init__(self):
        self.df_atoms = []
        self.sysName = ''
        self.atNum = ''
        self.boxSize = ''
    
    def setGroInfo(self, atoms_df, sysName, atNum, boxSize):
        self.df_atoms = atoms_df
        self.sysName = sysName
        self.atNum = atNum
        self.boxSize = boxSize
    
    def mergeRow(self, x):
        molNum1 = x['molNum']; molName1 = x['molName']; atomName1 = x['atomName']; globalIdx1 = x['globalIdx']
        posX1 = x['posX']; posY1 = x['posY']; posZ1 = x['posZ']
        if int(globalIdx1) > 99999:
                str1 = '{:>5}{:>5}{:>5}{:>5}{:>8}{:>8}{:>8}'.format(molNum1, molName1, atomName1, 
                        globalIdx1 - 99999, posX1, posY1, posZ1)
        else:
            str1 = '{:>5}{:>5}{:>5}{:>5}{:>8}{:>8}{:>8}'.format(molNum1, molName1, atomName1, 
                    globalIdx1 , posX1, posY1, posZ1)        
        return str1
    
    def updateCoord(self, inGro): # inGro is a gro object
        self.boxSize = inGro.boxSize
        newX = inGro.df_atoms['posX']
        newY = inGro.df_atoms['posY']
        newZ = inGro.df_atoms['posZ']
        self.df_atoms['posX'] = newX
        self.df_atoms['posY'] = newY
        self.df_atoms['posZ'] = newZ
    
    @countTime
    def initRctInfo(self, basicParameters):
        monR_list = basicParameters.monR_list
        croR_list = basicParameters.croR_list

        self.df_atoms['rct'] = 'False'
        self.df_atoms['rctNum'] = '0'

        for i in monR_list.keys():
            atNames = []
            for n in monR_list[i]:
                atNames.append(n[0])

            self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName.isin(atNames)), 'rct'] = 'True'
            self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName.isin(atNames)), 'rctNum'] = 1 # TODO: need to change when rctNum is different
        
        for i in croR_list.keys():
            atNames = []
            for n in croR_list[i]:
                atNames.append(n[0])
            self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName.isin(atNames)), 'rct'] = 'True'
            self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName.isin(atNames)), 'rctNum'] = 1 # TODO: need to change when rctNum is different

    def outDf(self, outName):
        df = self.df_atoms
        sysName = self.sysName
        atNum = self.atNum
        boxSize = self.boxSize
        outName = '{}.gro'.format(outName)
        df1 = pd.DataFrame(); df2 = pd.DataFrame()
        df1 = df1.append([sysName]); df1 = df1.append([atNum])
        df2.loc[:, 0] = df.apply(lambda x: self.mergeRow(x), axis=1)
        df2 = df2.append([boxSize])
        df1.to_csv(outName, mode = 'w', index=False, header=None)
        df2.to_csv(outName, mode = 'a', index=False, header=None)
        # return df2