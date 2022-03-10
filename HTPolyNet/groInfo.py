# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 11:01:01 2020

@author: huangming
"""

import pandas as pd
from HTPolyNet.countTime import *
from copy import deepcopy

class gro(object):
    def __init__(self):
        self.df_atoms = []
        self.sysName = ''
        self.atNum = ''
        self.boxSize = ''

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        result.__dict__.update(deepcopy(self.__dict__))
        return result

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

    def getProp(self, molName, atName, infoMap):
        info = {}
        for i in infoMap[molName]:
            if i[0] == atName:
                info[atName] = [i[1], i[2], i[3]] # rctNum & at prop
        return info

    @countTime
    def initRctInfo(self, basicParameters):
        monR_list = basicParameters.monR_list
        croR_list = basicParameters.croR_list
        infoMap = {**monR_list, **croR_list}

        self.df_atoms['rct'] = 'None'
        self.df_atoms['rctNum'] = '0'
        self.df_atoms['prop'] = 'N'
        self.df_atoms['rctGroup'] = '0'

        self.df_atoms['groupCon'] = ''
        self.df_atoms['chain'] = ''

        for i in monR_list.keys():
            atNames = []
            for n in monR_list[i]:
                info = self.getProp(i, n[0], infoMap)
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rct'] = 'True'
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rctNum'] = int(info[n[0]][0])
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'prop'] = info[n[0]][1]
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rctGroup'] = info[n[0]][2]

        for i in croR_list.keys():
            atNames = []
            for n in croR_list[i]:
                info = self.getProp(i, n[0], infoMap)
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rct'] = 'True'
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rctNum'] = int(info[n[0]][0])
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'prop'] = info[n[0]][1]
                self.df_atoms.loc[(self.df_atoms.molName == i) & (self.df_atoms.atomName == n[0]), 'rctGroup'] = info[n[0]][2]

    @countTime
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