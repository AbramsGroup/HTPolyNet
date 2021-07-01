# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:51:30 2020

@author: huang
"""
import pandas as pd
import numpy as np
import re

class mol2Info(object):
    def __init__(self):
        self.initIdx = 0
        self.basicInfo = ''
        self.atoms = ''
        self.bonds = ''
        self.com = ''
        self.mainResname = ''
        self.seq = []
    def setInfo(self, inObj):
        self.basicInfo = inObj.basicInfo
        self.atoms = inObj.atoms.copy(deep=True)
        self.bonds = inObj.bonds.copy(deep=True)
        self.com = inObj.com
    
    def updateResName(self, name):
        self.atoms.resname = name
        
    def updateAtIdx(self, df): # Used to update system atId after merge several systems
        a = []
        for index, value in df.iterrows():
            newIdx = str(int(value.atomId) + self.initIdx)
            a.append(newIdx)
        return a
    
    def changeIdx(self, x, key='atoms'):
        if key == 'atoms':
            x.atomId = str(x.atomId + 1)
        elif key == 'bonds':
            x.bondId = str(x.bondId + 1)
        return x
    
    def updateAtIdx2(self):
        self.atoms = self.atoms.reset_index(drop=True)
        self.atoms.atomId = self.atoms.index
        self.atoms = self.atoms.apply(lambda x: self.changeIdx(x), axis=1)
    
    def getNewIdx(self, x):
        idx1 = self.atoms.loc[(self.atoms.oriId == x.ai), 'atomId'].values
        idx2 = self.atoms.loc[(self.atoms.oriId == x.aj), 'atomId'].values
        # TODO: didn't work when reset charge. "Index 0 is out of bounds for axis 0 with size 0"
        x.ai = idx1[0]; x.aj = idx2[0]
        return x
    
    def updateBdIdx(self):
        df_new = self.bonds.copy(deep=True)
        df_new = df_new.reset_index(drop=True)
        df_new = df_new.apply(lambda x: self.getNewIdx(x), axis=1)
        df_new.bondId = df_new.index
        df_new = df_new.apply(lambda x: self.changeIdx(x, key='bonds'), axis=1)
        self.bonds = df_new
        
    def updateIdx(self):
        atIdx = self.updateAtIdx(self.atoms)
        self.atoms.loc[:, 'atomId'] = atIdx
        self.updateBdIdx()
    
    def name2Idx(self, lst):
        resname = lst[0]; atname = lst[1]
        idx = self.atoms[(self.atoms.resname.str.contains(resname)) & (self.atoms.atomName == atname)].atomId.to_list()
        return idx
    
    def getPos(self, atIdx):
        a = self.atoms[(self.atoms.atomId == atIdx)]
        return [a.x, a.y, a.z]
    
    def avgDis(self, lst):
        a = []
        for i in lst:
            a.append(float(i))
        
        return np.average(a)
        
    def getCOM(self, name):
        df = self.atoms[self.atoms.resname == name]
        x1 = self.avgDis(df.x)
        y1 = self.avgDis(df.y)
        z1 = self.avgDis(df.z)
        com = [x1, y1, z1]
        return com
    
    def searchHCon(self, atIdx, inKey, inHLst=[]):
        df_con = self.bonds[(self.bonds.ai == atIdx) | (self.bonds.aj == atIdx)]
        atIdx = list(set(df_con.ai.to_list() + df_con.aj.to_list()))
        con = []
        for idx in atIdx:
            if 'H' in self.atoms[(self.atoms['atomId'] == idx)].atomName.values[0]:
                con.append(idx)
        print('HCon: ', con)
        con1 = ''; atomName_ori = ''
        for idx in con:
            if idx in inHLst:
                continue
            else:
                atName = self.atoms[self.atoms.loc[:, 'atomId'] == idx]['atomName'].values[0]
                if len(con1) == 0:
                    con1 = idx
                    atomName_ori = atName
                else:
                    # take care the H sequence, start from small number
                    num0 = re.findall(r'\d+', atomName_ori)
                    num1 = re.findall(r'\d+', atName)
                    if len(num0) == 0:
                        break
                    else:
                        if len(num1) == 0:
                            con1 = idx
                            atomName_ori = atName
                        elif int(num1[0]) < int(num0[0]):
                            con1 = idx
                            atomName_ori = atName
                        elif int(num1[0]) > int(num0[0]):
                            continue
        return con1
    
    def addBonds(self, pairs):
        columns = ['bondId', 'ai', 'aj', 'type']
        for p in pairs:
            bdNum = len(self.bonds)
            bondId = str(bdNum + 1)
            ai = str(p[0]); aj = str(p[1]); atype = '1'
            tmp = pd.DataFrame([[bondId, ai, aj, atype]])
            tmp.columns = columns
            self.bonds = self.bonds.append(tmp, ignore_index=True)
    
    def shiftCoord(self, x, dist, resname, num):
        if x.resname == resname:
            x.x = float(x.x) - dist[0] - num
            x.y = float(x.y) - dist[1] - num
            x.z = float(x.z) - dist[2] - num
        return x

    def rotateMol(self, x, theta, resname):
        if x.resname == resname:
            theta = np.deg2rad(theta)
            rot = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
            v = np.array([x.x, x.y]).astype(float)
            v2 = np.dot(rot, v)
            x.x = float(v2[0])
            x.y = float(v2[1])
            x.z = x.z
        return x

    def shiftMol(self, pair, num):
        resname = self.atoms[(self.atoms.atomId == pair[0])].resname.values[0]
#        com = self.getCOM(resname);
        atPos0 = self.getPos(pair[0]) # atPos0 need to move to the atPos1
        atPos1 = self.getPos(pair[1])
        dist = [float(atPos0[0]) - float(atPos1[0]), 
                float(atPos0[1]) - float(atPos1[1]),
                float(atPos0[2]) - float(atPos1[2])]
        self.atoms = self.atoms.apply(lambda x: self.shiftCoord(x, dist, resname, num), axis=1)
        self.atoms = self.atoms.apply(lambda x: self.rotateMol(x, 30, resname), axis=1)

    def getSeq(self):
        df = self.atoms[(self.atoms.resname.str.contains(self.mainResname))]
        self.seq = df.atomName.to_list()
        
    def genBonds(self, inConInfo): 
        monH = []
        croH = []
        newBonds = []
        usedAtoms = []
        for conInfo in inConInfo:
            at1 = conInfo[0]; at2 = conInfo[1] # name of the atoms
            self.mainResname = conInfo[1][0]
            monlst = self.name2Idx(at1)
            crolst = self.name2Idx(at2) 
            print('monlst: ', monlst)
            for i in monlst:
                hatom = self.searchHCon(i, at1[0])
                monH.append(hatom)
            
            for i in crolst:
                hatom = self.searchHCon(i, at2[0], croH)
                croH.append(hatom)
            
            for i in range(len(crolst)):
                at1 = crolst[i]
                idx = 0
                at2 = monlst[idx]
                while(at2 in usedAtoms and idx < len(monlst)):
                    idx += 1
                    at2 = monlst[idx]
                    
                newBonds.append([at1, at2])
                usedAtoms.append(monlst[idx])
        print('monH: ', monH)
        print('croH: ', croH)
        print('newBonds: ', newBonds)
        self.addBonds(newBonds)

        # shift molecules to right position
        num = 3
        for pair in newBonds:
            self.shiftMol(pair, num)
            num += 1
            
        for i in monH:
            self.atoms.drop(self.atoms[self.atoms['atomId'] == i].index, inplace=True)
            self.bonds.drop(self.bonds[self.bonds['ai'] == i].index, inplace=True)
            self.bonds.drop(self.bonds[self.bonds['aj'] == i].index, inplace=True)
        
        for i in croH:
            self.atoms.drop(self.atoms[self.atoms['atomId'] == i].index, inplace=True)
            self.bonds.drop(self.bonds[self.bonds['ai'] == i].index, inplace=True)
            self.bonds.drop(self.bonds[self.bonds['aj'] == i].index, inplace=True)
            
        self.updateAtIdx2()
        self.updateBdIdx()
        self.basicInfo['atnum'] = str(len(self.atoms))
        self.basicInfo['bondnum'] = str(len(self.bonds))
        self.getSeq()

    def mergeRow(self, x, key='atoms'):
        if key == 'atoms':
            outStr = '{:>7} {:<4}{:>15}{:>10}{:>10}{:>5}{:>5}{:>5}{:>14}'.format(x.atomId, x.atomName, round(float(x.x), 2), round(float(x.y), 2), round(float(x.z), 2), 
                                                                     x.atype, x.set, x.resname, x.charge)
        elif key == 'bonds':
            outStr = '{:>6}{:>6}{:>6}{:>5}'.format(x.bondId, x.ai, x.aj, x.type)
        
        return outStr
    
    def addTopRow(self, df, inStr):
        df.loc[-1] = inStr
        df.index = df.index + 1
        df = df.sort_index().reset_index(drop=True)
        return df
    
    def outMol2(self, name):
        keys = ['@<TRIPOS>MOLECULE', '@<TRIPOS>ATOM', '@<TRIPOS>BOND']
        molStr0 = [self.basicInfo['name'], ' {} {} {} {} {}'.format(self.basicInfo['atnum'], self.basicInfo['bondnum'],
                                           self.basicInfo['setnum'], self.basicInfo['ssetnum'], self.basicInfo['featnum']),
                   self.basicInfo['molType'], self.basicInfo['chargeType']]
        
        atStr = self.atoms.apply(lambda x: self.mergeRow(x), axis=1)
        bdStr = self.bonds.apply(lambda x: self.mergeRow(x, key='bonds'), axis=1)
        molStr = pd.DataFrame({'0': molStr0})
        atStr = pd.DataFrame(atStr, columns=['0'])
        bdStr = pd.DataFrame(bdStr, columns=['0'])
        df_lst0 = [molStr, atStr, bdStr]
        df0 = []
        for i in range(len(keys)):
            df_tmp = self.addTopRow(df_lst0[i], keys[i])
            df0.append(df_tmp)
        
        df_out = pd.concat(df0)
        df_out.to_csv(name, mode='w', index=False, header=None)
        return df0
    