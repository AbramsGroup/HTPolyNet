# -*- coding: utf-8 -*-
"""
Searching potential bonds

@author: huang
"""
import pandas as pd
import numpy as np
from countTime import *

class searchBonds(object):
    def __init__(self, basicParameters, generatedBonds, inGro, inTop):
        self.monInfo = basicParameters.monInfo
        self.croInfo = basicParameters.croInfo
        self.cutoff = basicParameters.cutoff
        self.bondsRatio = basicParameters.bondsRatio
        self.boxSize = basicParameters.boxSize # one dimension box length
        self.maxBonds = basicParameters.maxBonds
        self.rctInfo = basicParameters.rctInfo
        self.generatedBonds = generatedBonds
        self.gro = inGro # gro object
        self.top = inTop # top object
        self.rctMon = []
        self.rctCro = []
        self.mol = []
        
    def calDist(self, aPos, bPos, pbc=True):
        boxSize = [self.boxSize, self.boxSize, self.boxSize]
        xlen = 0.5 * float(boxSize[0])
        ylen = 0.5 * float(boxSize[1])
        zlen = 0.5 * float(boxSize[2])
        
        x1 = float(aPos[0])-float(bPos[0])
        y1 = float(aPos[1])-float(bPos[1])
        z1 = float(aPos[2])-float(bPos[2])
        if pbc:
            if x1 > xlen:
                x1 =  x1 - float(boxSize[0])
            elif x1 <= -xlen:
                x1 = x1 + float(boxSize[0])
                
            elif y1 > ylen:
                y1 = y1 - float(boxSize[1])
            elif y1 <= -ylen:
                y1 = y1 + float(boxSize[1])
                
            elif z1 > zlen:
                z1 = z1 - float(boxSize[2])
            elif z1 <= -zlen:
                z1 = z1 + float(boxSize[2])
            
        x2 = np.power(x1,2)
        y2 = np.power(y1,2)
        z2 = np.power(z1,2)
    
        dist = np.sqrt(x2+y2+z2)
        return dist

    def filterRctPairs(self, x, df_mon):
        lst = []
        for index, row in df_mon.iterrows():
            if x.rctNum > 0 and row.rctNum > 0:
                dist = self.calDist([x.posX, x.posY, x.posZ], [row.posX, row.posY, row.posZ])
                if dist < self.cutoff and dist > 0.05:
                    lst.append([x.globalIdx, row.globalIdx, round(dist, 2), x.rctNum, row.rctNum, x.molCon])
        return lst

    @countTime
    def getRctDf(self):
        rctFunc = self.rctInfo
        atoms = self.gro.df_atoms
        top = self.top
        rctDf = []
        if 'molCon' in atoms.keys():
            pass
        else:
            atoms['molCon'] = '0'
        
        top.atoms['molCon'] = atoms['molCon']
        top.atoms['molNum'] = atoms['molNum']
        atLst = []
        for i in range(len(rctFunc)):
            croName = rctFunc[i][0][0]; monName = rctFunc[i][1][0]; rctPct = rctFunc[i][2][0]
            df_mon = atoms[(atoms.molName == monName) & (atoms.rct == 'True')]
            df_cro = atoms[(atoms.molName == croName) & (atoms.rct == 'True')]

            atLst.append(df_mon)
            atLst.append(df_cro)
            lst_tmp = df_cro.apply(lambda x: self.filterRctPairs(x, df_mon), axis=1).values.tolist()
            for ii in lst_tmp:
                if len(ii) == 0:
                    continue
                else:
                    for idx in ii:
                        idx.append(rctPct)
                rctDf = rctDf + ii
            
        names = ['acro', 'amon', 'dist', 'rctNum1', 'rctNum2', 'molCon', 'p']
        df = pd.DataFrame(rctDf, columns=names)
        df_atoms = pd.concat(atLst).drop_duplicates().reset_index(drop=True)
        return df_atoms, df
    
    def checkRepeat(self, lst, inLst):
        for i in lst:
            if set(i) == set(inLst):
                return False
        return True
    
    def checkCircuit(self, df_rctAtoms, row):
        mol1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molNum'].to_list()
        mol2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molNum'].to_list()
        molCon1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molCon'].to_list()[0].split(',')
        molCon2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molCon'].to_list()[0].split(',')
        
        if mol1[0] in molCon2 or mol2[0] in molCon1 or mol1[0] == mol2[0]:
            return False
        else:
            return True
    
    def checkAtomsRepeat(self, lst, inAtoms):
        for a in inAtoms:
            if a in lst:
                return False
            else:
                continue
        return True

    def updateMolCon(self, atomsDf, a1, a2):
        mol1 = atomsDf[atomsDf.globalIdx == str(a1)].molNum.values[0]
        mol2 = atomsDf[atomsDf.globalIdx == str(a2)].molNum.values[0]
        
        con1 = list(atomsDf[atomsDf.molNum == mol1].molCon)[0]
        con2 = list(atomsDf[atomsDf.molNum == mol2].molCon)[0]
        
        tmp = '{},{}'.format(con1, mol2)
        atomsDf.loc[atomsDf['molNum'] == mol1, 'molCon'] = tmp
        
        tmp = '{},{}'.format(con2, mol1)
        atomsDf.loc[atomsDf['molNum'] == mol2, 'molCon'] = tmp
    
    def finalRctPairs(self, df_rctAtoms, df_pairs):
        # Sort by distance between potential atoms
        df_pairs = df_pairs.sort_values(by=['dist'])
        atomsDf = self.gro.df_atoms
        # Remove repeat rows
        lst = []; rowList = []
        for index, row in df_pairs.iterrows():
            if self.checkRepeat(lst, [row.acro, row.amon]):
                lst.append([row.acro, row.amon]); rowList.append(row)
            else:
                continue
        df_tmp1 = pd.DataFrame(rowList) 
    
        # check circuit connection. Molecules cannot connect to the same molecules
        rowList = []; atomsList = []
        for index, row in df_tmp1.iterrows():
            if self.checkCircuit(atomsDf, row):
                if self.checkAtomsRepeat(atomsList, [row.acro, row.amon]):
                    atomsList.append(row.acro); atomsList.append(row.amon); rowList.append(row)
                    self.updateMolCon(atomsDf, row.acro, row.amon)
            else:
                continue
        df_tmp2 = pd.DataFrame(rowList)
        
        return df_tmp2
    
    def idx2Mol(self, pairs):
        # pairs is a dataframe
        atoms = []; mols = []
        for index, row in pairs.iterrows():
            atoms.append(row.acro)
            atoms.append(row.amon)
        
        groDf = self.gro.df_atoms
        print('atoms: ', atoms)
        for a in atoms:
            index = list(groDf[groDf.globalIdx == str(a)].molNum)[0]
            mols.append(index)
        
        self.mol = mols
    
    @countTime
    def main(self):
        pairs = []; count = 10
        while(len(pairs) == 0):
            for i in range(count): # max trial times
                df_rctAtoms, df_pairs = self.getRctDf()
                if len(df_pairs) > 0:
                    break
            a1 = self.finalRctPairs(df_rctAtoms, df_pairs)
            if len(a1) == 0:
                print(self.cutoff)
                self.cutoff += 0.5
                if self.cutoff > 0.5 * float(self.boxSize):
                    break
            pairs = a1
        
        self.idx2Mol(pairs)
        return a1, self.mol
        
        