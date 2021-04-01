# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 13:49:19 2020

@author: huang
"""
import pandas as pd
import mol2Info
from copy import deepcopy

class readMol(object):
    def __init__(self, name):
        self.name = name
        self.resname = ''
        self.mol2 = ''
        
    def getMolInfo(self):
        df1 = pd.read_csv(self.name, names=['0'], comment=';', header=None, sep='\n', skip_blank_lines=True)
        dil_indx = list(df1.loc[df1['0'].str.startswith('@')].index)
        df_sep = []
        for i in range(len(dil_indx)):
            if i == 0:
                continue
            else:
                df_tmp = df1.iloc[dil_indx[i-1] + 1:dil_indx[i], :].reset_index(drop=True)
                if '#' in df_tmp.to_string():
                    continue
                else:
                    df_sep.append(df_tmp)
        df_sep.append(df1.iloc[dil_indx[i] + 1:, :])
        return df_sep
    
    def splitRow(self, row, idx=0):
        x = deepcopy(row)
        data = list(x.str.split())[0]
        return data[idx]
    
    def extractInfo(self, df, key='basic'):
        if key == 'basic':
            basicInfo = {}
            basicInfo['name'] = df.loc[0, '0']
            basicInfo['atnum'] = df.loc[1, '0'].split()[0]
            basicInfo['bondnum'] = df.loc[1, '0'].split()[1]
            basicInfo['setnum'] = df.loc[1, '0'].split()[2]
            basicInfo['ssetnum'] = df.loc[1, '0'].split()[3]
            basicInfo['featnum'] = df.loc[1, '0'].split()[4]
            basicInfo['molType'] = df.loc[2, '0']
            basicInfo['chargeType'] = df.loc[3, '0']
            return basicInfo
            
        elif key == 'atoms':
            cNames = ['atomId', 'atomName', 'x', 'y', 'z', 'atype', 'set', 'resname', 'charge']
            df_new = pd.DataFrame(columns=cNames)
            
            for i in range(len(cNames)):
                df_new.loc[:, cNames[i]] = df.apply(lambda x: self.splitRow(x, i), axis=1)
            df_new['oriId'] = df_new['atomId']
            self.resname = df_new.loc[1, 'resname']
            return df_new
            
        elif key == 'bonds':
            cNames = ['bondId', 'ai', 'aj', 'type']
            df_new = pd.DataFrame(columns=cNames)
            
            for i in range(len(cNames)):
                df_new.loc[:, cNames[i]] = df.apply(lambda x: self.splitRow(x, i), axis=1)
            return df_new
    
    def calCOM(self, atoms):
        x = round(pd.to_numeric(atoms.x).sum()/len(atoms.x), 2)
        y = round(pd.to_numeric(atoms.y).sum()/len(atoms.y), 2)
        z = round(pd.to_numeric(atoms.z).sum()/len(atoms.z), 2)
        return [x, y, z]
        
    def main(self):
        a = self.getMolInfo()
        mol2 = mol2Info.mol2Info()
        
        info = self.extractInfo(a[0], 'basic')
        df_atoms = self.extractInfo(a[1], 'atoms')
        df_bonds = self.extractInfo(a[2], 'bonds')
        
        mol2.basicInfo = info
        mol2.atoms = df_atoms
        mol2.bonds = df_bonds
        mol2.com = self.calCOM(df_atoms)
        
        self.mol2 = mol2
        
if __name__ == '__main__':
    a = readMol('systems/unrctSystem/mST.mol2')
    a.main()
    