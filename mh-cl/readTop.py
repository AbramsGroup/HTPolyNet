# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 14:48:12 2020
This function will be replaced by readTop-2 script
@author: huangming
"""
import pandas as pd
import sys
from countTime import *

class initTop(object):
    def __init__(self):
        self.topName = ''
        self.itpName = ''
        self.aTypes = ''
        self.mTypes = ''
        self.bTypes = ''
        self.angTypes = ''
        self.dihTypes = ''
        self.impTypes = ''
        self.atoms = ''
        self.bonds = ''
        self.pairs = ''
        self.angles = ''
        self.dihs = ''
        self.imps = ''
#        self.system = ''
#        self.molecules = ''
        self.dupDihTypeKey = []
        self.topInfo = ''
        self.sumTop = []
        
    def setName(self, name1, name2):
        self.topName = name1
        self.itpName = name2
    
    @countTime
    def getTopInfo(self, name):
        print('name: ', name)
        df1 = pd.read_csv(name, names=['0'], comment=';', header=None, sep='\n', skip_blank_lines=True)
        print(df1)
        return df1

        dil_indx = list(df1.loc[df1['0'].str.startswith('[')].index)
        df_sep = []
        for i in range(len(dil_indx)):
            if i == 0:
                continue
            else:
                df_tmp = df1.iloc[dil_indx[i-1] + 1:dil_indx[i], :]
                if '#' in df_tmp.to_string():
                    continue
                else:
                    df_sep.append(df_tmp)
        df_sep.append(df1.iloc[dil_indx[i] + 1:, :])
        return df_sep
    
    def sepData(self, row, length, idx=0):
        data = list(row.str.split())[0]
        if len(data) != length:
            try:
                a = data[idx]
                return a
            except:
                return ' '
            return row.to_string()
        else:
            return data[idx]
    
    @countTime
    def initSession(self, df_new, df_ori, cNames):
        i = 0
        for c in cNames:
            df_new.loc[:, c] = df_ori.apply(lambda x: self.sepData(x, len(cNames), idx=i), axis=1)
            i += 1
        return df_new
    
    def subAtom2Atypes(self, row, df_atoms, length=2, search=False):
        if length == 2:
            a1 = row.ai; a2 = row.aj
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']
            row.ai = a1Type; row.aj = a2Type
        elif length == 3:
            a1 = row.ai; a2 = row.aj; a3 = row.ak
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']   
            a3Type = df_atoms.loc[int(a3)-1, 'type']  
            row.ai = a1Type; row.aj = a2Type; row.ak = a3Type
        elif length == 4:
            a1 = row.ai; a2 = row.aj; a3 = row.ak; a4 = row.al
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']   
            a3Type = df_atoms.loc[int(a3)-1, 'type']  
            a4Type = df_atoms.loc[int(a4)-1, 'type']
            if search:
                return [a1Type, a2Type, a3Type, a4Type]
            else:
                row.ai = a1Type; row.aj = a2Type; row.ak = a3Type; row.al = a4Type
        return row

    def extractType(self, df, df_atoms, keys='bonds'):
        if keys == 'bonds':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=2), axis=1)
        elif keys == 'angles':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=3), axis=1)
        elif keys == 'dih':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=4), axis=1)        
        return df_new
    
    def rmDihType(self, df, dupKey=False): # Since some dih type are mulitple defined, they don't need to add to the dih type section. Detail parameters are list in the dih section
        keys = []; dup_keys = []
        if dupKey:
            '''
            In this option, passed df is the df_dih
            '''
            for index, row in df.iterrows():
                if row.c0 != ' ':
                    a1Type, a2Type, a3Type, a4Type = self.subAtom2Atypes(row, self.atoms, length=4, search=True)
                    key = '{}-{}-{}-{}'.format(a1Type, a2Type, a3Type, a4Type)
                    dup_keys.append(key)
            self.dupDihTypeKey = dup_keys
        else:
            '''
            In this option, passed df is the df_dihTypes
            '''
            for index, row in df.iterrows():
                key = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
                if key in keys:
                    dup_keys.append(key)
                else:
                    keys.append(key)
            self.dupDihTypeKey = dup_keys
    
            for index, row in df.iterrows():
                key = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
                if key in dup_keys:
                    df.drop(index, inplace=True)    
        return df

    @countTime        
    def genTopSession(self):
#        df_lst0 = self.getTopInfo(self.topName)
        df_lst = self.getTopInfo(self.itpName)
        return df_lst
        atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        moltypeNames = ['name', 'nrexcl']
        atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
        bNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        pNames = ['ai', 'aj', 'funct']
        angNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        
        df_atypes = pd.DataFrame(columns=atypeNames)
        df_mtypes = pd.DataFrame(columns=moltypeNames)
        df_atoms = pd.DataFrame(columns=atNames)
        df_bonds = pd.DataFrame(columns=bNames)
        df_pairs = pd.DataFrame(columns=pNames)
        df_angles= pd.DataFrame(columns=angNames)
        df_dih = pd.DataFrame(columns=dihNames)
        df_imp = pd.DataFrame(columns=impNames)
        
        if len(df_lst) == 7:
            df_atypes = self.initSession(df_atypes, df_lst[0], atypeNames).reset_index(drop=True)
            df_mtypes = self.initSession(df_mtypes, df_lst[1], moltypeNames).reset_index(drop=True)
            df_atoms = self.initSession(df_atoms, df_lst[2], atNames).reset_index(drop=True)
            df_bonds = self.initSession(df_bonds, df_lst[3], bNames).reset_index(drop=True)
            df_pairs = self.initSession(df_pairs, df_lst[4], pNames).reset_index(drop=True)
            df_angles= self.initSession(df_angles, df_lst[5], angNames).reset_index(drop=True)
            df_dih = self.initSession(df_dih, df_lst[6], dihNames).reset_index(drop=True)
        
            df_bTypes = self.extractType(df_bonds, df_atoms, keys='bonds').drop_duplicates().reset_index(drop=True)
            df_angTypes = self.extractType(df_angles, df_atoms, keys='angles').drop_duplicates().reset_index(drop=True)
            df_dihTypes = self.extractType(df_dih, df_atoms, keys='dih').drop_duplicates().reset_index(drop=True)
            df_dihTypes = self.rmDihType(df_dihTypes)
        
            self.aTypes = df_atypes
            self.mTypes = df_mtypes
            self.bTypes = df_bTypes
            self.angTypes = df_angTypes
            self.dihTypes = df_dihTypes
            self.impTypes = pd.DataFrame(impNames)
            self.atoms = df_atoms
            self.bonds = df_bonds
            self.pairs = df_pairs
            self.angles = df_angles
            self.dihs = df_dih
            self.imps = df_imp
            self.topInfo = df_lst0
            
        if len(df_lst) == 8:
            df_atypes = self.initSession(df_atypes, df_lst[0], atypeNames).reset_index(drop=True)
            df_mtypes = self.initSession(df_mtypes, df_lst[1], moltypeNames).reset_index(drop=True)
            df_atoms = self.initSession(df_atoms, df_lst[2], atNames).reset_index(drop=True)
            df_bonds = self.initSession(df_bonds, df_lst[3], bNames).reset_index(drop=True)
            df_pairs = self.initSession(df_pairs, df_lst[4], pNames).reset_index(drop=True)
            df_angles= self.initSession(df_angles, df_lst[5], angNames).reset_index(drop=True)
            df_dih = self.initSession(df_dih, df_lst[6], dihNames).reset_index(drop=True)
            df_imp = self.initSession(df_imp, df_lst[7], dihNames).reset_index(drop=True)
        
            df_bTypes = self.extractType(df_bonds, df_atoms, keys='bonds').drop_duplicates().reset_index(drop=True)
            df_angTypes = self.extractType(df_angles, df_atoms, keys='angles').drop_duplicates().reset_index(drop=True)
            df_dihTypes = self.extractType(df_dih, df_atoms, keys='dih').drop_duplicates().reset_index(drop=True)
            df_impTypes = self.extractType(df_imp, df_atoms, keys='dih').drop_duplicates().reset_index(drop=True)
            df_dihTypes = self.rmDihType(df_dihTypes)
        
            self.aTypes = df_atypes
            self.mTypes = df_mtypes
            self.bTypes = df_bTypes
            self.angTypes = df_angTypes
            self.dihTypes = df_dihTypes
            self.impTypes = df_impTypes
            self.atoms = df_atoms
            self.bonds = df_bonds
            self.pairs = df_pairs
            self.angles = df_angles
            self.dihs = df_dih
            self.imps = df_imp
            self.topInfo = df_lst0
            
        elif len(df_lst) == 12: # already extract the corresponding types
            btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
            angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
            dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
            impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']

            df_bTypes = pd.DataFrame(columns=btypeNames)
            df_angTypes = pd.DataFrame(columns=angTypeNames)
            df_dihTypes = pd.DataFrame(columns=dihTypeNames)
            df_impTypes = pd.DataFrame(columns=impTypeNames)
            df_atypes = self.initSession(df_atypes, df_lst[0], atypeNames).reset_index(drop=True)
            df_bTypes = self.initSession(df_bTypes, df_lst[1], btypeNames).reset_index(drop=True)
            df_angTypes = self.initSession(df_angTypes, df_lst[2], angTypeNames).reset_index(drop=True)
            df_dihTypes = self.initSession(df_dihTypes, df_lst[3], dihTypeNames).reset_index(drop=True)
            df_impTypes = self.initSession(df_impTypes, df_lst[4], impTypeNames).reset_index(drop=True)
            df_mtypes = self.initSession(df_mtypes, df_lst[5], moltypeNames).reset_index(drop=True)
            df_atoms = self.initSession(df_atoms, df_lst[6], atNames).reset_index(drop=True)
            df_bonds = self.initSession(df_bonds, df_lst[7], bNames).reset_index(drop=True)
            df_pairs = self.initSession(df_pairs, df_lst[8], pNames).reset_index(drop=True)
            df_angles= self.initSession(df_angles, df_lst[9], angNames).reset_index(drop=True)
            df_dih = self.initSession(df_dih, df_lst[10], dihNames).reset_index(drop=True)
            df_imp = self.initSession(df_imp, df_lst[11], dihNames).reset_index(drop=True)
            self.aTypes = df_atypes
            self.mTypes = df_mtypes
            self.bTypes = df_bTypes
            self.angTypes = df_angTypes
            self.dihTypes = df_dihTypes
            self.impTypes = df_impTypes
            self.atoms = df_atoms
            self.bonds = df_bonds
            self.pairs = df_pairs
            self.angles = df_angles
            self.dihs = df_dih
            self.imps = df_imp
            self.topInfo = df_lst0
            self.rmDihType(df_dih, dupKey=True)

        self.sumTop = [df_lst0, self.aTypes, self.mTypes, self.bTypes, self.angTypes, self.dihTypes, self.impTypes,
                       self.atoms, self.bonds, self.pairs, self.angles, self.dihs, self.imps, self.dupDihTypeKey]
        
#    def addRctInfo(df_atoms, rctLst, rctNum):
#        df_atoms['rct']= df_atoms['nr'].apply(lambda x: True if x in rctLst else False)
#        df_atoms['rctNum'] = 0
#        for i in range(len(rctLst)):
#            r = rctLst[i]; rN = rctNum[i]
#            df_atoms['rctNum'] = df_atoms.apply(lambda x: rN if x.nr == r else x.rctNum, axis=1)
#        return df_atoms
    
if __name__ == '__main__':
    a = initTop()
#    a.setName('systems/VEA.top', 'systems/VEA.itp')
    a.setName('init.top', 'init.itp')
    a1 = a.genTopSession()  
#    import topInfo
#    b = topInfo.top()
#    b.setInfo(a.sumTop)
#    b.outDf('tmp-1')
    # b.checkCharge()
    # c = b.outDf('tmp-1.top')