# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 15:43:00 2020

@author: huang
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

    def getTopInfo(self, name):
        df1 = pd.read_csv(name, names=['0'], comment=';', header=None, sep='\n', skip_blank_lines=True)

        dil_indx = list(df1.loc[df1['0'].str.startswith('[')].index)
        df_sep = []
        for i in range(len(dil_indx)):
            if i == 0:
                continue
            else:
                df_tmp = df1.iloc[dil_indx[i - 1] + 1:dil_indx[i], :]
                if '#' in df_tmp.to_string():
                    continue
                else:
                    df_sep.append(df_tmp)
        df_sep.append(df1.iloc[dil_indx[i] + 1:, :])
        return df_sep

    @countTime
    def genTopSession(self):
        atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        moltypeNames = ['name', 'nrexcl']
        
        atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
        bNames = ['ai', 'aj', 'funct']
        pNames = ['ai', 'aj', 'funct']
        angNames = ['ai', 'aj', 'ak', 'funct']
        dihNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
    
        with open(self.itpName, 'r') as f:
            lines = f.read().split('\n')
            lst0 = []; lst_tmp = []
            for l in lines:
                if l.startswith('['):
                    lst0.append(lst_tmp)
                    lst_tmp = []
                elif l.startswith(';'):
                    continue
                else:
                    if len(l) > 1:
                        lst_tmp.append(l.split(';')[0].split())
                    
            lst0.append(lst_tmp) # last section (improper section) to the end of the list
            lst0 = lst0[1:]
        
            if len(lst0) == 12:
                names = [atypeNames, btypeNames, angTypeNames, dihTypeNames, impTypeNames, moltypeNames, 
             atNames, bNames, pNames, angNames, dihNames, impNames]
                
                info = ['aTypes', 'bTypes', 'angTypes', 'dihTypes', 'impTypes', 'mTypes',
                        'atoms', 'bonds', 'pairs', 'angles', 'dihs', 'imps']
                
                for i in range(len(lst0)):
                    if i == 11:
                        try:
                            df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                        except:
                            impNames = ['ai', 'aj', 'ak', 'al', 'funct']
                            df_tmp = pd.DataFrame(lst0[i], columns=impNames)
                    else:
                        print('names: ', names[i])
                        df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                    setattr(self, info[i], df_tmp)

            if len(lst0) == 7: # Itp file didn't contain the type section
                names = [moltypeNames, atNames, bNames, pNames, angNames, dihNames, impNames]
                info = ['mTypes', 'atoms', 'bonds', 'pairs', 'angles', 'dihs', 'imps']
                for i in range(len(lst0)):
                    if i == 6:
                        try:
                            df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                        except:
                            impNames = ['ai', 'aj', 'ak', 'al', 'funct']
                            df_tmp = pd.DataFrame(lst0[i], columns=impNames)
                    else:
                        df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                    setattr(self, info[i], df_tmp)

            elif len(lst0) == 8:
                names = [atypeNames, moltypeNames, atNames, bNames, pNames, angNames, dihNames, impNames]
                info = ['aTypes', 'mTypes', 'atoms', 'bonds', 'pairs', 'angles', 'dihs', 'imps']
                for i in range(len(lst0)):
                    if i == 7:
                        try:
                            df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                        except:
                            impNames = ['ai', 'aj', 'ak', 'al', 'funct']
                            df_tmp = pd.DataFrame(lst0[i], columns=impNames)
                    else:
                        df_tmp = pd.DataFrame(lst0[i], columns=names[i])
                    setattr(self, info[i], df_tmp)

        df_lst0 = self.getTopInfo(self.topName)
        self.sumTop = [df_lst0, self.aTypes, self.mTypes, self.bTypes,
                       self.angTypes, self.dihTypes, self.impTypes,
                       self.atoms, self.bonds, self.pairs, self.angles,
                       self.dihs, self.imps, self.dupDihTypeKey]

if __name__ == '__main__':
    a = initTop()
    a.setName('tmp.top', 'tmp.itp')
    a1 = a.genTopSession()
    
