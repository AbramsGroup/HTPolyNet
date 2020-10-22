# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 10:53:35 2020

@author: huang
"""

import pandas as pd
import os
import glob
import readTop

class molRctInfo(object):
    def __init__(self):
        self.mon = {}
        self.cro = {}
    
    def extractSeq(self, name):
        a1 = self.readMol2('{}.mol2'.format(name))
        a2 = self.readTop(name)
        a2.loc[:, 'residue'] = a1
        out = a2[a2.residue == name]
        return out.type, out.charge
    
    def getNames(self, nameList):
        data = {}
        for n in nameList:
            for files in glob.glob('{}*.top'.format(n)):
                tmp = files.strip('.top')
                aType, aCharge = self.extractSeq(tmp)
                data[tmp] = [aType.to_list(), aCharge.to_list()]
        
        return data
    
    def sepData(self, row):
        data = list(row.str.split())[0][7]
        return data
    
    def readMol2(self, name):
        df = pd.read_csv(name, names=['0'], header=None, sep='\n', skip_blank_lines=True)
        dil_indx = list(df.loc[df['0'].str.startswith('@')].index)
        df_sep = []
        for i in range(len(dil_indx)):
            if i == 0:
                continue
            else:
                df_tmp = df.iloc[dil_indx[i-1] + 1:dil_indx[i], :]
                if '#' in df_tmp.to_string():
                    continue
                else:
                    df_sep.append(df_tmp)
        df_sep.append(df.iloc[dil_indx[i] + 1:, :])
        a = []
        for index, row in df_sep[1].iterrows():
            a.append(self.sepData(row))
        
        return a
    
    def readTop(self, name):
        a = readTop.initTop()
        a.setName('{}.top'.format(name), '{}.itp'.format(name))
        a.genTopSession()
        atomsDf = a.atoms
        return atomsDf
    
    def main(self, nameList):
        a = self.getNames(nameList)
        return a
    
if __name__ == "__main__":
    os.chdir('systems/rctSystem')
    a = molRctInfo()
    nameList = ['VEA']
    b = a.main(nameList)
#    aa = a.readMol2('VEA1.mol2')
#    bb = a.readTop('VEA1')
#    bb.loc[:, 'residue'] = aa
#    cc = bb[bb.residue == 'VEA1']
#    dd = [cc.type, cc.charge]
#    