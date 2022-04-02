# -*- coding: utf-8 -*-
"""
Created on Sun Nov  1 08:57:03 2020

@author: huang
"""

import os
import glob
#import subprocess
import pandas as pd
#import parmed

import HTPolyNet.configuration as configuration
import HTPolyNet.readMol as readMol
import HTPolyNet.createRctMol as createRctMol
import HTPolyNet.ambertools as ambertools

class generateTypeInfo(object):
    def __init__(self,fs,cfg):
#        self.srcPath = os.getcwd()
        self.fs=fs
        self.cfg=cfg
        self.unrctPath = fs.unrctPath
        self.typePath = fs.typePath
        self.basicPath = fs.basicPath
        self.type = {}
    
    def getRctNum(self, keys, param):
        rctNum = []
        tmp_para = param.croInfo + param.monInfo
        for k in keys:
            for i in tmp_para:
                if i[1] == k:
                    rctNum.append(len(i[3]))
        return rctNum
    
    def createTemplate(self):
        if len(os.listdir(self.typePath)) > 0:
            pass
        else:
            # a = configuration.configuration()
            # a.setName('{}/options.txt'.format(self.basicPath))
            # a.readParam()
            # basicMol = a.unrctStruct
            #
            # # copy basic unreact mol2 file to the type folder
            # for molName in basicMol:
            #     n1 = '{}/{}.mol2'.format(self.unrctPath, molName)
            #     shutil.copy(n1, '{}/{}.mol2'.format(self.typePath, molName))

            mainMol = createRctMol.createRctMol()
            mainMol.getRctInfo(self.cfg) # assume only contain 1 cro
            keys = mainMol.getKeys('mon') + mainMol.getKeys('cro')
            mainMol.croResName = 'CRO'
            rctNum = self.getRctNum(keys, a)
            molList = {}
            
            for i in range(len(keys)):
                name1 = keys[i]
                molName1 = '{}/{}.mol2'.format(self.unrctPath, name1)
                mol1 = readMol.readMol(molName1)
                mol1.main()
                for name2 in mainMol.getKeys('mon'):
                    molName2 = '{}/{}.mol2'.format(self.unrctPath, name2)
                    mol2 = readMol.readMol(molName2)
                    mol2.main()
                    mol2.mol2.updateResName(mainMol.croResName)
                    mainMol.base = mol1.mol2
                    mainMol.connect = mol2.mol2
                    
                    unrctMol = mainMol.mergeMol(rctNum[i])
                    a1 = mainMol.creatMol(rctNum[i], unrctMol, name1)
                    key = '{}-{}'.format(mol1.resname, mol2.resname)
                    molList[key] = a1

            mainMol.mol2List = molList
            mainMol.outMolLst(self.typePath)
#        return molList
    
    def sepData(self, row, length, idx=0):
        data = row['0'].split()
        if len(data) != length:
            try:
                a = data[idx]
                return a
            except:
                return ' '
        else:
            return data[idx]
    
    def initSession(self, df_new, df_ori, cNames):
        i = 0
        for c in cNames:
            df_new.loc[:, c] = df_ori.apply(lambda x: self.sepData(x, len(cNames), idx=i), axis=1)
            i += 1
        return df_new
    
    def extractFF(self, nameList):
        db = {}
        btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
            
        for name in nameList:
            df1 = pd.read_csv(name, names=['0'], comment=';', header=None, sep='\n', skip_blank_lines=True)
            dil_indx = list(df1.loc[df1['0'].str.startswith('[')].index)
            df_sep = []; df_tmp = []
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

            df_bTypes = pd.DataFrame(columns=btypeNames)
            df_angTypes = pd.DataFrame(columns=angTypeNames)
            df_dihTypes = pd.DataFrame(columns=dihTypeNames)
            df_impTypes = pd.DataFrame(columns=impTypeNames)
            
            df_bTypes = self.initSession(df_bTypes, df_sep[1], btypeNames).reset_index(drop=True)
            df_angTypes = self.initSession(df_angTypes, df_sep[2], angTypeNames).reset_index(drop=True)
            df_dihTypes = self.initSession(df_dihTypes, df_sep[3], dihTypeNames).reset_index(drop=True)
            df_impTypes = self.initSession(df_impTypes, df_sep[4], impTypeNames).reset_index(drop=True)
            df_tmp = [df_bTypes, df_angTypes, df_dihTypes, df_impTypes]
            db[name] = df_tmp
            
        return db
    
    # maybe filter NaN's out
    def genFF(self, row, keys='bond'):
        if keys == 'bond':
            outKey = '{}-{}'.format(row.ai, row.aj)
            outValue = '{}, {}, {}'.format(row.funct, row.c0, row.c1)
        elif keys == 'angle':
            outKey = '{}-{}-{}'.format(row.ai, row.aj, row.ak)
            outValue = '{}, {}, {}'.format(row.funct, row.c0, row.c1)
        elif keys == 'dih':
            outKey = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
            outValue = '{}, {}, {}, {}'.format(row.funct, row.c0, row.c1, row.c2)
        elif keys == 'imp':
            outKey = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
            outValue = '{}, {}, {}, {}'.format(row.funct, row.c0, row.c1, row.c2)
        return [outKey, outValue]
    
    def filterFF(self, db):
        db_out = []; bTypes = {}; angTypes = {}; dihTypes = {}; impTypes = {}
        for keys, values in db.items():
            a1 = values[0].apply(lambda x: self.genFF(x, keys='bond'), axis=1)
            a2 = values[1].apply(lambda x: self.genFF(x, keys='angle'), axis=1)
            a3 = values[2].apply(lambda x: self.genFF(x, keys='dih'), axis=1)
            a4 = values[3].apply(lambda x: self.genFF(x, keys='imp'), axis=1)
            
            for values in a1:
                if values[0] not in bTypes.keys():
                    bTypes[values[0]] = values[1]
                    
            for values in a2:
                if values[0] not in angTypes.keys():
                    angTypes[values[0]] = values[1]
            
            for values in a3:
                if values[0] not in dihTypes.keys():
                    dihTypes[values[0]] = [values[1]]
                else:
                    if values[1] in dihTypes[values[0]]:
                        continue
                    else:
                        print('----dih dup values: ', values)
                        print('----old dih dup values: ', dihTypes[values[0]])
                        dihTypes[values[0]].append(values[1])
                        print('----new dih dup values: ', dihTypes[values[0]])

            for values in a4:
                if values[0] not in dihTypes.keys():
                    impTypes[values[0]] = values[1]
        db_out = [bTypes, angTypes, dihTypes, impTypes]
        return db_out

    # should not do this!! very bad!! parameters used should
    # be saved in a working directory, not the source!!
    # def updateParamFile(self, db):
    #     # e=ExtraGAFFParams()
    #     # newBondType = e.extra_bonds
    #     # newAngType = e.extra_angles
    #     # newDihType = e.extra_dihedrals
    #     db_BondType = db[0]
    #     db_AngType = db[1]
    #     db_DihType = db[2]
    #     db_ImpType = db[3]
        
    #     for keys, values in db_BondType.items():
    #         if keys not in newBondType:
    #             newBondType[keys] = values
        
    #     for keys, values in db_AngType.items():
    #         if keys not in newAngType:
    #             newAngType[keys] = values
            
    #     for keys, values in db_DihType.items():
    #         newDihType[keys] = values
    #     tmpList = [newBondType, newAngType, newDihType]
    #     name = ['dictBond', 'dictAngle', 'dictDihedral']
    #     if os.path.isfile('tmp.py'):
    #         os.remove('tmp.py')
    #     for i in range(len(tmpList)):
    #         with open('tmp.py', 'a') as f:
    #             print('{} = '.format(name[i]), tmpList[i], file=f)
    #     if os.path.isfile('parameters.py'):
    #         os.remove('parameters.py')
    #     os.rename('tmp.py', 'parameters.py')
        
    def obtainParam(self):
        ''' 
        1. Uses GAFF/Ambertools to generate itp file for each mol2 file
        2. Creates a dictionary keyed by {}-{}(-{}-{}) atom type -> interaction type keys and returns
        '''
        os.chdir(self.typePath)
        fileList = glob.glob('*.mol2')
        itpFile = glob.glob('*.itp')
        itpNeeded=[]
        nameList=[]
        for n in fileList:
            p=n.split('.')[0]+'.itp'
            if not p in itpFile:
                itpNeeded.append(p)
            else:
                nameList.append(p)
        A=ambertools.Parameterization()
        for f in itpNeeded:
            name = f.split('.')[0]
            print(f'---> Getting parameters for {name}.mol2...')
            # out1 = name + '-min'
            out2 = name + '-type'
            A.GAFFParameterize(f,out2,extra_antechamber_params='-eq 1 -pl 10',parmed_save_inline=False)
            # cmd1 = 'obabel {}.mol2 -O {}.mol2 --minimize --sd --c 1e-5'.format(name, out1)
            nameList.append(f'{out2}.itp')
    
        # read all itp's in and create a single merged Topology?

        db = self.extractFF(nameList)
        db = self.filterFF(db)
#        os.chdir(self.srcPath)
        # should not do this!! very bad!! parameters used should
        # be saved in a working directory, not the source!!
#        self.updateParamFile(db)
        return db  # nothing is done with this?
    
    def main(self, unrctPath, typePath):
        self.unrctPath = unrctPath
        self.typePath = typePath
        self.createTemplate()
        self.obtainParam()
        
if __name__ == '__main__':
    a = generateTypeInfo()
    srcPath = os.getcwd()
    a1 = a.main('{}/systems/unrctSystem/'.format(srcPath), '{}/systems/tmp/'.format(srcPath))
    