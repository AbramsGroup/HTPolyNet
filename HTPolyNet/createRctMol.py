# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 14:25:35 2020

@author: huang
"""
import HTPolyNet.readMol as readMol
import HTPolyNet.mol2Info as mol2Info

import pandas as pd
import itertools
import os

class createRctMol(object):
    def __init__(self):
        self.base = []
        self.connect = ''
        self.mol2List = []
        self.rctInfo = []
        self.croResName = 'CRO'
        
    def genMol(self, inObj, startIdx, resname):
        a = mol2Info.mol2Info()
        a.basicInfo = inObj.basicInfo
        a.atoms = inObj.atoms.copy(deep=True)
        a.bonds = inObj.bonds.copy(deep=True)
        a.com = inObj.com
        a.initIdx = startIdx
        a.atoms.loc[:, 'resname'] = resname
        a.updateIdx()
        
        return a
    
    def mergeAtoms(self, inLst):
        df = pd.concat(inLst)
        df = df.reset_index(drop=True)
        idx = []
        for i in df.index:
            idx.append(str(i + 1))
        df.oriId = idx
        df.atomId = idx
        
        return df
    
    def mergeBonds(self, inLst):
        df = pd.concat(inLst)
        df = df.reset_index(drop=True)
        idx = []
        for i in df.index:
            idx.append(str(i + 1))
        df.bondId = idx
        return df
    
    def mergeMol(self, num): # num is the number of the connect
        mol2List = [self.base]
        startIdx = int(self.base.atoms.atomId.to_list()[-1])
        
        for i in range(num):
            resname = self.connect.atoms.resname + str(i)
            tmp = self.genMol(self.connect, startIdx, resname) # idx start from 1 and first is the base mol
            startIdx += int(len(self.connect.atoms))
            mol2List.append(tmp)
        
        outMol = mol2Info.mol2Info()
        info = {}
        info['name'] = mol2List[0].basicInfo['name']
        info['setnum'] = mol2List[0].basicInfo['setnum']
        info['ssetnum'] = mol2List[0].basicInfo['ssetnum']
        info['featnum'] = mol2List[0].basicInfo['featnum']
        info['molType'] = mol2List[0].basicInfo['molType']
        info['chargeType'] = mol2List[0].basicInfo['chargeType']
        
        atNum = 0; bondsNum = 0; atDf = []; bondsDf = []
        for i in mol2List:
            atNum += int(i.basicInfo['atnum'])
            bondsNum += int(i.basicInfo['bondnum'])
            atDf.append(i.atoms); bondsDf.append(i.bonds)
        
        info['atnum'] = atNum; info['bondnum'] = bondsNum
        outMol.basicInfo = info
        outMol.atoms = self.mergeAtoms(atDf)
        outMol.bonds = self.mergeBonds(bondsDf)
        return outMol
    
    def getRctInfo(self, paramObj):
        monR_list = paramObj.monR_list
        croR_list = paramObj.croR_list
        rctInfo = {'mon': [], 'cro': []}
        for i in monR_list:
            tmp = {}; tmp[i] = []
            rctTimes = 0
            for ii in monR_list[i]:
                for idx in range(int(ii[1])):
                    tmp[i].append(ii[0])
                rctTimes += int(ii[1])
            tmp[i].append(rctTimes)
            rctInfo['mon'].append(tmp)

        for i in croR_list:
            tmp = {}; tmp[i] = []
            rctTimes = 0
            for ii in croR_list[i]:
                for idx in range(int(ii[1])):
                    tmp[i].append(ii[0])
                rctTimes += int(ii[1])
            tmp[i].append(rctTimes)
            rctInfo['cro'].append(tmp)

        self.rctInfo = rctInfo

    def creatMol(self, rctTimes, molObj, inKey, croKey='mon'):
        # get cro mol
        # print('croKey: ', croKey)
        for key, value in self.rctInfo[croKey][0].items(): # sort of hard code, need to take care
            at1 = [self.croResName, value[0]]
        
        # get mon mol
        for mol in self.rctInfo['mon'] + self.rctInfo['cro']:
            for key, value in mol.items():
                if key == inKey:
                    at2 = [key, value]

        mol2List = []
        a = list(itertools.combinations(at2[1][:-1], rctTimes))
        for ii in a:
            conInfo = []
            for info in ii:
                conInfo.append([at1, [inKey, info]])
            # print('conInfo: ', conInfo)
            a = mol2Info.mol2Info()
            a.setInfo(molObj)
            a.genBonds(conInfo)
            mol2List.append(a)
        return mol2List
    
    def getKeys(self, k):
        keys = []
        for i in self.rctInfo[k]:
            keys.append(list(i.keys())[0])
        
        return keys
    
    def symmetry(self):
        seq = []
        molList = {}
        for keys, mols in self.mol2List.items():
            molList[keys] = []
            for m in mols:
                if m.seq not in seq or m.seq.reverse() not in seq:
                    seq.append(m.seq)
                    molList[keys].append(m)
                else:
                    pass
        self.mol2List = molList
    
    def outMolLst(self, path):
        idx = 0
        for keys, value in self.mol2List.items():
            name = os.path.join(path, ''.join([i for i in keys if not i.isdigit()]))
            # name = path + ''.join([i for i in keys if not i.isdigit()])
            for i in range(len(value)):
                n = name + str(idx) + '.mol2'
                value[i].outMol2(n)
                idx += 1

# if __name__ == '__main__':
#     path = '../systems/unrctSystem/'
#     rctPath = '../systems/tmp/'
#     # Get connect info from options.txt
#     import HTPolyNet.configuration as configuration

#     a = configuration.parameters()
#     a.setName('../basic/options.txt')
#     a.readParam()
#     mainMol = createRctMol()
#     mainMol.getRctInfo(a) # assume only contain 1 cro
#     keys = mainMol.getKeys('mon') + mainMol.getKeys('cro')
#     print('keys: ', keys)
#     croName = '{}{}.mol2'.format(path, list(mainMol.rctInfo['cro'][0].keys())[0])
#     mainMol.croResName = 'CRO'

#     croMol = readMol.readMol(croName)
#     croMol.main()
#     croMol.mol2.updateResName(mainMol.croResName)
#     molList = {}

#     # generateChargeDb.py
#     # for i in range(len(keys)):
#     #     name = '{}{}.mol2'.format(path, keys[i])
#     #     monMol = readMol.readMol(name)
#     #     monMol.main()
#     #     mainMol.base = monMol.mol2
#     #     mainMol.connect = croMol.mol2
#     #     rctTimes = 1
#     #     for ii in range(rctTimes):
#     #         unrctMol = mainMol.mergeMol(ii + 1)
#     #         a1 = mainMol.creatMol(ii + 1, unrctMol, keys[i])
#     #         key = '{}{}'.format(monMol.resname, ii + 1)
#     #         molList[key] = a1
#     # mainMol.mol2List = molList
#     # mainMol.symmetry()
#     # mainMol.outMolLst(rctPath)

#     rctNum = []
#     tmp_para = a.croInfo + a.monInfo
#     for k in keys:
#         for i in tmp_para:
#             if i[1] == k:
#                 rctNum.append(len(i[3]))

#     # generateTypeInfo
#     for i in range(len(keys)):
#         name1 = keys[i]
#         molName1 = '{}/{}.mol2'.format(path, name1)
#         mol1 = readMol.readMol(molName1)
#         mol1.main()
#         for name2 in mainMol.getKeys('mon'):
#             molName2 = '{}/{}.mol2'.format(path, name2)
#             mol2 = readMol.readMol(molName2)
#             mol2.main()
#             print('mainMol.croResName: ', mainMol.croResName)
#             mol2.mol2.updateResName(mainMol.croResName)
#             mainMol.base = mol1.mol2
#             mainMol.connect = mol2.mol2

#             unrctMol = mainMol.mergeMol(rctNum[i])
#             a1 = mainMol.creatMol(rctNum[i], unrctMol, name1)
#             key = '{}-{}'.format(mol1.resname, mol2.resname)
#             molList[key] = a1
#     mainMol.mol2List = molList
#     mainMol.outMolLst(rctPath)

#     seqList = {}
#     for keys, value in molList.items():
#         seq = []
#         for m in value:
#             seq.append(m.seq)
#         seqList[keys] = seq
            
    
