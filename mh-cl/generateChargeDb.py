# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:34:00 2020

@author: huang
"""

import readParameters
import createRctMol
import readMol
import subprocess
import os
import glob
from decimal import Decimal

class generateChargeDb(object):
    def __init__(self):
        self.srcPath = os.getcwd()
        self.unrctPath = ''
        self.rctPath = ''
        self.basicPath = self.srcPath + '/basic'
        self.charge = {}
        
    def createTemplate(self, rctTimes):
        a = readParameters.parameters()
        a.setName('{}/options.txt'.format(self.basicPath))
        a.readParam()
        mainMol = createRctMol.createRctMol()
        mainMol.getRctInfo(a) # assume only contain 1 cro
        keys = mainMol.getKeys('mon') + mainMol.getKeys('cro')
        croName = '{}{}.mol2'.format(self.unrctPath, list(mainMol.rctInfo['cro'][0].keys())[0])
        mainMol.croResName = 'CRO'
        
        croMol = readMol.readMol(croName)
        croMol.main()
        croMol.mol2.updateResName(mainMol.croResName)
        molList = {}
        for i in range(len(keys)):
            name = '{}{}.mol2'.format(self.unrctPath, keys[i])
            monMol = readMol.readMol(name)
            monMol.main()
            mainMol.base = monMol.mol2
            mainMol.connect = croMol.mol2
            for ii in range(rctTimes):
                unrctMol = mainMol.mergeMol(ii + 1)
                a1 = mainMol.creatMol(ii + 1, unrctMol, keys[i])
                key = '{}{}'.format(monMol.resname, ii + 1)
                molList[key] = a1
        mainMol.mol2List = molList
        mainMol.symmetry()
        mainMol.outMolLst(self.rctPath)
    
    def obtainParam(self):
        os.chdir(self.rctPath)
        fileList = glob.glob('*.mol2')
        cfileList = glob.glob('*-charge.mol2')
        if len(cfileList) > 0:
            nameList = cfileList
        else:
            nameList = []
            for f in fileList:
                name = f.split('.')[0]
                out1 = name + '-min'
                out2 = name + '-charge'
                nameList.append('{}.mol2'.format(out2))
                if os.path.isfile('{}.mol2'.format(out2)):
                    continue
                else:
                    cmd1 = 'obabel {}.mol2 -O {}.mol2 --minimize --sd --c 1e-5'.format(name, out1)
                    cmd2 = 'antechamber -j 4 -fi mol2 -fo mol2 -c gas -at gaff -i {}.mol2 -o {}.mol2 -pf Y -nc 0 -eq 1 -pl 10'.format(out1, out2)
                    subprocess.call(cmd1, shell=True)
                    subprocess.call(cmd2, shell=True)
        
        self.extractCharge(nameList) # Extract and save
        os.chdir(self.srcPath)

    def avgCharge(self, inList):
        outList = []
        cc = 0
        for l in inList:
            cc += Decimal(l)

        print('cc1: ', cc)
        if cc > 0:
            idx = inList.index(max(inList))
        else:
            idx = inList.index(min(inList))

        inList[idx] = Decimal(inList[idx]) - Decimal(cc)
        return inList

    def getSeq(self, df):
        atomType = df[(df.resname.str.contains('CRO') == False)].atype.to_list()
        charge = df[(df.resname.str.contains('CRO') == False)].charge.to_list()
        atomSeq = ''; chargeSeq = ''
        charge = self.avgCharge(charge)
        for i in range(len(atomType)):
            if i != len(atomType):
                atomSeq += '{}/'.format(atomType[i])
                chargeSeq += '{}/'.format(charge[i])
            else:
                atomSeq += '{}'.format(atomType[i])
                chargeSeq += '{}'.format(charge[i])
        return atomSeq, chargeSeq
    
    def saveMap(self, charge):
        with open('{}/charges.txt'.format(self.basicPath), 'a') as f:
            for key, values in charge.items():
                f.write('{}:{}\n'.format(key, values))
            
    def extractCharge(self, nameList):
        charge = {}
        cc = 0
        for f in nameList:
            a = readMol.readMol(f)
            a.main()
            atomsDf = a.mol2.atoms
            atomSeq, chargeSeq = self.getSeq(atomsDf)
            charge[atomSeq] = chargeSeq
            for c in chargeSeq.split('/'):
                try:
                    cc += Decimal(c)
                except:
                    continue
            print('cc: ', cc)
        self.charge = charge
        self.saveMap(charge)
    
    def main(self, unrctPath, rctPath, rctTimes):
        self.unrctPath = unrctPath
        self.rctPath = rctPath
        files = os.listdir(self.rctPath)
        if len(files) > 0:
            pass
        else:
            self.createTemplate(rctTimes)
        self.obtainParam()
        return self.charge
        
if __name__ == '__main__':
    a = generateChargeDb()
    nameList = ['systems/unrctSystem/VEB.mol2']
    a.extractCharge(nameList)
    a1 = a.charge