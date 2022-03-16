# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:34:00 2020

@author: huang
"""

from re import T
from typing import Type
import subprocess
import os
import glob
from decimal import Decimal

import HTPolyNet.configuration as configuration
import HTPolyNet.createRctMol as createRctMol
import HTPolyNet.readMol as readMol

class generateChargeDb(object):
    def __init__(self):
        self.srcPath = os.getcwd()
        self.unrctPath = ''
        self.rctPath = ''
        self.basicPath = self.srcPath + '/basic'
        self.charge = {}
    
    def setTemplate(self, mainMol, mainKey='mon'):
        # a = readParameters.parameters()
        # a.setName('{}/options.txt'.format(self.basicPath))
        # a.readParam()
        # mainMol = createRctMol.createRctMol()
        # mainMol.getRctInfo(a) # assume only contain 1 cro
        keys = mainMol.getKeys(mainKey)
        if mainKey == 'mon':
            tmpKey = 'cro'
        else:
            tmpKey = 'mon'

        croName = '{}/{}.mol2'.format(self.unrctPath, list(mainMol.rctInfo[tmpKey][0].keys())[0])
        mainMol.croResName = 'CRO'
        croMol = readMol.readMol(croName)
        croMol.main()
        croMol.mol2.updateResName(mainMol.croResName)
        molList = {}
        for i in range(len(keys)):
            name = '{}/{}.mol2'.format(self.unrctPath, keys[i])
            monMol = readMol.readMol(name)
            monMol.main()
            mainMol.base = monMol.mol2
            mainMol.connect = croMol.mol2
            rctTimes = mainMol.rctInfo[mainKey][0][keys[i]][-1]

            for ii in range(rctTimes):
                unrctMol = mainMol.mergeMol(ii + 1)
                a1 = mainMol.creatMol(ii + 1, unrctMol, keys[i], croKey=tmpKey)
                key = '{}{}'.format(monMol.resname, ii + 1)
                molList[key] = a1
        return molList

    def createTemplate(self, rctTimes):
        
        a = configuration.configuration()
        a.setName('{}/options.txt'.format(self.basicPath))
        a.readParam()
        mainMol = createRctMol.createRctMol()
        mainMol.getRctInfo(a) # assume only contain 1 cro
        molList = {}
        molListMon = self.setTemplate(mainMol, 'mon')
        molListCro = self.setTemplate(mainMol, 'cro')
        molList.update(molListMon)
        molList.update(molListCro)

        # keys = mainMol.getKeys('mon') + mainMol.getKeys('cro')
        # croName = '{}/{}.mol2'.format(self.unrctPath, list(mainMol.rctInfo['mon'][0].keys())[0])
        # mainMol.croResName = 'CRO'
        
        # croMol = readMol.readMol(croName)
        # croMol.main()
        # croMol.mol2.updateResName(mainMol.croResName)
        # molList = {}
        # for i in range(len(monKeys)):
        #     name = '{}/{}.mol2'.format(self.unrctPath, keys[i])
        #     monMol = readMol.readMol(name)
        #     monMol.main()
        #     mainMol.base = monMol.mol2
        #     mainMol.connect = croMol.mol2
        #     try:
        #         rctTimes = mainMol.rctInfo['mon'][0][keys[i]][-1]
        #     except:
        #         rctTimes = mainMol.rctInfo['cro'][0][keys[i]][-1]

        #     for ii in range(rctTimes):
        #         unrctMol = mainMol.mergeMol(ii + 1)
        #         a1 = mainMol.creatMol(ii + 1, unrctMol, keys[i])
        #         key = '{}{}'.format(monMol.resname, ii + 1)
        #         molList[key] = a1
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
                    # cmd1 = 'obabel {}.mol2 -O {}.mol2 --minimize --sd --c 1e-5'.format(name, out1)
                    # using this cmd on epoxy will change the force field parm
                    cmd2 = 'antechamber -j 4 -fi mol2 -fo mol2 -c gas -at gaff -i {}.mol2 -o {}.mol2 -pf Y -nc 0 -eq 1 -pl 10'.format(name, out2)
                    print('--> Getting parameters from {}.mol2...'.format(name))
                    # a1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    # out, err = a1.communicate()
                    a2 = subprocess.Popen(cmd2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                    out, err = a2.communicate()
        self.extractCharge(nameList) # Extract and save
        os.chdir(self.srcPath)

    def avgCharge(self, inList):
        outList = []
        cc = 0
        for l in inList:
            cc += Decimal(l)

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
                atomSeq += '{}/'.format(atomType[i][0])
                chargeSeq += '{}/'.format(charge[i])
            else:
                atomSeq += '{}'.format(atomType[i][0])
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
            try:
                a = readMol.readMol(f)
                a.main()
            except:
                continue

            atomsDf = a.mol2.atoms
            atomSeq, chargeSeq = self.getSeq(atomsDf)
            charge[atomSeq] = chargeSeq
            # for c in chargeSeq.split('/'):
            #     try:
            #         cc += Decimal(c)
            #     except:
            #         continue
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