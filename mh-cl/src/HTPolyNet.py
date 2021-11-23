#!/usr/bin/python3
"""
step 1: read needed parameters from basic/option.txt
step 2: init systems 
        - merge molecules 
          - Using GROMACS to insert the molecules
          - Merge top to generate a total topology dataframe
        - Export the gro and top file to the folder
        - copy the needed files to the folder (.mdp file)
        - update the coordinate in the gro df

@author: huang
"""
import os
from shutil import copyfile
from shutil import move
from shutil import rmtree
from copy import deepcopy
import subprocess
import sys
HTPATH = os.getenv('HTPOLYPATH')
sys.path.append('{}/src'.format(HTPATH))

import readParameters
import mergeTop
import readTop2
import readGro
import groInfo
import topInfo
import md
import searchBonds
import genBonds
import generateChargeDb
import generateTypeInfo
import processTop
import endCapping
import getCappingParam
from countTime import *

class main(object):
    def __init__(self, re=False):
        self.cpu = ''
        self.gpu = ''
        self.trials = ''
        self.reProject = ''
        self.stepwise = ''

        # end capping
        self.cappingBonds = [] # potential bonds for capping
        self.unrctMap = {}

        # layer limit
        self.boxLimit = 1
        self.layerConvLimit = 1

        self.HTPolyPath = HTPATH
        self.topPath = ''
        self.projPath = ''

        self.basicFolder = ''
        self.mdpFolder = ''
        self.resFolder = '' # result folder
        self.systemsFolder = ''
        self.unrctFolder = ''
        self.rctFolder = ''
        self.typeFolder = ''

        self.topMap = {}
        self.initGro = ''
        self.initTop = ''
        self.dumpPairs = {} # pair map when check circuit

        self.basicParameter = ''
        self.basicFFType = []
        self.molNames = []
        self.chargeMap = {}

        self.prevGro = ''
        self.prevTop = ''

        self.gro = '' # save gro object
        self.top = '' # save top object

        self.chains = []
        self.old_pairs = []
        self.pairs_detail = {}
        # needed parameters
        self.maxBonds = 100
        self.conv = 0
        self.desConv = 0
        self.desBonds = 0

    def initFolder(self):
        if self.reProject == '':
            i = 0
            while(os.path.isdir(os.path.join(self.topPath, 'proj{}'.format(i)))):
                i += 1
            self.projPath = os.path.join(self.topPath, 'proj{}'.format(i))
            os.mkdir('proj{}'.format(i))
        else:
            self.projPath = os.path.join(self.topPath, self.reProject)

        self.basicFolder = os.path.join(self.projPath, 'basic')
        self.mdpFolder = os.path.join(self.projPath, 'mdp')
        self.resFolder = os.path.join(self.projPath, 'results')

        self.systemsFolder = os.path.join(self.projPath, 'systems')
        self.unrctFolder = os.path.join(self.systemsFolder, 'unrctSystem')
        self.rctFolder = os.path.join(self.systemsFolder, 'rctSystem')
        self.typeFolder = os.path.join(self.systemsFolder, 'typeSystem')
        if self.reProject == '':
            os.mkdir(self.basicFolder)
            os.mkdir(self.mdpFolder)
            os.mkdir(self.resFolder)
            os.mkdir(self.systemsFolder)
            os.mkdir(self.unrctFolder)
            os.mkdir(self.rctFolder)
            os.mkdir(self.typeFolder)

            cmd1 = 'cp {}/options.txt {}/options.txt'.format(self.topPath, self.basicFolder)
            cmd2 = 'cp {}/*mol2 {}/'.format(self.topPath, self.unrctFolder)
            path = os.path.join(self.HTPolyPath, 'mdp')
            cmd3 = 'cp -r {}/* {}'.format(path, self.mdpFolder)

            subprocess.call(cmd1, shell=True)
            subprocess.call(cmd2, shell=True)
            subprocess.call(cmd3, shell=True)

        else:
            cmd0 = 'rm -r {}/*'.format(self.resFolder)
            subprocess.call(cmd0, shell=True)

    def setParam(self, name):
        a = readParameters.parameters()
        a.setName(name)
        a.readParam()
        self.basicParameter = a
        self.cpu = int(a.CPU)
        self.gpu = int(a.GPU)
        self.trials = int(a.trials)
        self.reProject = a.reProject
        self.stepwise = a.stepwise
        self.cappingBonds = a.cappingBonds
        self.boxLimit = float(a.boxLimit)
        self.layerConvLimit = float(a.layerConvLimit)

    def getGroInfo(self, name):
        a = readGro.initGro()
        a.setName(name)
        df_init, sysName, atNum, boxSize = a.readGRO()
        m = groInfo.gro()
        m.setGroInfo(df_init, sysName, atNum, boxSize)
        self.gro = m

    def getTopInfo(self, topName, itpName):
        a = readTop2.initTop()
        a.setName(topName, itpName)
        a.genTopSession()
        b = topInfo.top()
        b.setInfo(a.sumTop)
        b.checkCharge()
        return b

    @countTime
    def updateCoord(self, name):
        a = readGro.initGro()
        a.setName(name)
        df_init, sysName, atNum, boxSize = a.readGRO()
        a1 = groInfo.gro()
        a1.setGroInfo(df_init, sysName, atNum, boxSize)
        self.gro.updateCoord(a1)
        self.initGro.updateCoord(a1)
        
    def initSys(self):
        print('-> Creating mixture...')
        param = self.basicParameter
        molInfo = {}
        nameList = []
        monInfo = param.monInfo
        croInfo = param.croInfo
        for i in monInfo:
            molInfo[i[1]] = i[2]
            nameList.append(i[1])
        
        for i in croInfo:
            molInfo[i[1]] = i[2]
            nameList.append(i[1])

        os.mkdir('init'); os.chdir('init')
        copyfile('{}/npt-1.mdp'.format(self.mdpFolder), 'npt-1.mdp')
        copyfile('{}/em.mdp'.format(self.mdpFolder), 'em.mdp')

        for n in nameList:
            copyfile('{}/{}.gro'.format(self.unrctFolder, n), '{}.gro'.format(n))
            copyfile('{}/{}.top'.format(self.unrctFolder, n), '{}.top'.format(n))
            copyfile('{}/{}.itp'.format(self.unrctFolder, n), '{}.itp'.format(n))
            
        # Insert molecules to systems
        import extendSys
        a = extendSys.extendSys('gmx_mpi')
        a.extendSys(param.monInfo, param.croInfo, param.boxSize, 'init')
        
        # Get df of gro file
        self.getGroInfo('init')
        
        # Get parameters from parameters file
        topList = []
        for n in nameList:
            a = self.topMap[n]
            # a = self.getTopInfo('{}.top'.format(n), '{}.itp'.format(n))
            nNum = int(molInfo[n])
            for i in range(nNum):
                topList.append(a)
                
        # Get sum of top
        topSum = mergeTop.mergeTopList(topList)
        sysTop = topSum.outDf('init')
        self.top = topSum
        self.initTop = deepcopy(self.top)

        print('-> Successful creating mixture!!')

        # EM and NPT to equilibrate the structure
        print('-> Conduct energy minization on the new mixture')
        a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
        a.emSimulation('init', 'init', 'min-1', size=False)
        print('-> Conduct NPT on the new mixture')
        boxSize = param.boxSize
        if boxSize[0] == boxSize[1] == boxSize[2]:
            a.NPTSimulation('min-1', 'init', 'npt-1', 'npt-1', check=True, re=False)
        else:
            copyfile('{}/npt-l.mdp'.format(self.mdpFolder), 'npt-l.mdp')
            a.NPTSimulation('min-1', 'init', 'npt-l', 'npt-l', check=True, re=False)
            a.NPTSimulation('npt-l', 'init', 'npt-1', 'npt-1', check=True, re=False)
        i = 0
        # TODO: can not ensure the NPT is finished well
        print('-> The mixture is good to go!!')
        while(not a.checkMDFinish('npt-1')):
            if i > 5:
                print('Still cannot converge NPT system, restart')
                sys.exit()
            elif i == 0:
                inName = 'npt-1'
            else:
                inName = 'npt-1' + str(i)
            a.extraRun(inName, 'init', i)
            if os.path.isfile('npt.gro'):
                move('npt.gro', 'npt-1.gro')
            i += 1

        # init rct info for potential atoms, add rct columns to the df in the gro object df
        self.gro.initRctInfo(self.basicParameter)
        self.initGro = deepcopy(self.gro)

        # Update coord to the gro df
        self.updateCoord('npt-1')

        # Back to the working directory and start crosslinking approach
        os.chdir(self.resFolder)

    def countConv(self):
        num1 = 0
        for i in self.old_pairs:
            num1 += len(i)

        conv = round(num1 / int(self.maxBonds), 2)
        return conv

    def stepCapping(self):
        os.mkdir('capping')
        os.chdir('capping')
        gro = deepcopy(self.gro)
        top = deepcopy(self.top)
        a = endCapping.endCapping(gro, top, self.basicFFType, self.unrctMap, self.cappingBonds)
        gro = a.gro
        top = a.top
        top.endCappingtopClean()

        gro.outDf('sys')
        top.topClean(key='bonds')
        top.outDf('sys')
        os.chdir('..')

    def finishSim(self, folderName, conv, step=0):
        os.chdir('..')
        #os.mkdir('Final'); os.chdir('Final')

        #if conv >= self.desConv:
        #    self.top.outDf('tmp')
            # a = endCapping.endCapping(self.gro, self.top, self.basicFFType, self.unrctMap, self.cappingBonds)
            #gro = a.gro
            #top = a.top
            #top.endCappingtopClean()

            #gro.outDf('sys')
            #top.topClean(key='bonds')
            #top.outDf('sys')
        #else:
        #    if step == 0:
        #        pass
        #    else:
                #a = endCapping.endCapping(self.gro, self.top, self.basicFFType, self.unrctMap, self.cappingBonds)
                #gro = a.gro
                #top = a.top
                #top.endCappingtopClean()

                #gro.outDf('sys')
                # self.prevTop.topClean(key='bonds')
                #top.outDf('sys')

        os.chdir(self.resFolder)
        conv = self.countConv()
        move(folderName, '{}-{}'.format(folderName, conv))
        self.old_pairs = []
        self.reInitSys()
        
    def reInitSys(self):
        path = '{}/init'.format(self.resFolder)
        os.chdir(path)
        # Get df of gro file
        gro = readGro.initGro()
        top = readTop2.initTop()
        
        gro.setName('init')
        df_init, sysName, atNum, boxSize = gro.readGRO()
        atomsDf = groInfo.gro()
        atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        
        self.gro = atomsDf
        self.gro.initRctInfo(self.basicParameter)
        if os.path.isfile('npt-init.gro'):
            self.updateCoord('npt-init')
        else:
            self.updateCoord('npt-1')
        topDf = topInfo.top()
        top.setName('init.top', 'init.itp')
        top.genTopSession()
        topDf.setInfo(top.sumTop)
        # topDf.checkCharge()
        self.top = topDf
        os.chdir(self.resFolder)
    
    def setupFolder(self, idx):
        if os.path.isdir('step{}'.format(idx)):
            rmtree('step{}'.format(idx))
        os.mkdir('step{}'.format(idx))
        copyfile('{}/em.mdp'.format(self.mdpFolder), '{}/em.mdp'.format('step{}'.format(idx)))
        copyfile('{}/npt-cl.mdp'.format(self.mdpFolder), '{}/npt-cl.mdp'.format('step{}'.format(idx)))
        copyfile('{}/npt-sw.mdp'.format(self.mdpFolder), '{}/npt-sw.mdp'.format('step{}'.format(idx)))
        return 'step{}'.format(idx)
    
    def calMaxBonds(self):
        maxRct = 0
        for i in self.basicParameter.monInfo:
            molNum = int(i[2])
            tmp = 0
            for ii in i[3]:
                tmp += int(ii[1])
            maxRct += molNum * tmp
            
        for i in self.basicParameter.croInfo:
            molNum = int(i[2])
            tmp = 0
            for ii in i[3]:
                tmp += int(ii[1])
            maxRct += molNum * tmp
                
        self.maxBonds = maxRct * 0.5
        print('-> Maximum number of the potential bonds is {}'.format(self.maxBonds))
        self.desConv = float(self.basicParameter.bondsRatio)
        self.desBonds = float(self.basicParameter.bondsRatio) * self.maxBonds

    def logBonds(self, step, cutoff):
        num1 = 0
        for i in self.old_pairs:
            num1 += len(i)
#        num1 = len(self.old_pairs)
        conv = num1/self.maxBonds
        self.conv = conv

        with open('../bond.txt', 'a') as f1:
#            str1 = 'step {} generate {} bonds. {} bonds left. Reach conversion {:.2f}\n'.format(step, 
#                         num1, self.maxBonds - num1, conv)
            str1 = 'step {}: {} bonds are formed, within cutoff {}nm. {} bonds left. Reach conversion {:.2f}\n'.format(step,
                         len(self.old_pairs[int(step)]), round(cutoff, 2), self.maxBonds - num1, conv)
            f1.write(str1)
        
        with open('../bonds_connect_Info.txt', 'w') as f2:
            str0 = 'Total bonds: {}\n'.format(self.maxBonds)
            # values = self.pairs_detail['step{}'.format(step)]
            f2.write(str0)
            for keys, values in self.pairs_detail.items():
                f2.write('{}: \n'.format(keys))
                for index, row in values.iterrows():
                    f2.write('atom1: {}\t{}\t atom2: {}\t{}\t mol1: {}\t mol2: {}\n'.format(
                        row.amon, row.monAtomName, row.acro, row.croAtomName, row.monMol, row.croMol))

    @countTime
    def stepwiseRelax(self):
        print('----> Start stepwise relaxation step to relax the system')
        k = self.stepwise
        outName = 'sw'
        for i in range(len(k)):
            groName = outName
            topName = '{}-{}'.format(outName, i)
            self.gro.outDf(groName)
            self.top.outDf(topName, float(k[i]), simple=False, stepRelax=True)
            a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
            cond0 = a.emSimulation(groName, topName, 'sw-min-{}'.format(i), size=False, check=False)
            if cond0 == False:
                print('----> Stepwised EM failed')
                return False

            cond1 = a.NPTSimulation('sw-min-{}'.format(i), topName,
                                    'sw-npt-{}'.format(i), 'npt-sw',
                                    check=False, re=True)
            if cond1 == False:
                print('----> Stepwised NPT failed')
                return False

            print('----> Stepwised step for k = {} is succuessful'.format(k[i]))
            self.updateCoord('sw-npt-{}'.format(i))

        return True

    def mainProcess(self, repeatTimes):
        # Init systems
        os.chdir(self.resFolder)
        self.initSys()
        conv = 0
        # calculate max bonds
        self.calMaxBonds()

        # Start crosslinking approach
        step = 0        
        for i in range(repeatTimes):
            print('--> Start crosslinking procedure on replica {}!!'.format(i))
            folderName = 'sim{}'.format(i)
            os.mkdir(folderName); os.chdir(folderName)
            print('---> A new NPT to mix the system at 300K to create new replica')
            os.mkdir('init'); os.chdir('init')

            self.top.outDf('init')
            self.gro.outDf('init')

            copyfile('{}/npt-init.mdp'.format(self.mdpFolder), 'npt-init.mdp')
            copyfile('{}/em.mdp'.format(self.mdpFolder), 'em.mdp')
            a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
            a.emSimulation('init', 'init', 'min-1', size=False)
            a.NPTSimulation('min-1', 'init', 'npt-init', 'npt-init', check=False, re=False)
            
            # complete structs
            tmpGro = readGro.initGro()
            tmpGro.setName('npt-init')
            tmp_df, tmp_name, tmp_atNum, tmp_boxSize = tmpGro.readGRO()
            cmd1 = f'echo 0 0 | gmx_mpi trjconv -s npt-init.tpr -f npt-init.gro -o npt-init.gro -center -pbc atom -box {tmp_boxSize}'
            subprocess.call(cmd1, shell=True)
            self.updateCoord('npt-init')
            os.chdir('..')
            print('---> New replica is good to go')
            boxLimit = self.boxLimit # setup layer curing condition
            while(len(self.old_pairs) < int(self.maxBonds)):
                print('---> step {}'.format(step))
                folderName1 = self.setupFolder(step)
                os.chdir(folderName1)
                print('     (Content can be found under folder {})'.format(os.getcwd()))
                # searching potential bonds
                print('----> Start searching bonds')
                if conv < boxLimit * self.layerConvLimit:
                    boxLimit = self.boxLimit
                else:
                    print('1st layer conversion reached desired {} conversion'.format(self.layerConvLimit))
                    boxLimit = 1
                    self.boxLimit = 1
                sbonds = searchBonds.searchBonds(self.cpu, self.basicParameter, self.old_pairs, self.gro, self.top,
                                                    self.conv, self.desBonds, self.chains, boxLimit)
                pairs, chains, rMols, cutoff = sbonds.sBonds()
                # intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']
                self.chains = chains
                if len(pairs) > 0:
                    print('----> Start generating bonds')
                    self.pairs_detail['step{}'.format(step)] = pairs

                    # generate bonds
                    gbonds = genBonds.genBonds(self.gro, self.top, pairs, self.chargeMap, rMols, cat='map')
                    gbonds.gBonds() # update atom's rct status

                    self.gro = gbonds.gro
                    self.top = gbonds.top
                    self.top.checkCharge()

                    cond = self.stepwiseRelax()
                    if cond == False:
                        print('----> Stepwised cannot relax the systems, some wired bonds may formed in the step, will start a new replica')
                        self.finishSim(folderName, 0, step=step)
                        step = 0
                        break

                    groName = 'cl-{}'.format(i); topName = 'init'
                    self.gro.outDf(groName)
                    self.top.outDf(topName)

                    # intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']

                    # Equilibrate system
                    print('----> Energy minimization on the normal system')
                    a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
                    cond0 = a.emSimulation(groName, topName, 'min-1', size=False, check=False)
                    if cond0 == False:
                        print('EM failed')
                        self.finishSim(folderName, 0, step=step)
                        step = 0
                        break
                    print('----> NPT on the normal system')
                    cond1 = a.NPTSimulation('min-1', topName, 'npt-cl', 'npt-cl', check=False, re=True)
                    if cond1 == False:
                        print('NPT failed')
                        self.finishSim(folderName, 0, step=step)
                        step = 0
                        break
                    
                    # Update coord
                    self.updateCoord('npt-cl')
                    self.old_pairs.append(pairs)
                    self.logBonds(step, cutoff)
                    # self.stepCapping()
                    conv = self.countConv()
                    print('----> step {} reaches {} conversion'.format(step, round(conv, 2)))
                    if conv >= self.desConv:
                        self.finishSim(folderName, conv, step=step)
                        step = 0
                        break

                    self.prevGro = deepcopy(self.gro)
                    self.prevTop = deepcopy(self.top)
                    os.chdir('..')
                    step += 1
                else:
                    self.finishSim(folderName, conv, step=step)
                    step = 0
                    break

            self.gro = deepcopy(self.initGro)
            self.top = deepcopy(self.initTop)

    def getMolNames(self):
        names = []
        for n in self.basicParameter.unrctStruct:
            names.append(n)

        for n in self.basicParameter.monInfo:
            names.append(n[1])
        
        for n in self.basicParameter.croInfo:
            names.append(n[1])
        
        self.molNames = names
    
    def getChargeMaps(self, name):
        maps = {}
        with open(name, 'r') as f:
            for i in f.readlines():
                key, value = i.split(':')
                maps[key] = value
        return maps

    def getUnrctPara(self):
        import pandas as pd
        atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']

        basicName = self.basicParameter.unrctStruct
        aTypes = pd.DataFrame(columns=atypeNames)
        bTypes = pd.DataFrame(columns=btypeNames)
        angTypes = pd.DataFrame(columns=angTypeNames)
        dihTypes = pd.DataFrame(columns=dihTypeNames)
        impTypes = pd.DataFrame(columns=impTypeNames)
        atoms = pd.DataFrame(columns=atNames)

        for name in basicName:
            top = readTop2.initTop()
            top.setName('{}.top'.format(name), '{}.itp'.format(name))
            top.genTopSession()
            topSum = top.sumTop
            aTypes = aTypes.append(topSum[1], ignore_index=True)
            bTypes = bTypes.append(topSum[3], ignore_index=True)
            angTypes = angTypes.append(topSum[4], ignore_index=True)
            dihTypes = dihTypes.append(top.fullDihTypes, ignore_index=True)
            impTypes = impTypes.append(topSum[6], ignore_index=True)
            topSum[7].residue = name.split('-')[0]
            atoms = atoms.append(topSum[7], ignore_index=True)

        aTypes.drop_duplicates(inplace=True, ignore_index=True)
        bTypes.drop_duplicates(inplace=True, ignore_index=True)
        angTypes.drop_duplicates(inplace=True, ignore_index=True)
        dihTypes.drop_duplicates(inplace=True, ignore_index=True)
        impTypes.drop_duplicates(inplace=True, ignore_index=True)
        self.basicFFType = [aTypes, bTypes, angTypes, dihTypes, impTypes, atoms]

        unrctMap = getCappingParam.genUnrctMapping(self.basicParameter)
        self.unrctMap = unrctMap

    def preparePara(self):
        import prepareParam

        path = os.getcwd()
        self.topPath = path

        paraFile = os.path.join(path, 'options.txt')
        if os.path.isfile(paraFile):
            pass
        else:
            print('Didn\'t find configeration file. Please check!')
            sys.exit()

        self.setParam(paraFile)
        self.getMolNames()
        self.initFolder()

        os.chdir(self.unrctFolder)
        for n in self.molNames:
            fileName = '{}.mol2'.format(n)
            if os.path.isfile('{}.gro'.format(n)) and os.path.isfile('{}.top'.format(n)):
                b = processTop.processTop(n, repeat=True)
                b.main()
            else:
                a = prepareParam.prepareParam()
                a.PrepareFile(fileName, n, n)
                b = processTop.processTop(n) # process the top file, make it standard
                b.main()
            self.topMap[b.name] = b.top
        self.getUnrctPara()
        os.chdir(self.projPath)
        
        if os.path.isfile('{}/charges.txt'.format(self.basicFolder)):
            self.chargeMap = self.getChargeMaps('{}/charges.txt'.format(self.basicFolder))
        else:
            print('--> Start generating charge data base')
            a1 = generateChargeDb.generateChargeDb()
            cc = a1.main(self.unrctFolder, self.rctFolder, 4) # could be more
            
            self.chargeMap = cc
        
        os.chdir(self.projPath)
        print('--> Start generating reacted molecules type data base')
        a = generateTypeInfo.generateTypeInfo(self.topPath)
        a.main(self.unrctFolder, self.typeFolder)
        
if __name__ == '__main__':
    a = main() # change name like gmx_cl ....
    a.preparePara()
    a.mainProcess(a.trials)
    print('All replicas are been tested')
    
    # TODO: need to check that charge been update as the template. 
