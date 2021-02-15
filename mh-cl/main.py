# -*- coding: utf-8 -*-
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

import os
from shutil import copyfile
from shutil import move
from shutil import rmtree
from copy import deepcopy
import sys

class main(object):
    def __init__(self, cpu):
        self.cpu = cpu
        self.srcPath = os.getcwd()
        
        self.basicFolder = ''
        self.systemsFolder = ''
        self.unrctFolder = ''
        self.rctFolder = ''
        self.typeFolder = ''
        self.mdpFolder = ''
        self.workingFolder = '' # current loop folder location
        
        self.topMap = {}
        self.initGro = ''
        self.initTop = ''

        self.basicParameter = ''
        self.molNames = []
        self.chargeMap = {}
        self.gro = '' # save gro object
        self.top = '' # save top object
        self.old_pairs = []
        self.pairs_detail = {}
        # needed parameters
        self.maxBonds = 100
        self.conv = 0
        self.desBonds = 0

    def setParam(self, name):
        path = '{}/{}'.format(self.basicFolder, name)
        a = readParameters.parameters()
        a.setName(path)
        a.readParam()
        self.basicParameter = a
    
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
    
    def updateCoord(self, name):
        a = readGro.initGro()
        a.setName(name)
        df_init, sysName, atNum, boxSize = a.readGRO()
        a1 = groInfo.gro()
        a1.setGroInfo(df_init, sysName, atNum, boxSize)
        self.gro.updateCoord(a1)
        self.initGro = self.gro
        
    def initSys(self, ig=False):
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
        
        self.workingFolder = os.getcwd()
        if not ig:
            os.mkdir('init')
        os.chdir('init')
        copyfile('{}/npt-init.mdp'.format(self.mdpFolder), 'npt-init.mdp')
#        copyfile('{}/nvt.mdp'.format(self.mdpFolder), 'nvt-1.mdp')
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

        if not ig:
            # EM and NPT to equilibrate the structure
            a = md.md('gmx_mpi', 'mpirun', self.cpu)
            a.emSimulation('init', 'init', 'min-1', size=False)
            a.NPTSimulation('min-1', 'init', 'npt-init', 'npt-init', check=False, re=False)
            i = 0
            while(not a.checkMDFinish('npt-init')):
                if i > 5:
                    print('Still cannot converge NPT system, restart')
                    sys.exit()
                elif i == 0:
                    inName = 'npt-init'
                else:
                    inName = 'npt-init' + str(i)
                a.extraRun(inName, 'init', i)
                if os.path.isfile('npt.gro'):
                    move('npt.gro', 'npt-init.gro')
                i += 1

#        a.NVTSimulation('npt-init', 'init', 'nvt-1', 'nvt-1', check=False)
        
        # Update coord to the gro df
        self.updateCoord('npt-init')
        
        # init rct info for potential atoms, add rct columns to the df in the gro object df
        self.gro.initRctInfo(self.basicParameter)
        self.initGro = deepcopy(self.gro)
        # Back to the working directory and start crosslinking approach
        os.chdir(self.workingFolder)
    
    def finishSim(self, folderName):
        os.chdir(self.workingFolder)
        num1 = 0
        for i in self.old_pairs:
            num1 += len(i)

        conv = round(num1/int(self.maxBonds), 2)
        move(folderName, '{}-{}'.format(folderName, conv))
        self.old_pairs = []
        self.reInitSys()
        
    def reInitSys(self):
        path = '{}/init'.format(self.workingFolder)
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
        self.updateCoord('npt-init')
        topDf = topInfo.top()
        top.setName('init.top', 'init.itp')
        top.genTopSession()
        topDf.setInfo(top.sumTop)
        topDf.checkCharge()
        self.top = topDf
        os.chdir(self.workingFolder)
    
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
                         len(self.old_pairs[int(step)]), cutoff, self.maxBonds - num1, conv)
            f1.write(str1)
        
        with open('../bonds_Info{}.txt'.format(step), 'w') as f2:
            str0 = 'Total bonds: {}\n'.format(self.maxBonds)
            f2.write(str0)
            for keys, values in self.pairs_detail.items():
                f2.write('{}: \n'.format(keys))
                print('values: ', values)
                for index, row in values.iterrows():
                    f2.write('atom1: {}\tatom2: {}\n'.format(row.amon, row.acro))

    def stepwiseRelax(self):
        k = [0.001, 0.01, 0.1, 1]
        outName = 'sw'
        for i in range(len(k)):
            groName = outName
            topName = '{}-{}'.format(outName, i)
            self.gro.outDf(groName)
            self.top.outDf(topName, k[i], simple=False, stepRelax=True)
            a = md.md('gmx_mpi', 'mpirun', self.cpu)
            cond0 = a.emSimulation(groName, topName, 'sw-min-{}'.format(i), size=False, check=False)
            if cond0 == False:
                print('EM failed')
                return False

            cond1 = a.NPTSimulation('sw-min-{}'.format(i), topName,
                                    'sw-npt-{}'.format(i), 'npt-sw',
                                    check=False, re=True)
            if cond1 == False:
                print('NPT failed')
                return False

            self.updateCoord('sw-npt-{}'.format(i))

        return True

    def mainProcess(self, repeatTimes, ig=False):
        # ig is the key to do the test. if ig == True, won't create the init system
        if not ig:
            if os.path.isdir('results'):
                rmtree('results')
            else:
                pass
            os.mkdir('results'); os.chdir('results')

        else:
            os.chdir('results')
            os.system('rm -rf sim*')
        # Init systems
        self.initSys(ig)
        
        # calculate max bonds
        self.calMaxBonds()
        print('maxBonds: ', self.maxBonds)
        # Start crosslinking approach
        step = 0        
        for i in range(repeatTimes):
            folderName = 'sim{}'.format(i)
            os.mkdir(folderName)
            os.chdir(folderName)

            while(len(self.old_pairs) < int(self.maxBonds)):
                intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']
                intDf.to_csv('int-df0.csv')

                folderName1 = self.setupFolder(step)
                os.chdir(folderName1)

                # searching potential bonds
                sbonds = searchBonds.searchBonds(self.basicParameter, self.old_pairs, self.gro, self.top,
                                                 self.conv, self.desBonds)
                pairs, rMols, cutoff = sbonds.main()

                # intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']

                if len(pairs) > 0:
#                    print('pairs.amon: ', pairs.amon.values)
#                    print('pairs.amon.to_string(): ', pairs.amon.to_string())

                    self.pairs_detail['step{}'.format(step)] = pairs

                    self.gro.outDf('init') # just for check!
                    # generate bonds
                    gbonds = genBonds.genBonds(self.gro, self.top, pairs, self.chargeMap, rMols, cat='map')
                    gbonds.main() # update atom's rct status

                    self.gro = gbonds.gro
                    self.top = gbonds.top
                    self.top.checkCharge()

                    cond = self.stepwiseRelax()
                    if cond == False:
                        self.finishSim(folderName)
                        step = 0
                        break

                    groName = 'cl-{}'.format(i); topName = 'init'
                    self.gro.outDf(groName)
                    self.top.outDf(topName)

                    intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']

                    # Equilibrate system
                    a = md.md('gmx_mpi', 'mpirun', self.cpu)
                    cond0 = a.emSimulation(groName, topName, 'min-1', size=False, check=False)
                    if cond0 == False:
                        print('EM failed')
                        self.finishSim(folderName)
                        step = 0
                        break
                    
                    cond1 = a.NPTSimulation('min-1', topName, 'npt-cl', 'npt-cl', check=False, re=True)
                    if cond1 == False:
                        print('NPT failed')
                        self.finishSim(folderName)
                        step = 0
                        break
                    
                    # Update coord
                    self.updateCoord('npt-cl')
                    step += 1
                    os.chdir('..')
                    self.old_pairs.append(pairs)
                    self.logBonds(step, cutoff)

                    if len(self.old_pairs) > 0.95 * int(self.maxBonds):
                        self.finishSim(folderName)
                        step = 0
                        break
                else:
                    self.finishSim(folderName) 
                    step = 0
                    break

            self.gro = deepcopy(self.initGro)
            self.top = deepcopy(self.initTop)

    def getMolNames(self):
        names = []
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
    
    def preparePara(self):
        import prepareParam
        
        path = self.srcPath
        basicFolder = '{}/{}'.format(path, 'basic')
        systemsFolder = '{}/{}'.format(path, 'systems')
        mdpFolder = '{}/{}'.format(path, 'mdp')
        self.basicFolder = basicFolder
        self.systemsFolder = systemsFolder
        self.mdpFolder = mdpFolder
        
        self.setParam('options.txt')
        self.getMolNames()
        
        self.unrctFolder = '{}/unrctSystem/'.format(self.systemsFolder)
        self.rctFolder = '{}/rctSystem/'.format(self.systemsFolder)
        self.typeFolder = '{}/typeSystem/'.format(self.systemsFolder)
        
        os.chdir(self.unrctFolder)
        for n in self.molNames: # Need test
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

        os.chdir(self.srcPath)
        
        if os.path.isfile('{}/charges.txt'.format(self.basicFolder)):
            self.chargeMap = self.getChargeMaps('{}/charges.txt'.format(self.basicFolder))
        else:
            a1 = generateChargeDb.generateChargeDb()
            cc = a1.main(self.unrctFolder, self.rctFolder, 4) # could be more
            
            self.chargeMap = cc
        
        os.chdir(self.srcPath)
        a = generateTypeInfo.generateTypeInfo()
        a.main(self.unrctFolder, self.typeFolder)
        
if __name__ == '__main__':
    a = main(4) # change name like gmx_cl ....
    a.preparePara()
    a.mainProcess(2, ig=False)

    
    # TODO: need to check that charge been update as the template. 
