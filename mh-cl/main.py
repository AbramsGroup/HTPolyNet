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
import readTop
import readGro
import groInfo
import topInfo
import md
import searchBonds
import genBonds

import os
from shutil import copyfile
from shutil import move
from shutil import rmtree


class main(object):
    def __init__(self):
        self.srcPath = os.getcwd()
        self.basicFolder = ''
        self.systemsFolder = ''
        self.mdpFolder = ''
        self.workingFolder = '' # current loop folder location
        
        self.basicParameter = ''
        self.molNames = []
        self.chargeMap = {}
        self.gro = '' # save gro object
        self.top = '' # save top object
        self.old_pairs = []
        
        # needed parameters
        self.cutoff = 0.2
        self.maxBonds = 100
        
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
        a = readTop.initTop()
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
        
    def initSys(self):
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
        os.mkdir('init'); os.chdir('init')
        copyfile('{}/npt-init.mdp'.format(self.mdpFolder), 'npt-init.mdp')
        copyfile('{}/nvt.mdp'.format(self.mdpFolder), 'nvt-1.mdp')
        copyfile('{}/em.mdp'.format(self.mdpFolder), 'em.mdp')

        for n in nameList:
            copyfile('{}/{}.gro'.format(self.systemsFolder, n), '{}.gro'.format(n))
            copyfile('{}/{}.top'.format(self.systemsFolder, n), '{}.top'.format(n))
            copyfile('{}/{}.itp'.format(self.systemsFolder, n), '{}.itp'.format(n))
            
        # Insert molecules to systems
        import extendSys
        a = extendSys.extendSys('gmx_mpi')
        a.extendSys(param.monInfo, param.croInfo, param.boxSize, 'init')
        
        # Get df of gro file
        self.getGroInfo('init')
        
        # Get parameters from parameters file
        topList = []
        for n in nameList:
            a = self.getTopInfo('{}.top'.format(n), '{}.itp'.format(n))
            nNum = int(molInfo[n])
            for i in range(nNum):
                topList.append(a)
                
        # Get sum of top
        topSum = mergeTop.mergeTopList(topList)
        sysTop = topSum.outDf('init')
        self.top = topSum
        
        # EM and NPT to equilibrate the structure
        a = md.md('gmx_mpi', 'mpirun', '16')
        a.emSimulation('init', 'init', 'min-1', size=False)
        a.NPTSimulation('min-1', 'init', 'npt-init', 'npt-init', check=False, re=True)
        a.NVTSimulation('npt-init', 'init', 'nvt-1', 'nvt-1', check=False)
        
        # Update coord to the gro df
        self.updateCoord('nvt-1')
        
        # init rct info for potential atoms, add rct columns to the df in the gro object df
        self.gro.initRctInfo(self.basicParameter)
        
        # Back to the working directory and start crosslinking approach
        os.chdir(self.workingFolder)
    
    def finishSim(self, folderName):
        os.chdir(self.workingFolder)
        conv = round(len(self.old_pairs)/int(self.basicParameter.maxBonds), 2)
        move(folderName, '{}-{}'.format(folderName, conv))
        self.old_pairs = []
        self.reInitSys()
        
    def reInitSys(self):
        path = '{}/init'.format(self.workingFolder)
        os.chdir(path)
        # Get df of gro file
        gro = readGro.initGro()
        top = readTop.initTop()
        
        gro.setName('init')
        df_init, sysName, atNum, boxSize = gro.readGRO()
        atomsDf = groInfo.gro()
        atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        
        self.gro = atomsDf
        self.gro.initRctInfo(self.basicParameter)
        self.updateCoord('nvt-1')
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
        return 'step{}'.format(idx)
    
    def mainProcess(self, repeatTimes):
        if os.path.isdir('results'):
            rmtree('results')
        else:
            pass
        os.mkdir('results'); os.chdir('results')  
        
        # Init systems
        self.initSys()
        
        # Start crosslinking approach
        step = 0        
        for i in range(repeatTimes):
            folderName = 'sim{}'.format(i)
            os.mkdir(folderName)
            os.chdir(folderName)

            while(len(self.old_pairs) < int(self.basicParameter.maxBonds)):
                
                # searching potential bonds
                sbonds, rMols = searchBonds.searchBonds(self.basicParameter, self.old_pairs, self.gro, self.top)
                pairs = sbonds.main()
                
                if len(pairs) > 0:
                    self.old_pairs.append(pairs)
                    folderName1 = self.setupFolder(step)  
                    os.chdir(folderName1)
                    # generate bonds
                    gbonds = genBonds.genBonds(self.gro, self.top, pairs, self.chargeMap)
                    gbonds.main()
                    self.gro = gbonds.gro
                    self.top = gbonds.top
                    groName = 'cl-{}'.format(i); topName = 'init'
                    self.gro.outDf(groName)
                    self.top.outDf(topName)
                    
                    # Equilibrate system
                    a = md.md('gmx_mpi', 'mpirun', '16')
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
                    if len(self.old_pairs) > 0.95 * int(self.basicParameter.maxBonds):
                        self.finishSim(folderName)
                        step = 0
                        break
                else:
                    self.finishSim(folderName) 
                    step = 0
                    break
    
    def getMolNames(self):
        names = []
        for n in self.basicParameter.monInfo():
            names.append(n[1])
        
        for n in self.basicParameter.croInfo():
            names.append(n[1])
        
        self.molNames = names
        
    def preparePara(self):
        import prepareParam
        import molRctInfo
        
        path = self.srcPath
        basicFolder = '{}/{}'.format(path, 'basic')
        systemsFolder = '{}/{}'.format(path, 'systems')
        mdpFolder = '{}/{}'.format(path, 'mdp')
        self.basicFolder = basicFolder
        self.systemsFolder = systemsFolder
        self.mdpFolder = mdpFolder
        
        self.setParam('options.txt')
        self.getMolNames()
        
        os.chdir('{}/unrctSystem'.format(self.systemsFolder))
        for n in self.molNames:
            fileName = '{}.mol2'.format(n)
            a = prepareParam.prepareParam()
            a.PrepareFile(fileName, n, n)
        
        os.chdir('{}/rctSystem'.format(self.systemsFolder))
        for n in self.molNames:
            fileName = '{}.mol2'.format(n)
            a = prepareParam.prepareParam()
            a.PrepareFile(fileName, n, n)
        
        a1 = molRctInfo.molRctInfo()
        a2 = a1.main(self.molNames)
        self.chargeMap = a2
        os.chdir(self.srcPath)
        
if __name__ == '__main__':
    a = main()
    a.preparePara()
    a.mainProcess(2)
    
    