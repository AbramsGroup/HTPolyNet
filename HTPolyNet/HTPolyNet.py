#!/usr/bin/env python
"""
@author: huang, abrams
"""

import os
import sys
from shutil import copyfile, copy, move, rmtree
from copy import deepcopy
import subprocess
import argparse as ap
import pandas as pd


''' intrapackage imports '''
from HTPolyNet.configuration import Configuration
from HTPolyNet.coordinates import Coordinates
import HTPolyNet.searchBonds as searchBonds
import HTPolyNet.genBonds as genBonds
import HTPolyNet.generateChargeDb as generateChargeDb
import HTPolyNet.generateTypeInfo as generateTypeInfo

from HTPolyNet.ambertools import GAFFParameterize
from HTPolyNet.topology import Topology

from HTPolyNet.software import Software
from HTPolyNet.countTime import *
from HTPolyNet.projectfilesystem import ProjectFileSystem, Library
from HTPolyNet.gromacs import insert_molecules, grompp_and_mdrun
from HTPolyNet.molecule import react_mol2, Oligomer

class HTPolyNet:
    ''' Class for a single HTPolyNet runtime session '''
    def __init__(self,cfgfile='',logfile='',restart=False):
        self._openlog(logfile)
        self.log('HTPolyNet runtime begins.\n')
        # will die if software requirements are not met
        self.software=Software()
        self.log(str(self.software))
        
        if cfgfile=='':
            raise RuntimeError('HTPolyNet requires a configuration file.')

        self.restart=restart
        
        self.cfg=Configuration.read(cfgfile)
        self.log(str(self.cfg))
        ''' Create the initial file system for the project.  If this is 
            a restart, reference the newest project directory in the current
            directory.  If this is not a restart, generate the *next* 
            project directory. '''
        self.pfs=ProjectFileSystem(root=os.getcwd(),reProject=self.restart)
        self.log('Finished initialization.\n')

    def _openlog(self,filename):
        self.logio=None
        if filename!='':
            self.logio=open(filename,'w')

    def log(self,msg):
        if self.logio:
            self.logio.write(msg)
            self.logio.flush()

    def initialize_topology(self,force_parameterize=False,force_capping=False):
        ''' Create a full gromacs topology that includes all directives necessary 
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined. '''
        self.pfs.cd(self.pfs.unrctPath)
        fetch=self.pfs.fetch
        exist=self.pfs.exist
        store=self.pfs.store
        ''' initialize an empty topology '''
        self.Topology=Topology(system=self.cfg.Title)
        if os.path.isfile('init.top'):
            self.log('init.top already exists in {self.pfs.unrctPath} but we will rebuild it anyway!')
        comp=self.cfg.parameters['composition']
        for n,m in self.cfg.monomers.items():
            if not fetch(f'{n}.mol2'):
                raise FileNotFoundError(f'{n}.mol2 not found in project root, parent, or library')
            already_parameterized = all([exist(f'{n}{ex}') for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']])
            msg=''
            if force_parameterize or not already_parameterized:
                fetch(f'{n}.mol2')
                self.log(f'Parameterizing {n}...\n')
                msg+=GAFFParameterize(n,f'{n}-p',force=False,parmed_save_inline=False)
                for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                    store(f'{n}{ex}',overwrite='if_newer')
            else:
                for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                    fetch(f'{n}{ex}')
            self.log(msg+'\n'+f'Reading {n}-p.top...\n')
            t=Topology.read_gro(f'{n}-p.top')
            m.Topology['active']=t
            self.log(f'...replicating\n')
            t.rep_ex(comp[n])
            self.Topology.merge(t)
            ''' just in case antechamber changed any user-defined atom names... '''
            userMol2=Coordinates.read_mol2(f'{n}.mol2')
            paramMol2=Coordinates.read_mol2(f'{n}-p.mol2')
            m.update_atom_specs(paramMol2,userMol2)
            m.Coords['active']=paramMol2
            if len(m.capping_bonds)>0:
                msg=''
                already_parameterized = all([exist(f'{n}-capped{ex}') for ex in  ['.mol2','.top','.itp','.gro']])
                if force_capping or not already_parameterized:
                    capped=paramMol2.cap(m.capping_bonds)
                    capped.write_mol2(f'{n}-capped-unminimized.mol2')
                    self.log(f'Parameterizing {n}-capped...\n')
                    msg=GAFFParameterize(f'{n}-capped-unminimized',f'{n}-capped',force=False,parmed_save_inline=False)
                    fetch('em-single-molecule.mdp')
                    grompp_and_mdrun(gro=f'{n}-capped',
                                    top=f'{n}-capped',
                                    mdp='em-single-molecule',
                                    out=f'{n}-capped',boxSize=[15.,15.,15.])
                    m.Coords['inactive']=Coordinates.read_gro(f'{n}-capped.gro')
                    capped.copy_coords(m.Coords['inactive'])
                    capped.write_mol2(f'{n}-capped.mol2')
                    for ex in ['-capped.mol2','-capped.top','-capped.itp','-capped.gro']:
                        store(f'{n}{ex}',overwrite='if_newer')
                else:
                    for ex in ['-capped.mol2','-capped.top','-capped.itp','-capped.gro']:
                        fetch(f'{n}{ex}')
                m.Coords['inactive']=Coordinates.read_gro(f'{n}-capped.gro')
                self.log(msg+'\n'+f'Reading {n}-capped.top...\n')
                m.Topology['inactive']=Topology.read_gro(f'{n}-capped.top')
                self.Topology.merge_types(t)
        
        self.log(f'Extended topology has {self.Topology.atomcount()} atoms.\n')
        self.log(f'Extended topology has {len(self.Topology.D["dihedraltypes"])} dihedraltypes.\n')
        assert 'defaults' in self.Topology.D, 'Error: lost defaults?'
        self.cfg.oligomers={}
        for r in self.cfg.reactions:
            self.log(f'Executing reaction(s) for {r.reactants}')
            mname,nname=r.reactants
            m=self.cfg.monomers[mname]
            n=self.cfg.monomers[nname]
            these_oligos=react_mol2(m,n)
            for o,mol in these_oligos.items():
                already_parameterized=all([exist(f'{o}{ex}') for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']])
                if force_parameterize or not already_parameterized:
                    self.log(f'Parameterizing {o}...')
                    mol.write_mol2(f'{o}.mol2')
                    msg=GAFFParameterize(f'{o}',f'{o}-p',force=False,parmed_save_inline=False)
                    fetch('em-single-molecule.mdp')
                    grompp_and_mdrun(gro=f'{o}-p',
                                    top=f'{o}-p',
                                    mdp='em-single-molecule',
                                    out=f'{o}-p',boxSize=[20.,20.,20.])
                    tmp=Coordinates.read_gro(f'{o}-p.gro')
                    mol.copy_coords(tmp)
                    mol.write_mol2(f'{o}-p.mol2')
                    for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                        store(f'{o}{ex}',overwrite='if_newer')
                else:
                    for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                        fetch(f'{o}{ex}')
                self.log(msg+'\n'+f'Reading {o}-parameterized.top...\n')
                self.cfg.oligomers[o]=Oligomer(Topology.read_gro(f'{o}-p.top'),Coordinates.read_gro(f'{o}-p.gro'))
                self.Topology.merge_types(t)

        # write the system topology
        self.Topology.to_file('init.top')
        self.log('Wrote init.top')

    def generate_liquid_simulation(self):
        # go to the results path, make the directory 'init', cd into it
        self.pfs.cd(self.pfs.next_results_dir())
        # fetch unreacted init.top amd all monomer gro's 
        # from parameterization directory
        self.pfs.fetch('init.top',altpath=self.pfs.unrctPath)
        for n in self.cfg.monomers.keys():
            self.pfs.fetch(f'{n}-p.gro',altpath=self.pfs.unrctPath)
        # fetch mdp files from library, or die if not found
        self.pfs.fetch('em.mdp')
        self.pfs.fetch('npt-1.mdp')
        # extend system, make gro file
        msg=insert_molecules(self.cfg.monomers,self.cfg.parameters['composition'],self.cfg.parameters['initial_boxsize'],'init',basename_modifier='-p')
        self.log(msg)
        self.Coordinates=Coordinates.read_gro('init.gro')
        assert self.Topology.atomcount()==self.Coordinates.atomcount(), 'Error: Atom count mismatch'
        self.log('Generated init.top and init.gro.\n')
        msg=grompp_and_mdrun(gro='init',top='init',out='min-1',mdp='em')
        self.log(msg)
        msg=grompp_and_mdrun(gro='min-1',top='init',out='npt-1',mdp='npt-1')
        self.log(msg)
        self.log('Final configuration in npt-1.gro\n')
        sacmol=Coordinates.read_gro('npt-1.gro')
        self.Coordinates.copy_coords(sacmol)
        self.pfs.cdroot()

    def make_oligomer_templates(self):

        # TODO: Typing and charging calculations for all oligomer templates
        #       Types added to topology
        #       Charges added to charge database
        #       1. create all oligomer templates -> mol2
        #       2. antechamber to get charges
        #       3. parmchk2 and tleap to get types

        # TODO: figure out what this is
        # generate the charge map
        # charges_file=f'{self.pfs.basicPath}/charges.txt'
        # if os.path.isfile(charges_file):
        #     self.log(f'--> Reading charge database from {charges_file}\n')
        #     self.chargeMap=self.getChargeMaps(charges_file)
        # else:
        #     self.log('--> Generating new charge database\n')
        #     a1=generateChargeDb.generateChargeDb()
        #     self.chargeMap=a1.main(self.pfs.unrctPath,self.pfs.rctPath,4) # could be more

        # TODO: figure out what this is
        self.pfs.goToProjectRoot()
        self.log('--> Generating reacted molecules type database\n')
        a=generateTypeInfo.generateTypeInfo(self.pfs,self.cfg)
        # TODO: catch the return of generateTypeInfo.main which is a database of all type information
        a.main(self.pfs.unrctFolder,self.pfs.typeFolder)


#       self.Types=getalltypes()
#       self.Topology=getalltopologies()

        # do we really need these attributes if they are all under 'cfg' anyway?
        # self.cpu = ''
        # self.gpu = ''
        # self.trials = ''
        # self.reProject = ''
        # self.stepwise = ''

        # end capping
#        self.cappingBonds = [] # potential bonds for capping
        # self.unrctMap = {}

        # self.topMap = {}
        # self.initGro = ''
        # self.initTop = ''
        # self.dumpPairs = {} # pair map when check circuit

        # cfg
        # self.basicParameter = ''

        # self.basicFFType = []
        # self.molNames = []
        # self.chargeMap = {}

        # self.prevGro = ''
        # self.prevTop = ''

        # self.gro = '' # save gro object
        # self.top = '' # save top object

        # self.chains = []
        # self.old_pairs = []
        # self.pairs_detail = {}

        # needed parameters
        # self.conv = 0
        # self.desConv = 0
        # self.desBonds = 0
        # self.layer_status = False # status whether layer reached desired conversion

    def initreport(self):
        print('Libraries:')
        print(self.Library)
        print()
        print(self.cfg)
        print()
        print(self.software)

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
        import HTPolyNet.extendSys as extendSys
        a = extendSys.extendSys('gmx_mpi')
        a.extendSys(param.monInfo, param.croInfo, param.boxSize, 'init')
        # creates init.gro ONLY
        
        # Get df of gro file
        self.getGroInfo('init')
        # my syntax:
        #self.gro=Coordinates.from_groFile('init.gro')
        
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

        conv = round(num1 / int(self.cfg.maxBonds), 2)
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
        self.initSys()  # returns to resFolder
        conv = 0
        # calculate max bonds
        # self.calMaxBonds()

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
                if not self.layer_status:
                    if conv < boxLimit * self.layerConvLimit:
                        boxLimit = self.boxLimit
                    else:
                        print('1st layer conversion reached desired {} conversion'.format(self.layerConvLimit))
                        boxLimit = 1
                        self.layer_status = True
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


    
    def getChargeMaps(self, name):
        maps = {}
        with open(name, 'r') as f:
            for i in f.readlines():
                key, value = i.split(':')
                maps[key] = value
        return maps

    def getUnrctPara(self):
        atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']

        basicName = self.cfg.unrctStruct
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

        unrctMap = getCappingParam.genUnrctMapping(self.cfg)
        self.unrctMap = unrctMap

    def main(self,**kwargs):
        force_capping=kwargs.get('force_capping',False)
        force_parameterize=kwargs.get('force_parameterize',False)
        self.initialize_topology(force_capping=force_capping,force_parameterize=force_parameterize)
        self.generate_liquid_simulation()
#        self.SCUR()

    # create self.Types and self.Topology, each is a dictionary of dataframes
    # def preparePara(self,log=True):
    #     os.chdir(self.pfs.unrctPath)
    #     if log:
    #         logf=open('parameterization.log','w')
    #         logf.write('Beginning parameterizations.\n')
    #         logf.write('self.cfg.molNames: '+','.join(self.cfg.molNames)+'\n')
    #     #A=ambertools.Parameterization()    print(args)
    #             if log:
    #                 logf.write(msg)
    #             T=processTop.Topology(n) # process the top file, make it standard
    #         T.generate()
    #         self.Topologies.append(T)
    #         self.topMap[T.name] = T.top

    #     self.getUnrctPara()

    #     os.chdir(self.pfs.projPath)
    #     charges_file=f'{self.pfs.basicPath}/charges.txt'
    #     if os.path.isfile(charges_file):
    #         self.chargeMap=self.getChargeMaps(charges_file)
    #     else:
    #         print('--> Generating new charge database')
    #         a1=generateChargeDb.generateChargeDb()
    #         self.chargeMap=a1.main(self.pfs.unrctPath,self.pfs.rctPath,4) # could be more
        
    #     os.chdir(self.pfs.projPath)
    #     print('--> Generating reacted molecules type database')
    #     a=generateTypeInfo.generateTypeInfo(self.pfs,self.cfg)
    #     # TODO: catch the return of generateTypeInfo.main which is a database of all type information
    #     a.main(self.pfs.unrctFolder,self.pfs.typeFolder)
    #     if log:
    #         logf.close()

# def init():
#     LibraryResourcePaths=IdentifyLibraryResourcePaths()
#     example_cfg='VEA-VEB-STY-example.cfg'
#     getme=os.path.join(LibraryResourcePaths["cfg"],example_cfg)
#     print(f'HTPolyNet is copying {example_cfg} from {LibraryResourcePaths["cfg"]}')
#     os.system(f'cp {getme} .')
#     print(f'After editing this file, you can launch using\n"htpolynet run -cfg <name-of-config-file>"')

# def run(a,cfg=''):
#     print(f'HTPolyNet is going to try to run in {os.getcwd()}...')
#     print(dir(a))
#     a.cfg=readCfg.configuration.readCfgFile(cfg)
#     a.parseCfg()
#     for k in dir(a.cfg):
#         if '__' not in k and k in a.cfg.__dict__:
#             print(f'{k} = {a.cfg.__dict__[k]}')
#     print(a.cfg.__dict__)
#     #print(a.cfg.baseDict,a.cfg.rctInfo)
#     # a.preparePara()
#     # a.mainProcess(a.trials)
#     pass

def info():
    print('This is some information on your installed version of HTPolyNet')
    l=Library()
    l.info()
    sw=Software()
    sw.info()

def cli():
    parser=ap.ArgumentParser()
    parser.add_argument('command',type=str,default=None,help='command (init, info, run)')
    parser.add_argument('-cfg',type=str,default='',help='input config file')
    parser.add_argument('-log',type=str,default='out.log',help='log file')
    parser.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    parser.add_argument('--force-capping',default=False,action='store_true',help='force GAFF reparameterization any monomer capping directives in config file')
    parser.add_argument('--force-parameterize',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    args=parser.parse_args()

    if args.command=='info':
        info()
    elif args.command=='run':
        a=HTPolyNet(cfgfile=args.cfg,logfile=args.log,restart=args.restart)
        a.main(force_capping=args.force_capping,force_parameterize=args.force_parameterize)
    else:
        print(f'HTPolyNet command {args.command} not recognized')
