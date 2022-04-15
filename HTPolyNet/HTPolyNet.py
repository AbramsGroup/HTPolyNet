#!/usr/bin/env python
"""
@author: huang, abrams
"""
import logging
import os
import argparse as ap

''' intrapackage imports '''
from HTPolyNet.configuration import Configuration
from HTPolyNet.coordinates import Coordinates
# import HTPolyNet.searchBonds as searchBonds
# import HTPolyNet.genBonds as genBonds
# import HTPolyNet.generateChargeDb as generateChargeDb
# import HTPolyNet.generateTypeInfo as generateTypeInfo

from HTPolyNet.ambertools import GAFFParameterize
from HTPolyNet.topology import Topology

from HTPolyNet.software import Software, Command
from HTPolyNet.countTime import *
from HTPolyNet.projectfilesystem import ProjectFileSystem, Library

from HTPolyNet.gromacs import insert_molecules, grompp_and_mdrun, analyze_sea
from HTPolyNet.molecule import oligomerize, Oligomer

class HTPolyNet:
    ''' Class for a single HTPolyNet runtime session '''
    def __init__(self,cfgfile='',restart=False):
        self.software=Software()
        logging.info(str(self.software))
        if cfgfile=='':
            logging.error('HTPolyNet requires a configuration file.\n')
            raise RuntimeError('HTPolyNet requires a configuration file.')
        self.cfg=Configuration.read(cfgfile)
        ''' Create the initial file system for the project.  If this is 
            a restart, reference the newest project directory in the current
            directory.  If this is not a restart, generate the *next* 
            project directory. '''
        self.pfs=ProjectFileSystem(root=os.getcwd(),verbose=True,reProject=restart)
        logging.info(f'Read configuration from {cfgfile}')
        ''' initialize an empty topology '''
        self.Topology=Topology(system=self.cfg.Title)
        self.local_data_searchpath=[self.pfs.rootPath,self.pfs.projPath]

    def checkout(self,filename,altpath=None):
        if not self.pfs.library.checkout(filename,nowarn=True):
            searchpath=self.local_data_searchpath
            if altpath:
                searchpath.append(altpath)
            for p in searchpath:
                if os.path.exists(os.path.join(p,filename)):
                    Command(f'cp {os.path.join(p,filename)} .').run()
                    return True
            return False
        return True

    def initialize_topology(self,force_parameterize=False,force_capping=False):
        ''' Create a full gromacs topology that includes all directives necessary 
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined. '''
        self.pfs.cd(self.pfs.unrctPath)
        exists=self.pfs.library.exists
        checkin=self.pfs.library.checkin
        if os.path.isfile('init.top'):
            logging.info(f'init.top already exists in {self.pfs.unrctPath} but we will rebuild it anyway!')
        comp=self.cfg.parameters['composition']
        for n,m in self.cfg.monomers.items():
            already_parameterized = all([exists(f'{n}{ex}') for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']])
            if force_parameterize or not already_parameterized:
                logging.info(f'Fetching input monomer structure {n}.mol2')
                if self.checkout(f'{n}.mol2'):
                    logging.info(f'Parameterizing monomer {n}')
                    GAFFParameterize(n,f'{n}-p',force=force_parameterize,parmed_save_inline=False)
                    for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                        checkin(f'{n}{ex}',overwrite=True)
                else:
                    logging.error(f'No {n}.mol2 found.')
                    raise FileNotFoundError(f'No {n}.mol2 found.')
            else:
                logging.info(f'Fetching parameterized monomer {n}')
                self.checkout(f'{n}.mol2')
                for ex in ['-p.mol2','-p.top','-p.itp','-p.gro']:
                    self.checkout(f'{n}{ex}')
            logging.info(f'Reading {n}-p.top, replicating, merging into global topology.')
            t=Topology.read_gro(f'{n}-p.top')
            m.Topology['active']=t
            t.rep_ex(comp[n])
            logging.debug('initialize_topology merging {comp[n]} copies of {n} into global topology')
            self.Topology.merge(t)
            ''' just in case antechamber changed any user-defined atom names... '''
            userMol2=Coordinates.read_mol2(f'{n}.mol2')
            paramMol2=Coordinates.read_mol2(f'{n}-p.mol2')
            m.update_atom_specs(paramMol2,userMol2)
            m.Coords['active']=paramMol2
            if len(m.capping_bonds)>0:
                already_parameterized = all([exists(f'{n}-capped.{ex}') for ex in  ['mol2','top','itp','gro']])
                if force_capping or not already_parameterized:
                    logging.info('Building capped version of monomer {m}')
                    capped=paramMol2.cap(m.capping_bonds)
                    capped.write_mol2(f'{n}-capped-unminimized.mol2')
                    logging.info(f'Parameterizing {n}-capped\n')
                    GAFFParameterize(f'{n}-capped-unminimized',f'{n}-capped',force=force_capping,parmed_save_inline=False)
                    logging.info(f'Minimizing {n}-capped')
                    self.checkout('em-single-molecule.mdp')
                    grompp_and_mdrun(gro=f'{n}-capped',top=f'{n}-capped',
                                    mdp='em-single-molecule',
                                    out=f'{n}-capped',boxSize=[15.,15.,15.])
                    m.Coords['inactive']=Coordinates.read_gro(f'{n}-capped.gro')
                    capped.copy_coords(m.Coords['inactive'])
                    capped.write_mol2(f'{n}-capped.mol2')
                    for ex in ['mol2','top','itp','gro']:
                        checkin(f'{n}-capped.{ex}',overwrite='if_newer')
                else:
                    logging.info(f'Fetching parameterization of capped version of monomer {n}')
                    for ex in ['mol2','top','itp','gro']:
                        self.checkout(f'{n}-capped.{ex}')
                m.Coords['inactive']=Coordinates.read_gro(f'{n}-capped.gro')
                logging.info(f'Reading {n}-capped.top and merging types into global topology')
                m.Topology['inactive']=Topology.read_gro(f'{n}-capped.top')
                self.Topology.merge_types(t)
        logging.info(f'Extended topology has {self.Topology.atomcount()} atoms.')

    def determine_monomer_sea(self,seamdp='nvt-sea',boxSize=[4,4,4],force_sea_calculation=False):
        ''' Determine symmetry-equivalent atoms in every monomer '''
        self.pfs.cd(self.pfs.unrctPath)
        #exists=self.pfs.library.exists
        checkin=self.pfs.library.checkin
        self.checkout(f'{seamdp}.mdp')
        for n,m in self.cfg.monomers.items():
            success=self.checkout(f'{n}-p.sea')
            logging.debug(f'Should be false {not success} or {force_sea_calculation}')
            if not success or force_sea_calculation:
                logging.info(f'Determining SEA for molecule {n}')
                for ex in ['top','itp','gro']:
                    self.checkout(f'{n}-p.{ex}')
                logging.info(f'  hot md running...output to {n}-p-sea')
                grompp_and_mdrun(gro=f'{n}-p',top=f'{n}-p',
                                mdp=seamdp,out=f'{n}-p-sea',boxSize=boxSize)
                m.Coords['active'].D['atoms']['sea-idx']=analyze_sea(f'{n}-p-sea')
                m.Coords['active'].write_sea(f'{n}-p.sea')
                checkin(f'{n}-p.sea')
            if len(m.capping_bonds)>0:
                if not self.checkout(f'{n}-capped.sea') or force_sea_calculation:
                    logging.info(f'Determining SEA for molecule {n}-capped')
                    for ex in ['mol2','top','itp','gro']:
                        self.checkout(f'{n}-capped.{ex}')
                    logging.info(f'  hot md running...output to {n}-capped-sea')
                    grompp_and_mdrun(gro=f'{n}-capped',top=f'{n}-capped',
                                    mdp=seamdp,out=f'{n}-capped-sea',boxSize=boxSize)
                    m.Coords['inactive'].D['atoms']['sea-idx']=analyze_sea(f'{n}-capped-sea')
                    m.Coords['inactive'].write_sea(f'{n}-capped.sea')
                    checkin(f'{n}-capped.sea')

    def make_oligomer_templates(self,force_parameterize=False):
        logging.info(f'Building oligomer templates in {self.pfs.rctPath}')
        exists=self.pfs.library.exists
        checkin=self.pfs.library.checkin
        self.pfs.cd(self.pfs.rctPath)
        self.cfg.oligomers={}
        for r in self.cfg.reactions:
            logging.info(f'Executing reaction(s) for {r.reactants}')
            mname,nname=r.reactants
            m=self.cfg.monomers[mname]
            n=self.cfg.monomers[nname]
            these_oligos=oligomerize(m,n)
            logging.debug(f'Preparing to parameterize {len(these_oligos)} oligomers')
            for o,olig in these_oligos.items():
                logging.debug(f'   olig {o}')
                already_parameterized=all([exists(f'{o}-p.{ex}') for ex in ['mol2','top','itp','gro']])
                if force_parameterize or not already_parameterized:
                    logging.info(f'Parameterizing oligomer template {o}')
                    olig.coord.write_mol2(f'{o}.mol2')
                    GAFFParameterize(f'{o}',f'{o}-p',force=force_parameterize,parmed_save_inline=False)
                    logging.info(f'Minimizing oligomer template {o}')
                    self.checkout('em-single-molecule.mdp')
                    grompp_and_mdrun(gro=f'{o}-p',
                                    top=f'{o}-p',
                                    mdp='em-single-molecule',
                                    out=f'{o}-p',boxSize=[20.,20.,20.])
                    tmp=Coordinates.read_gro(f'{o}-p.gro')
                    olig.coord.copy_coords(tmp)
                    olig.coord.write_mol2(f'{o}-p.mol2')
                    for ex in ['mol2','top','itp','gro']:
                        checkin(f'{o}-p.{ex}',overwrite='if_newer')
                else:
                    logging.info(f'Oligomer {o} is already parameterized.')
                    for ex in ['mol2','top','itp','gro']:
                        self.checkout(f'{o}-p.{ex}')
                logging.info(f'Reading {o}-p.top and merging types into global topology')
                t=Topology.read_gro(f'{o}-p.top')
                olig.topology=t
                olig.unlinked_topology=Topology()
                logging.info(f'Oligo {o} stoich: {olig.stoich}')
                for m,c in olig.stoich:
                    logging.info(f'Providing oligo {o} base topology from {m.name}')
                    self.checkout(f'{m.name}-p.top')
                    self.checkout(f'{m.name}-p.itp')
                    t=Topology.read_gro(f'{m.name}-p.top')
                    t.rep_ex(c)
                    olig.unlinked_topology.merge(t)
                
                #olig.analyze()

                self.cfg.oligomers[o]=olig

                self.Topology.merge_types(olig.topology)

    def write_global_topology(self,filename='init.top',dest=''):
        if dest=='':
            dest=self.pfs.unrctPath
        self.pfs.cd(dest)
        self.Topology.to_file(filename)
        logging.info(f'Wrote {filename} to {dest}')

    def do_liquid_simulation(self):
        # go to the results path, make the directory 'init', cd into it
        self.pfs.cd(self.pfs.next_results_dir())
        # fetch unreacted init.top amd all monomer gro's 
        # from parameterization directory
        self.checkout('init.top',altpath=self.pfs.unrctPath)
        for n in self.cfg.monomers.keys():
            self.checkout(f'{n}-p.gro',altpath=self.pfs.unrctPath)
        # fetch mdp files from library, or die if not found
        self.checkout('em.mdp')
        self.checkout('npt-1.mdp')
        if 'initial_boxsize' in self.cfg.parameters:
            boxsize=self.cfg.parameters['initial_boxsize']
        elif 'initial_density' in self.cfg.parameters:
            mass_kg=self.Topology.total_mass(units='SI')
            V0_m3=mass_kg/self.cfg.parameters['initial_density']
            L0_m=V0_m3**(1./3.)
            L0_nm=L0_m*1.e9
            logging.info(f'Initial density {self.cfg.parameters["initial_density"]} kg/m^3 and total mass {mass_kg} kg dictate an initial box side length of {L0_nm} nm')
            boxsize=[L0_nm,L0_nm,L0_nm]
        # extend system, make gro file
        msg=insert_molecules(self.cfg.monomers,self.cfg.parameters['composition'],boxsize,'init',basename_modifier='-p')
        logging.info(msg)
        self.Coordinates=Coordinates.read_gro('init.gro')
        assert self.Topology.atomcount()==self.Coordinates.atomcount(), 'Error: Atom count mismatch'
        logging.info('Generated init.top and init.gro.')
        msg=grompp_and_mdrun(gro='init',top='init',out='min-1',mdp='em')
        logging.info(msg)
        # TODO: modify this to run in stages until volume is equilibrated
        msg=grompp_and_mdrun(gro='min-1',top='init',out='npt-1',mdp='npt-1')
        logging.info(msg)
        logging.info('Final configuration in npt-1.gro\n')
        sacmol=Coordinates.read_gro('npt-1.gro')
        self.Coordinates.copy_coords(sacmol)
        self.pfs.cdroot()

    def SCUR(self):
        # Search - Connect - Update - Relax
        self.initialize_reactive_topology()
        max_nxlinkbonds=self.D['atoms']['z'].sum()/2
        scur_complete=False
        scur_search_radius=self.cfg.parameters['SCUR_cutoff']
        desired_conversion=self.cfg.parameters['conversion']
        radial_increment=self.cfg.parameters.get('SCUR_radial_increment',0.5)
        maxiter=self.cfg.parameters.get('maxSCURiter',20)
        iter=0
        while not scur_complete:
            logging.info(f'SCUR iteration {iter} begins')
            scur_complete=True
            # TODO: everything -- identify bonds less than radius
            # make bonds, relax
            num_newbonds=self.scur_make_bonds(scur_search_radius)
            curr_nxlinkbonds=max_nxlinkbonds-self.D['atoms']['z'].sum()
            curr_conversion=curr_nxlinkbonds/max_nxlinkbonds
            scur_complete=curr_conversion>desired_conversion
            scur_complete=scur_complete or iter>maxiter
            if not scur_complete:
                if num_newbonds==0:
                    logging.info(f'No new bonds in SCUR iteration {iter}')
                    logging.info(f'-> updating search radius to {scur_search_radius}')
                    scur_search_radius += radial_increment
                    logging.info(f'-> updating search radius to {scur_search_radius}')
            iter+=1
            logging.info(f'SCUR iteration {iter} ends')
            logging.info(f'Current conversion: {curr_conversion}')
            logging.info(f'   SCUR complete? {scur_complete}')
        logging.info(f'SCUR iterations complete.')

    def make_scur_bonds(self,radius):
        adf=self.Coord.D['atoms']
        raset=adf[(adf['rctvty'].isin('HT'))&(adf['z']>0)]
        return 0

    def initialize_reactive_topology(self):
        ''' adds the 'rctvty' and 'z' attributes to each Coord atom 
            rctvty: reactivity flag, 'H' or 'T' or 'N'
            H's may only bond to 'T's, and 'N''s cannot bond.
            z: number of crosslink bonding positions available at 
            this atom
        '''
        adf=self.Coords.D['atoms']
        rctvty=[]
        atz=[]
        ''' first pass -- assign independently using
            monomer templates '''
        for i,r in adf.iterrows():
            molname=r['resName']
            if not molname in self.cfg.monomers:
                logging.error(f'Molecule {molname} is not in the monomer list.')
                raise Exception('bug')
            m=self.cfg.monomers[molname]
            atomname=r['atomName']
            if atomname in m.reactive_atoms:
                ht=m.reactive_atoms[atomname].ht
                z=m.reactive_atoms[atomname].z
            else:
                ht='N'
                z=0
            rctvty.append(ht)
            atz.append(z)
        adf['rctvty']=rctvty
        adf['z']=atz

        ''' second pass -- update individual atoms rctvty and z based
            on current bonding topology '''
        # for i in range(len(adf)):
        #     idx=i+1
        #     iz=adf.iloc[i]['z']
        #     ir=adf.iloc[i]['rctvty']
        #     imolname=adf.iloc[i]['resName']
        #     iatomname=adf.iloc[i]['atomName']
        #     im=self.cfg.monomers[imolname]
        #     imaxz=0
        #     if iatomname in im.reactive_atoms:
        #         imaxz=im.reactive_atoms[iatomname].z
        #     jni=bondlist.partners_of(idx)
        #     for jdx in jni:
        #         j=jdx-1
        #         jr=adf.iloc[j]['rctvty']
        #         if list(sorted([ir,jr]))==['H','T']:
        #             adf.iloc[i]['z']-=1
        #             jmolname=adf.iloc[j]['resName']
        #             jatomname=adf.iloc[j]['atomName']
        #             jm=self.cfg.monomers[jmolname]
        #             jmaxz=0
        #             if jatomname in jm.reactive_atoms:
        #                 jmaxz=jm.reactive_atoms[jatomname].z




    def initreport(self):
        print(self.cfg)
        print()
        print(self.software)

    # def getGroInfo(self, name):
    #     a = readGro.initGro()
    #     a.setName(name)
    #     df_init, sysName, atNum, boxSize = a.readGRO()
    #     m = groInfo.gro()
    #     m.setGroInfo(df_init, sysName, atNum, boxSize)
    #     self.gro = m

    # def getTopInfo(self, topName, itpName):
    #     a = readTop2.initTop()
    #     a.setName(topName, itpName)
    #     a.genTopSession()
    #     b = topInfo.top()
    #     b.setInfo(a.sumTop)
    #     b.checkCharge()
    #     return b

    # @countTime
    # def updateCoord(self, name):
    #     a = readGro.initGro()
    #     a.setName(name)
    #     df_init, sysName, atNum, boxSize = a.readGRO()
    #     a1 = groInfo.gro()
    #     a1.setGroInfo(df_init, sysName, atNum, boxSize)
    #     self.gro.updateCoord(a1)
    #     self.initGro.updateCoord(a1)
        
    # def initSys(self):
    #     print('-> Creating mixture...')
    #     param = self.basicParameter
    #     molInfo = {}
    #     nameList = []
    #     monInfo = param.monInfo
    #     croInfo = param.croInfo
    #     for i in monInfo:
    #         molInfo[i[1]] = i[2]
    #         nameList.append(i[1])
        
    #     for i in croInfo:
    #         molInfo[i[1]] = i[2]
    #         nameList.append(i[1])

    #     os.mkdir('init'); os.chdir('init')
    #     copyfile('{}/npt-1.mdp'.format(self.mdpFolder), 'npt-1.mdp')
    #     copyfile('{}/em.mdp'.format(self.mdpFolder), 'em.mdp')

    #     for n in nameList:
    #         copyfile('{}/{}.gro'.format(self.unrctFolder, n), '{}.gro'.format(n))
    #         copyfile('{}/{}.top'.format(self.unrctFolder, n), '{}.top'.format(n))
    #         copyfile('{}/{}.itp'.format(self.unrctFolder, n), '{}.itp'.format(n))
            
    #     # Insert molecules to systems
    #     import HTPolyNet.extendSys as extendSys
    #     a = extendSys.extendSys('gmx_mpi')
    #     a.extendSys(param.monInfo, param.croInfo, param.boxSize, 'init')
    #     # creates init.gro ONLY
        
    #     # Get df of gro file
    #     self.getGroInfo('init')
    #     # my syntax:
    #     #self.gro=Coordinates.from_groFile('init.gro')
        
    #     # Get parameters from parameters file
    #     topList = []
    #     for n in nameList:
    #         a = self.topMap[n]
    #         # a = self.getTopInfo('{}.top'.format(n), '{}.itp'.format(n))
    #         nNum = int(molInfo[n])
    #         for i in range(nNum):
    #             topList.append(a)
                
    #     # Get sum of top
    #     topSum = mergeTop.mergeTopList(topList)
    #     sysTop = topSum.outDf('init')
    #     self.top = topSum
    #     self.initTop = deepcopy(self.top)

    #     print('-> Successful creating mixture!!')

    #     # EM and NPT to equilibrate the structure
    #     print('-> Conduct energy minization on the new mixture')
    #     a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
    #     a.emSimulation('init', 'init', 'min-1', size=False)
    #     print('-> Conduct NPT on the new mixture')
    #     boxSize = param.boxSize
    #     if boxSize[0] == boxSize[1] == boxSize[2]:
    #         a.NPTSimulation('min-1', 'init', 'npt-1', 'npt-1', check=True, re=False)
    #     else:
    #         copyfile('{}/npt-l.mdp'.format(self.mdpFolder), 'npt-l.mdp')
    #         a.NPTSimulation('min-1', 'init', 'npt-l', 'npt-l', check=True, re=False)
    #         a.NPTSimulation('npt-l', 'init', 'npt-1', 'npt-1', check=True, re=False)
    #     i = 0
    #     # TODO: can not ensure the NPT is finished well
    #     print('-> The mixture is good to go!!')
    #     while(not a.checkMDFinish('npt-1')):
    #         if i > 5:
    #             print('Still cannot converge NPT system, restart')
    #             sys.exit()
    #         elif i == 0:
    #             inName = 'npt-1'
    #         else:
    #             inName = 'npt-1' + str(i)
    #         a.extraRun(inName, 'init', i)
    #         if os.path.isfile('npt.gro'):
    #             move('npt.gro', 'npt-1.gro')
    #         i += 1

    #     # init rct info for potential atoms, add rct columns to the df in the gro object df
    #     self.gro.initRctInfo(self.basicParameter)
    #     self.initGro = deepcopy(self.gro)

    #     # Update coord to the gro df
    #     self.updateCoord('npt-1')

    #     # Back to the working directory and start crosslinking approach
    #     os.chdir(self.resFolder)

    # def countConv(self):
    #     num1 = 0
    #     for i in self.old_pairs:
    #         num1 += len(i)

    #     conv = round(num1 / int(self.cfg.maxBonds), 2)
    #     return conv

    # def stepCapping(self):
    #     os.mkdir('capping')
    #     os.chdir('capping')
    #     gro = deepcopy(self.gro)
    #     top = deepcopy(self.top)
    #     a = endCapping.endCapping(gro, top, self.basicFFType, self.unrctMap, self.cappingBonds)
    #     gro = a.gro
    #     top = a.top
    #     top.endCappingtopClean()

    #     gro.outDf('sys')
    #     top.topClean(key='bonds')
    #     top.outDf('sys')
    #     os.chdir('..')

    # def finishSim(self, folderName, conv, step=0):
    #     os.chdir('..')
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

#         os.chdir(self.resFolder)
#         conv = self.countConv()
#         move(folderName, '{}-{}'.format(folderName, conv))
#         self.old_pairs = []
#         self.reInitSys()
        
#     def reInitSys(self):
#         path = '{}/init'.format(self.resFolder)
#         os.chdir(path)
#         # Get df of gro file
#         gro = readGro.initGro()
#         top = readTop2.initTop()
        
#         gro.setName('init')
#         df_init, sysName, atNum, boxSize = gro.readGRO()
#         atomsDf = groInfo.gro()
#         atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
        
#         self.gro = atomsDf
#         self.gro.initRctInfo(self.basicParameter)
#         if os.path.isfile('npt-init.gro'):
#             self.updateCoord('npt-init')
#         else:
#             self.updateCoord('npt-1')
#         topDf = topInfo.top()
#         top.setName('init.top', 'init.itp')
#         top.genTopSession()
#         topDf.setInfo(top.sumTop)
#         # topDf.checkCharge()
#         self.top = topDf
#         os.chdir(self.resFolder)
    
#     def setupFolder(self, idx):
#         if os.path.isdir('step{}'.format(idx)):
#             rmtree('step{}'.format(idx))
#         os.mkdir('step{}'.format(idx))
#         copyfile('{}/em.mdp'.format(self.mdpFolder), '{}/em.mdp'.format('step{}'.format(idx)))
#         copyfile('{}/npt-cl.mdp'.format(self.mdpFolder), '{}/npt-cl.mdp'.format('step{}'.format(idx)))
#         copyfile('{}/npt-sw.mdp'.format(self.mdpFolder), '{}/npt-sw.mdp'.format('step{}'.format(idx)))
#         return 'step{}'.format(idx)
    


#     def logBonds(self, step, cutoff):
#         num1 = 0
#         for i in self.old_pairs:
#             num1 += len(i)
# #        num1 = len(self.old_pairs)
#         conv = num1/self.maxBonds
#         self.conv = conv

#         with open('../bond.txt', 'a') as f1:
# #            str1 = 'step {} generate {} bonds. {} bonds left. Reach conversion {:.2f}\n'.format(step, 
# #                         num1, self.maxBonds - num1, conv)
#             str1 = 'step {}: {} bonds are formed, within cutoff {}nm. {} bonds left. Reach conversion {:.2f}\n'.format(step,
#                          len(self.old_pairs[int(step)]), round(cutoff, 2), self.maxBonds - num1, conv)
#             f1.write(str1)
        
#         with open('../bonds_connect_Info.txt', 'w') as f2:
#             str0 = 'Total bonds: {}\n'.format(self.maxBonds)
#             # values = self.pairs_detail['step{}'.format(step)]
#             f2.write(str0)
#             for keys, values in self.pairs_detail.items():
#                 f2.write('{}: \n'.format(keys))
#                 for index, row in values.iterrows():
#                     f2.write('atom1: {}\t{}\t atom2: {}\t{}\t mol1: {}\t mol2: {}\n'.format(
#                         row.amon, row.monAtomName, row.acro, row.croAtomName, row.monMol, row.croMol))

#     @countTime
#     def stepwiseRelax(self):
#         print('----> Start stepwise relaxation step to relax the system')
#         k = self.stepwise
#         outName = 'sw'
#         for i in range(len(k)):
#             groName = outName
#             topName = '{}-{}'.format(outName, i)
#             self.gro.outDf(groName)
#             self.top.outDf(topName, float(k[i]), simple=False, stepRelax=True)
#             a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
#             cond0 = a.emSimulation(groName, topName, 'sw-min-{}'.format(i), size=False, check=False)
#             if cond0 == False:
#                 print('----> Stepwised EM failed')
#                 return False

#             cond1 = a.NPTSimulation('sw-min-{}'.format(i), topName,
#                                     'sw-npt-{}'.format(i), 'npt-sw',
#                                     check=False, re=True)
#             if cond1 == False:
#                 print('----> Stepwised NPT failed')
#                 return False

#             print('----> Stepwised step for k = {} is succuessful'.format(k[i]))
#             self.updateCoord('sw-npt-{}'.format(i))

#         return True

#     def mainProcess(self, repeatTimes):
#         # Init systems
#         os.chdir(self.resFolder)
#         self.initSys()  # returns to resFolder
#         conv = 0
#         # calculate max bonds
#         # self.calMaxBonds()

#         # Start crosslinking approach
#         step = 0        
#         for i in range(repeatTimes):
#             print('--> Start crosslinking procedure on replica {}!!'.format(i))
#             folderName = 'sim{}'.format(i)
#             os.mkdir(folderName); os.chdir(folderName)
#             print('---> A new NPT to mix the system at 300K to create new replica')
#             os.mkdir('init'); os.chdir('init')

#             self.top.outDf('init')
#             self.gro.outDf('init')

#             copyfile('{}/npt-init.mdp'.format(self.mdpFolder), 'npt-init.mdp')
#             copyfile('{}/em.mdp'.format(self.mdpFolder), 'em.mdp')
#             a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
#             a.emSimulation('init', 'init', 'min-1', size=False)
#             a.NPTSimulation('min-1', 'init', 'npt-init', 'npt-init', check=False, re=False)
            
#             # complete structs
#             tmpGro = readGro.initGro()
#             tmpGro.setName('npt-init')
#             tmp_df, tmp_name, tmp_atNum, tmp_boxSize = tmpGro.readGRO()
#             cmd1 = f'echo 0 0 | gmx_mpi trjconv -s npt-init.tpr -f npt-init.gro -o npt-init.gro -center -pbc atom -box {tmp_boxSize}'
#             subprocess.call(cmd1, shell=True)
#             self.updateCoord('npt-init')
#             os.chdir('..')
#             print('---> New replica is good to go')
#             boxLimit = self.boxLimit # setup layer curing condition
            
#             while(len(self.old_pairs) < int(self.maxBonds)):
#                 print('---> step {}'.format(step))
#                 folderName1 = self.setupFolder(step)
#                 os.chdir(folderName1)
#                 print('     (Content can be found under folder {})'.format(os.getcwd()))
#                 # searching potential bonds
#                 print('----> Start searching bonds')
#                 if not self.layer_status:
#                     if conv < boxLimit * self.layerConvLimit:
#                         boxLimit = self.boxLimit
#                     else:
#                         print('1st layer conversion reached desired {} conversion'.format(self.layerConvLimit))
#                         boxLimit = 1
#                         self.layer_status = True
#                 sbonds = searchBonds.searchBonds(self.cpu, self.basicParameter, self.old_pairs, self.gro, self.top,
#                                                     self.conv, self.desBonds, self.chains, boxLimit)
#                 pairs, chains, rMols, cutoff = sbonds.sBonds()
#                 # intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']
#                 self.chains = chains
#                 if len(pairs) > 0:
#                     print('----> Start generating bonds')
#                     self.pairs_detail['step{}'.format(step)] = pairs

#                     # generate bonds
#                     gbonds = genBonds.genBonds(self.gro, self.top, pairs, self.chargeMap, rMols, cat='map')
#                     gbonds.gBonds() # update atom's rct status

#                     self.gro = gbonds.gro
#                     self.top = gbonds.top
#                     self.top.checkCharge()

#                     cond = self.stepwiseRelax()
#                     if cond == False:
#                         print('----> Stepwised cannot relax the systems, some wired bonds may formed in the step, will start a new replica')
#                         self.finishSim(folderName, 0, step=step)
#                         step = 0
#                         break

#                     groName = 'cl-{}'.format(i); topName = 'init'
#                     self.gro.outDf(groName)
#                     self.top.outDf(topName)

#                     # intDf = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']

#                     # Equilibrate system
#                     print('----> Energy minimization on the normal system')
#                     a = md.md('gmx_mpi', 'mpirun', self.cpu, nGPU=self.gpu)
#                     cond0 = a.emSimulation(groName, topName, 'min-1', size=False, check=False)
#                     if cond0 == False:
#                         print('EM failed')
#                         self.finishSim(folderName, 0, step=step)
#                         step = 0
#                         break
#                     print('----> NPT on the normal system')
#                     cond1 = a.NPTSimulation('min-1', topName, 'npt-cl', 'npt-cl', check=False, re=True)
#                     if cond1 == False:
#                         print('NPT failed')
#                         self.finishSim(folderName, 0, step=step)
#                         step = 0
#                         break
                    
#                     # Update coord
#                     self.updateCoord('npt-cl')
#                     self.old_pairs.append(pairs)
#                     self.logBonds(step, cutoff)
#                     # self.stepCapping()
#                     conv = self.countConv()
#                     print('----> step {} reaches {} conversion'.format(step, round(conv, 2)))
#                     if conv >= self.desConv:
#                         self.finishSim(folderName, conv, step=step)
#                         step = 0
#                         break

#                     self.prevGro = deepcopy(self.gro)
#                     self.prevTop = deepcopy(self.top)
#                     os.chdir('..')
#                     step += 1
#                 else:
#                     self.finishSim(folderName, conv, step=step)
#                     step = 0
#                     break

#             self.gro = deepcopy(self.initGro)
#             self.top = deepcopy(self.initTop)


    
#     def getChargeMaps(self, name):
#         maps = {}
#         with open(name, 'r') as f:
#             for i in f.readlines():
#                 key, value = i.split(':')
#                 maps[key] = value
#         return maps

#     def getUnrctPara(self):
#         atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
#         btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
#         angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
#         dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
#         impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
#         atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']

#         basicName = self.cfg.unrctStruct
#         aTypes = pd.DataFrame(columns=atypeNames)
#         bTypes = pd.DataFrame(columns=btypeNames)
#         angTypes = pd.DataFrame(columns=angTypeNames)
#         dihTypes = pd.DataFrame(columns=dihTypeNames)
#         impTypes = pd.DataFrame(columns=impTypeNames)
#         atoms = pd.DataFrame(columns=atNames)

#         for name in basicName:
#             top = readTop2.initTop()
#             top.setName('{}.top'.format(name), '{}.itp'.format(name))
#             top.genTopSession()
#             topSum = top.sumTop
#             aTypes = aTypes.append(topSum[1], ignore_index=True)
#             bTypes = bTypes.append(topSum[3], ignore_index=True)
#             angTypes = angTypes.append(topSum[4], ignore_index=True)
#             dihTypes = dihTypes.append(top.fullDihTypes, ignore_index=True)
#             impTypes = impTypes.append(topSum[6], ignore_index=True)
#             topSum[7].residue = name.split('-')[0]
#             atoms = atoms.append(topSum[7], ignore_index=True)

#         aTypes.drop_duplicates(inplace=True, ignore_index=True)
#         bTypes.drop_duplicates(inplace=True, ignore_index=True)
#         angTypes.drop_duplicates(inplace=True, ignore_index=True)
#         dihTypes.drop_duplicates(inplace=True, ignore_index=True)
#         impTypes.drop_duplicates(inplace=True, ignore_index=True)
#         self.basicFFType = [aTypes, bTypes, angTypes, dihTypes, impTypes, atoms]

#         unrctMap = getCappingParam.genUnrctMapping(self.cfg)
#         self.unrctMap = unrctMap

    def main(self,**kwargs):
        force_capping=kwargs.get('force_capping',False)
        force_parameterize=kwargs.get('force_parameterize',False)
        force_sea_calculation=kwargs.get('force_sea_calculation',False)
        self.initialize_topology(force_capping=force_capping,force_parameterize=force_parameterize)
        self.determine_monomer_sea(force_sea_calculation=force_sea_calculation)
        self.make_oligomer_templates(force_parameterize=force_parameterize)
        self.write_global_topology(dest=self.pfs.unrctPath)
        self.do_liquid_simulation()
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
    #     a=generateTypeInfo.generateTypeInfo(pfs,self.cfg)
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
    parser.add_argument('-log',type=str,default='htpolynet_runtime.log',help='log file')
    parser.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    parser.add_argument('--force-capping',default=False,action='store_true',help='force GAFF reparameterization any monomer capping directives in config file')
    parser.add_argument('--force-parameterize',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    parser.add_argument('--force-sea-calculation',default=False,action='store_true',help='force calculation of symmetry-equivalent atoms in any input mol2 structures')
    args=parser.parse_args()

    logging.basicConfig(filename=args.log,encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)
    logging.info('HTPolyNet Runtime begins.')


    if args.command=='info':
        info()
    elif args.command=='run':
        a=HTPolyNet(cfgfile=args.cfg,restart=args.restart)
        a.main(force_capping=args.force_capping,force_parameterize=args.force_parameterize,force_sea_calculation=args.force_sea_calculation)
    else:
        print(f'HTPolyNet command {args.command} not recognized')
    
    logging.info('HTPolynet Runtime ends.')
