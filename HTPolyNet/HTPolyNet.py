#!/usr/bin/env python
"""
@author: huang, abrams
"""
import logging
import os
import shutil
import pandas as pd
import argparse as ap
import numpy as np
from copy import deepcopy
from itertools import product
from multiprocessing import Pool
from functools import partial

''' intrapackage imports '''
from HTPolyNet.configuration import Configuration
from HTPolyNet.topocoord import TopoCoord, BTRC
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.gromacs import insert_molecules, grompp_and_mdrun, density_trace
from HTPolyNet.checkpoint import CPstate, Checkpoint
from HTPolyNet.plot import trace

class HTPolyNet:
    ''' Class for a single HTPolyNet runtime session '''
    def __init__(self,cfgfile='',restart=False):
        logging.info(software.to_string())
        if cfgfile=='':
            logging.error('HTPolyNet requires a configuration file.\n')
            raise RuntimeError('HTPolyNet requires a configuration file.')
        logging.info(f'Configuration: {cfgfile}')
        self.cfg=Configuration.read(os.path.join(pfs.root(),cfgfile))
        self.TopoCoord=TopoCoord()
        self.cfg.parameters['restart']=restart
        if self.cfg.parameters['restart']:
            logging.info(f'***** THIS IS A RESTART *****')

    def checkout(self,filename,altpath=None):
        if not pfs.checkout(filename):
            searchpath=pfs.local_data_searchpath()
            logging.info(f'No {filename} found in libraries; checking local data searchpath {searchpath}')
            if altpath:
                searchpath.append(altpath)
                logging.info(f'and alternative path {altpath}')
            for p in searchpath:
                fullfilename=os.path.join(p,filename)
                if os.path.exists(fullfilename):
                    basefilename=os.path.basename(filename)
                    shutil.copyfile(fullfilename,basefilename)
                    logging.info(f'Found at {fullfilename}')
                    return True
            return False
        return True

    def generate_molecules(self,force_parameterization=False,force_sea_calculation=False,force_checkin=False):
        logging.info('*'*10+' GENERATING MOLECULE TEMPLATES '+'*'*10)
        self.molecules={}
        for mname,M in self.cfg.molecules.items():
            M.set_origin('unparameterized')
        cwd=pfs.go_to('molecules')
        # checkin=pfs.checkin
        # exists=pfs.exists
        ''' Each molecule implied by the cfg is 'generated' here, either by
            reading from the library or direct parameterization.  In some cases,
            the molecule is to be generated by a reaction; if so, it's 
            `generator` attribute will be a Reaction instance '''
        msg=', '.join([m for m in self.cfg.molecules])
        logging.debug(f'Will generate {msg}')
        all_made=all([(m in self.molecules and M.get_origin()!='unparameterized') for m,M in self.cfg.molecules.items()])
        iter=1
        while not all_made:
            logging.debug(f'Pass number {iter} through molecules in cfg')
            iter+=1
            ''' We have to generate molecules that act as precursors to other molecules before
                we create the molecules that need those precursors '''
            for mname,M in self.cfg.molecules.items():
                ''' a molecule with no generator specified is a monomer and can be generated by
                    parameterizing an existing mol2 file or reading in a previous parameterization
                    from the library '''
                if mname not in self.molecules:
                    self.generate_molecule(M,force_parameterization=force_parameterization,force_sea_calculation=force_sea_calculation,force_checkin=force_checkin)
                    if M.get_origin()=='unparameterized':
                        # could not make since waiting on precursors
                        continue
                    self.molecules[mname]=M
                    logging.info(f'Generated {mname}')
            all_made=all([(m in self.molecules and M.get_origin()!='unparameterized') for m,M in self.cfg.molecules.items()])
            logging.info(f'Done making molecules: {all_made}')
            for m in self.cfg.molecules:
                logging.info(f'Mol {m} made? {m in self.molecules} -- origin: {M.get_origin()}')
            
        ''' We need to copy all symmetry info down to atoms in each molecule based on the reactants
            used to generate them.  We then must use this symmetry information to expand the list
            of Reactions '''
        new_molecules=self.cfg.symmetry_expand_reactions(self.molecules)
        for mname,M in new_molecules.items():
            ''' a molecule with no generator specified is a monomer and can be generated by
                parameterizing an existing mol2 file or reading in a previous parameterization
                from the library '''
            if mname not in self.molecules:
                logging.debug(f'Generating {mname}')
                self.generate_molecule(M,force_parameterization=force_parameterization,force_sea_calculation=force_sea_calculation,force_checkin=force_checkin)
                assert M.get_origin()!='unparameterized'
                self.molecules[mname]=M
                logging.info(f'Generated {mname}')
        self.molecules.update(new_molecules)
        self.cfg.maxconv=self.cfg.calculate_maximum_conversion()
        logging.info(f'Maximum conversion is {self.cfg.maxconv} bonds.')
        for M in self.molecules.values():
            logging.debug(f'Ring detector for {M.name}')
            M.label_ring_atoms()
        logging.debug(f'Reaction bond(s) in each molecular template:')
        for M in self.molecules.values():
            if len(M.reaction_bonds)>0:
                logging.debug(f'Template {M.name}:')
                for b in M.reaction_bonds:
                    (i,j),(ri,rj),(A,B)=b
                    logging.debug(f'   {i}({ri}:{A})---{j}({rj}:{B})')

    def generate_molecule(self,M,**kwargs):
        mname=M.name
        checkin=pfs.checkin
        exists=pfs.exists
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        force_sea_calculation=kwargs.get('force_sea_calculation',False)
        if force_parameterization or not M.previously_parameterized():
            logging.debug(f'Parameterization of {mname} requested -- can we generate {mname}?')
            generatable=(not M.generator) or (all([m in self.molecules for m in M.generator.reactants.values()]))
            if generatable:
                logging.info(f'Yes -- calling {mname}.generate()')
                M.generate(available_molecules=self.molecules,**self.cfg.parameters)
                for ex in ['mol2','top','itp','gro']:
                    checkin(f'molecules/parameterized/{mname}.{ex}',overwrite=force_checkin)
                M.set_origin('newly parameterized')
            else:
                logging.debug(f'...no, did not generate {mname}.')
                logging.debug(f'not ({mname}.generator) {bool(not M.generator)}')
                if M.generator:
                    logging.debug(f'reactants {list(M.generator.reactants.values())}')
                return
        else:
            logging.info(f'Fetching parameterized {mname}')
            for ex in ['mol2','top','itp','gro']:
                self.checkout(f'molecules/parameterized/{mname}.{ex}')
            M.load_top_gro(f'{mname}.top',f'{mname}.gro')
            M.set_sequence()
            M.set_reaction_bonds(self.molecules)
            M.TopoCoord.set_gro_attribute('reactantName',M.name)
            M.set_origin('previously parameterized')
        ''' The cfg allows user to indicate whether or not to determine and use
            symmetry-equivalent atoms in any molecule. '''
        if mname in self.cfg.use_sea:
            if force_sea_calculation or not exists(f'molecules/parameterized/{mname}.sea'):
                logging.info(f'Doing SEA calculation on {mname}')
                M.calculate_sea()
                M.analyze_sea_topology()
                M.write_gro_attributes(['sea-idx'],f'{M.name}.sea')
                checkin(f'molecules/parameterized/{mname}.sea',overwrite=force_checkin)
            else:
                logging.debug(f'Reading sea data into {M.name}')
                self.checkout(f'molecules/parameterized/{mname}.sea')
                M.read_gro_attributes(f'{mname}.sea',attribute_list=['sea-idx'])
                M.analyze_sea_topology()
        else:
            M.inherit_attribute_from_reactants('sea-idx',available_molecules=self.molecules)
            M.write_gro_attributes(['sea-idx'],f'{M.name}.sea')
        return True

    def initialize_topology(self,inpfnm='init'):
        """Create a full gromacs topology that includes all directives necessary 
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined.

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        """
        cwd=pfs.go_to('systems/init')
        if os.path.isfile(f'{inpfnm}.top'):
            logging.info(f'{inpfnm}.top already exists in {cwd} but we will rebuild it anyway!')
        ''' for each monomer named in the cfg, either parameterize it or fetch its parameterization '''
        already_merged=[]
        for item in self.cfg.initial_composition:
            M=self.molecules[item['molecule']]
            N=item['count']
            t=deepcopy(M.TopoCoord.Topology)
            t.adjust_charges(0)
            t.rep_ex(N)
            logging.info(f'Merging {N} copies of {M.name}\'s topology into global topology')
            self.TopoCoord.Topology.merge(t)
            already_merged.append(M.name)
        for othermol,M in self.molecules.items():
            if not othermol in already_merged:
                self.TopoCoord.Topology.merge_types(M.TopoCoord.Topology)
        logging.info(f'Global topology has {self.TopoCoord.Topology.atomcount()} atoms.')
        self.TopoCoord.write_top(f'{inpfnm}.top')
        logging.info(f'Wrote {inpfnm}.top to {cwd}')

    def setup_liquid_simulation(self,inpfnm='init'):
        """Builds initial top and gro files for initial liquid simulation

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        """
        cwd=pfs.go_to('systems/init')
        if 'initial_boxsize' in self.cfg.parameters:
            boxsize=self.cfg.parameters['initial_boxsize']
        elif 'initial_density' in self.cfg.parameters:
            mass_kg=self.TopoCoord.total_mass(units='SI')
            V0_m3=mass_kg/self.cfg.parameters['initial_density']
            L0_m=V0_m3**(1./3.)
            L0_nm=L0_m*1.e9
            logging.info(f'Initial density {self.cfg.parameters["initial_density"]} kg/m^3 and total mass {mass_kg:.3e} kg dictate an initial box side length of {L0_nm:.3f} nm')
            boxsize=[L0_nm,L0_nm,L0_nm]
        # extend system, make gro file
        clist=self.cfg.initial_composition
        c_togromacs={}
        for cc in clist:
            c_togromacs[cc['molecule']]=cc['count']
        m_togromacs={}
        for mname,M in self.cfg.molecules.items():
            if mname in c_togromacs:
                m_togromacs[mname]=M
                self.checkout(f'molecules/parameterized/{mname}.gro',altpath=pfs.subpath('molecules'))
        # for m,M in m_togromacs.items():
        #     logging.info(f'Molecule to gromacs: {m} ({M.name})')
        # for m,c in c_togromacs.items():
        #     logging.info(f'Composition to gromacs: {m} {c}')
        if not os.path.exists(f'{inpfnm}.gro') or not os.path.exists(f'{inpfnm}.top'): 
            msg=insert_molecules(c_togromacs,boxsize,inpfnm,**self.cfg.parameters)
            logging.info(f'Generated {inpfnm}.top and {inpfnm}.gro.')
        else:
            logging.info(f'Found {inpfnm}.gro.')
        self.TopoCoord.read_gro(f'{inpfnm}.gro')
        self.TopoCoord.atom_count()
        self.TopoCoord.inherit_attributes_from_molecules(['cycle-idx','reactantName'],self.cfg.molecules)
        self.TopoCoord.set_z(self.cfg.reactions,self.molecules)
        self.TopoCoord.write_gro_attributes(['z','cycle-idx','reactantName'],f'{inpfnm}.grx')
        self.TopoCoord.make_ringlist()
        self.TopoCoord.make_resid_graph()

    def do_liquid_simulation(self,inpfnm='init',deffnm='npt-1'):
        """do_liquid_simulation Manages execution of gmx mdrun to perform minimization
            and NPT MD simulation of the initial liquid system.  Final coordinates are 
            loaded into the global TopoCoord.

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        :param deffnm: deffnm prefix fed to gmx mdrun, defaults to 'npt-1'
        :type deffnm: str, optional
        """
        cwd=pfs.go_to('systems/init')
        if os.path.exists(f'{deffnm}.gro') and self.cfg.parameters['restart']:
            logging.info(f'{deffnm}.gro exists in {os.getcwd()}; skipping initial NPT md.')
        else:
            self.checkout('mdp/em.mdp')
            self.checkout('mdp/npt-1.mdp')
            logging.info(f'Conducting initial NPT MD simulation of liquid')
            msg=grompp_and_mdrun(gro=inpfnm,top=inpfnm,out=f'{inpfnm}-min-1',mdp='em',**self.cfg.parameters)
            msg=grompp_and_mdrun(gro=f'{inpfnm}-min-1',top=inpfnm,out=deffnm,mdp=deffnm,**self.cfg.parameters)
            logging.info(f'Generated configuration {deffnm}.gro\n')
        density_trace(deffnm,**self.cfg.parameters)
        trace(['Density'],[deffnm],outfile='../../plots/init-density.png')
        # update coordinates; will wrap upon reading in
        self.TopoCoord.copy_coords(TopoCoord(grofilename=f'{deffnm}.gro'))

    def CURE(self):
        # Connect - Update - Relax - Equilibrate
        logging.info('*'*10+' CONNECT - UPDATE - RELAX - EQUILIBRATE (CURE) begins '+'*'*10)
        max_nxlinkbonds=self.cfg.maxconv
        desired_conversion=self.cfg.parameters['desired_conversion']
        cure_search_radius=self.cfg.parameters['CURE_initial_search_radius']
        checkpoint_file=self.cfg.parameters.get('checkpoint_file','checkpoint.yaml')
        bonds_file=self.cfg.parameters.get('bonds_file','bonds.csv')
        radial_increment=self.cfg.parameters.get('CURE_radial_increment',0.5)
        maxiter=self.cfg.parameters.get('max_CURE_iterations',20)
        n_stages=self.cfg.parameters.get('max_bond_relaxation_stages',6)
        
        dragging=False
        drag_limit_nm=self.cfg.parameters.get('drag_limit',0.3)
        n_dragstages=self.cfg.parameters.get('max_drag_stages',0)
        if n_dragstages>0:
            dragging=True

        curr_nxlinkbonds=0
        max_search_radius=min(self.TopoCoord.Coordinates.box.diagonal()/2)
        max_radidx=int((max_search_radius-cure_search_radius)/radial_increment)
        cure_finished=False
        CP=Checkpoint(checkpoint_file=checkpoint_file,bonds_file=bonds_file)
        CP.iter=1
        while not CP.state==CPstate.finished:
            cwd=pfs.go_to(f'systems/iter-{CP.iter}')
            CP.read_checkpoint(self)
            if CP.state==CPstate.fresh:
                logging.info(f'CURE iteration {CP.iter}/{maxiter} begins.')
                CP.current_radidx=0 # radius index
                CP.current_stage=0
                CP.set_state(CPstate.bondsearch) # no need to write a checkpoing
            if CP.state==CPstate.bondsearch:
                if dragging:
                    next_stage=CPstate.drag
                else:
                    next_stage=CPstate.update
                CP.radius=cure_search_radius+CP.current_radidx*radial_increment
                while CP.state==CPstate.bondsearch and CP.current_radidx<max_radidx:
                    nbdf=self.make_bonds(CP.radius,header=['ai','aj','reactantName'])
                    if nbdf.shape[0]>0:
                        CP.register_bonds(nbdf,bonds_are='unrelaxed')
                        CP.write_checkpoint(self,next_stage,prefix='0-connect')
                    else:
                        logging.debug(f'CURE iteration {CP.iter}: increasing bondsearch radius to {CP.radius+radial_increment}')
                        CP.current_radidx+=1
                        CP.radius+=radial_increment
                if CP.state==CPstate.bondsearch: # loop exited on radius violation
                    logging.debug(f'CURE iteration {CP.iter} failed to find bonds.')
                    cure_finished=True
            if CP.state==CPstate.drag:
                CP.read_checkpoint(self)
                self.checkout('mdp/drag-em.mdp')
                self.checkout('mdp/drag-nvt.mdp')
                self.checkout('mdp/drag-npt.mdp')
                CP.bonds['initial-distance']=self.TopoCoord.return_bond_lengths(CP.bonds)
                self.TopoCoord.add_restraints(CP.bonds,typ=6)
                begin_dragstage=CP.current_dragstage
                for i in range(begin_dragstage,n_dragstages):
                    self.TopoCoord.attenuate_bond_parameters(CP.bonds,i,n_dragstages,minimum_distance=drag_limit_nm)
                    stagepref=f'1-drag-stage-{i}'
                    CP.write_checkpoint(self,CPstate.drag,prefix=stagepref)
                    msg=grompp_and_mdrun(gro=stagepref,top=stagepref,out=stagepref+'-min',mdp='drag-em',rdd=CP.radius,**self.cfg.parameters)
                    msg=grompp_and_mdrun(gro=stagepref+'-min',top=stagepref,out=stagepref+'-nvt',mdp='drag-nvt',rdd=CP.radius,**self.cfg.parameters)
                    msg=grompp_and_mdrun(gro=stagepref+'-nvt',top=stagepref,out=stagepref+'-npt',mdp='drag-npt',rdd=CP.radius,**self.cfg.parameters)
                    self.TopoCoord.copy_coords(TopoCoord(grofilename=stagepref+'-npt.gro'))
                    current_lengths=np.array(self.TopoCoord.return_bond_lengths(CP.bonds))
                    logging.debug(f'-> avg new pair separation distance: {current_lengths.mean():.3f}')
                    CP.current_dragstage+=1
                CP.current_dragstage-=1
                self.TopoCoord.remove_restraints(CP.bonds)
                CP.write_checkpoint(self,CPstate.update,prefix='1-drag')
            if CP.state==CPstate.update:
                CP.read_checkpoint(self)
                CP.bonds=self.TopoCoord.update_topology_and_coordinates(CP.bonds,template_dict=self.molecules)
                CP.current_stage=0
                CP.write_checkpoint(self,CPstate.relax,prefix='2-update')
            if CP.state==CPstate.relax:
                CP.read_checkpoint(self)
                self.checkout('mdp/relax-em.mdp')
                self.checkout('mdp/relax-nvt.mdp')
                self.checkout('mdp/relax-npt.mdp')
                CP.bonds['initial-distance']=self.TopoCoord.return_bond_lengths(CP.bonds)
                begin_stage=CP.current_stage
                for i in range(begin_stage,n_stages):
                    saveT=self.TopoCoord.copy_bond_parameters(CP.bonds)
                    self.TopoCoord.attenuate_bond_parameters(CP.bonds,i,n_stages)
                    stagepref=f'3-relax-stage-{i}'
                    CP.write_checkpoint(self,CPstate.relax,prefix=stagepref)
                    msg=grompp_and_mdrun(gro=stagepref,top=stagepref,out=stagepref+'-min',mdp='relax-em',rdd=CP.radius,**self.cfg.parameters)
                    msg=grompp_and_mdrun(gro=stagepref+'-min',top=stagepref,out=stagepref+'-nvt',mdp='relax-nvt',rdd=CP.radius,**self.cfg.parameters)
                    msg=grompp_and_mdrun(gro=stagepref+'-nvt',top=stagepref,out=stagepref+'-npt',mdp='relax-npt',rdd=CP.radius,**self.cfg.parameters)
                    self.TopoCoord.copy_coords(TopoCoord(grofilename=stagepref+'-npt.gro'))
                    self.TopoCoord.restore_bond_parameters(saveT)
                    current_lengths=np.array(self.TopoCoord.return_bond_lengths(CP.bonds))
                    logging.debug(f'-> avg new bond length: {current_lengths.mean():.3f}')
                    CP.current_stage+=1
                CP.current_stage-=1
                CP.bonds_are='relaxed'
                CP.write_checkpoint(self,CPstate.equilibrate,prefix='4-equilibrate')
            if CP.state==CPstate.equilibrate:
                self.checkout('mdp/npt-equilibrate.mdp')
                gro,ext=os.path.splitext(CP.gro)
                top,ext=os.path.splitext(CP.top)
                msg=grompp_and_mdrun(gro=gro,top=top,out=gro+'-post',mdp='equilibrate-npt',**self.cfg.parameters)
                self.TopoCoord.copy_coords(TopoCoord(grofilename=gro+'-post.gro'))
                CP.write_checkpoint(self,CPstate.post_equilibration,prefix='final')
            if CP.state==CPstate.post_equilibration:
                curr_nxlinkbonds+=CP.bonds.shape[0]
                curr_conversion=curr_nxlinkbonds/max_nxlinkbonds
                logging.debug(f'Current conversion: {curr_conversion} ({curr_nxlinkbonds}/{max_nxlinkbonds})')
                cure_finished=curr_conversion>desired_conversion
                cure_finished=cure_finished or CP.iter>=maxiter
                CP.set_state(CPstate.fresh)
                CP.iter+=1
            if cure_finished:
                CP.set_state(CPstate.finished)

    def set_system(self,CP=None):
        if CP:
            top=CP.top
            gro=CP.gro
            grx=CP.grx
            logging.debug(f'RESETTING SYSTEM FROM {top} {gro} {grx}')
            self.TopoCoord=TopoCoord(top,gro)
            if (grx):
                self.TopoCoord.read_gro_attributes(grx)
        
    def register_system(self,CP=None,extra_attributes=['z','cycle-idx','reactantName']):
        if CP:
            logging.debug(f'WRITING SYSTEM TO {CP.top} {CP.gro} {CP.grx}')
            self.TopoCoord.write_top(CP.top)
            self.TopoCoord.write_gro(CP.gro)
            self.TopoCoord.write_gro_attributes(extra_attributes,CP.grx)

    def make_bonds(self,radius,header=['ai','aj','reactantName']):
        adf=self.TopoCoord.gro_DataFrame('atoms')
        ncpu=self.cfg.parameters['ncpu']
        self.TopoCoord.linkcell_initialize(radius,ncpu=ncpu)
        # self.Coordinates.linkcell.make_memberlists(self.Coordinates.A)
        raset=adf[adf['z']>0]  # this view will be used for downselecting to potential A-B partners
        #logging.debug(f'make_bonds: there are {raset.shape[0]} reactive atoms')
        # self.Coordinates.show_z_report()
        newbonds=[]
        ''' generate the list of new bonds to make '''
        for R in self.cfg.reactions:
            if R.stage=='post-cure':
                continue
            logging.debug(f'*** BONDS from reaction {R.name}')
            for bond in R.bonds:
                A=R.atoms[bond['atoms'][0]]
                B=R.atoms[bond['atoms'][1]]
                aname=A['atom']
                areactantname_template=R.reactants[A['reactant']]
                aresid_template=A['resid']
                aresname=self.molecules[areactantname_template].get_resname(aresid_template)
                az=A['z']
                bname=B['atom']
                breactantname_template=R.reactants[B['reactant']]
                bresid_template=B['resid']
                bresname=self.molecules[breactantname_template].get_resname(bresid_template)
                bz=B['z']
                Aset=raset[(raset['atomName']==aname)&(raset['resName']==aresname)&(raset['z']==az)&(raset['reactantName']==areactantname_template)]
                Bset=raset[(raset['atomName']==bname)&(raset['resName']==bresname)&(raset['z']==bz)&(raset['reactantName']==breactantname_template)]
                Pbonds=list(product(Aset['globalIdx'].to_list(),Bset['globalIdx'].to_list()))
                logging.debug(f'Examining {Aset.shape[0]}x{Bset.shape[0]}={len(Pbonds)} {aresname}:{aname}({az})-{bresname}:{bname}({bz}) pairs')
                passbonds=[]
                if len(Pbonds)>0:
                    bondtestoutcomes={k:0 for k in BTRC}
                    logging.debug(f'Bond search will use {ncpu} processors')
                    p=Pool(processes=ncpu)
                    Pbonds_split=np.array_split(Pbonds,ncpu)
                    results=p.map(partial(self.TopoCoord.bondtest_par,radius=radius), Pbonds_split)
                    p.close()
                    p.join()
                    rc=[]
                    for i,l in enumerate(results):
                        rc.extend(l)
                    for i,R2 in enumerate(rc):
                        RC,rij=R2 
                        bondtestoutcomes[RC]+=1
                        if RC==BTRC.passed:
                            passbonds.append((Pbonds[i],rij,R.product))
                    logging.debug(f'*** {len(passbonds)} out of {len(Pbonds)} bonds pass initial filter')
                    logging.debug(f'Bond test outcomes:')
                    for k,v in bondtestoutcomes.items():
                        logging.debug(f'   {str(k)}: {v}')
                    newbonds.extend(passbonds)

        ''' Sort new potential bonds by length (ascending) and claim each one
            in order so long as both its atoms are available '''
        newbonds.sort(key=lambda x: x[1])
        logging.debug(f'*** Pruning {len(newbonds)} bonds...')
        atomset=list(set(list([x[0][0] for x in newbonds])+list([x[0][1] for x in newbonds])))
        resid_pairs=[]
        allowed_bond=[True for x in newbonds]
        for k,b in enumerate(newbonds):
            bb,p,t=b 
            i,j=bb
            i_resid=self.TopoCoord.get_gro_attribute_by_attributes('resNum',{'globalIdx':i})
            j_resid=self.TopoCoord.get_gro_attribute_by_attributes('resNum',{'globalIdx':j})
            # need to be careful about intramolecular bonds (capping)
            rp=(i_resid,j_resid)
            if not rp in resid_pairs:
                resid_pairs.append(rp)
                resid_pairs.append(rp[::-1])
            else:
                allowed_bond[k]=False
        # keep shortest bonds up to point atoms are all used up; do not repeat resid pairs
        # so that only one bond is allowed to connect any two residues
        keepbonds=[]
        disallowed=0
        for k,n in enumerate(newbonds):
            if not allowed_bond[k]:
                disallowed+=1
                continue
            b=n[0]
            if b[0] in atomset and b[1] in atomset:
                atomset.remove(b[0])
                atomset.remove(b[1])
                keepbonds.append((b,n[2]))
        logging.debug(f'*** accepted the {len(keepbonds)} shortest non-competing bonds')
        logging.debug(f'    {disallowed} bonds that repeat resid pairs thrown out.')

        kdf=pd.DataFrame({header[0]:[x[0][0] for x in keepbonds],
                          header[1]:[x[0][1] for x in keepbonds],
                          header[2]:[x[1] for x in keepbonds]})
        return kdf

    def initreport(self):
        print(self.cfg)
        print()
        print(self.software)

    def main(self,**kwargs):
        force_parameterization=kwargs.get('force_parameterization',False)
        force_sea_calculation=kwargs.get('force_sea_calculation',False)
        force_checkin=kwargs.get('force_checkin',False)
        self.generate_molecules(
            force_parameterization=force_parameterization,force_sea_calculation=force_sea_calculation,force_checkin=force_checkin
        )
        self.initialize_topology()
        self.setup_liquid_simulation()
        self.do_liquid_simulation()
        self.CURE()
        # self.finalize()

def info():
    print('This is some information on your installed version of HTPolyNet')
    pfs.info()
    software.info()

def cli():
    parser=ap.ArgumentParser()
    parser.add_argument('command',type=str,default=None,help='command (init, info, run)')
    parser.add_argument('-cfg',type=str,default='',help='input config file')
    parser.add_argument('-lib',type=str,default='',help='local library, assumed flat')
    parser.add_argument('-log',type=str,default='htpolynet_runtime.log',help='log file')
    parser.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    parser.add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    parser.add_argument('--force-sea-calculation',default=False,action='store_true',help='force calculation of symmetry-equivalent atoms in any input mol2 structures')
    parser.add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    parser.add_argument('--loglevel',type=str,default='info',help='Log level; info, debug')
    args=parser.parse_args()

    ''' set up logging '''
    loglevel=args.loglevel
    loglevel_numeric=getattr(logging, loglevel.upper())
    logging.basicConfig(filename=args.log,encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=loglevel_numeric)
    logging.info('HTPolyNet runtime begins.')
    ''' set up the project file system and access to HTPolyNet libraries '''
    userlib=None
    if args.lib!='':
        userlib=args.lib
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,reProject=args.restart,userlibrary=userlib)
    software.sw_setup()

    if args.command=='info':
        info()
    elif args.command=='run':
        a=HTPolyNet(cfgfile=args.cfg,restart=args.restart)
        a.main(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization,force_sea_calculation=args.force_sea_calculation)
    else:
        print(f'HTPolyNet command {args.command} not recognized')
    
    logging.info('HTPolynet runtime ends.')
