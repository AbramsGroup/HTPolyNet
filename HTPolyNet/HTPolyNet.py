#!/usr/bin/env python
"""
@author: huang, abrams
"""
import logging
import os
import shutil
import yaml
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
        new_molecules=self.cfg.symmetry_expand_reactions()
        for mname,M in new_molecules.items():
            ''' a molecule with no generator specified is a monomer and can be generated by
                parameterizing an existing mol2 file or reading in a previous parameterization
                from the library '''
            if mname not in self.molecules:
                self.generate_molecule(M,force_parameterization=force_parameterization,force_sea_calculation=force_sea_calculation,force_checkin=force_checkin)
                assert M.get_origin()!='unparameterized'
                self.molecules[mname]=M
                logging.info(f'Generated {mname}')
        self.molecules.update(new_molecules)
        self.cfg.maxconv=self.cfg.calculate_maximum_conversion()
        logging.info(f'Maximum conversion is {self.cfg.maxconv} bonds.')
        precursors=[M for M in self.molecules.values() if not M.generator]
        for M in precursors:
            M.propagate_z(self.cfg.reactions,self.molecules)
        reaction_products=[M for M in self.molecules.values() if M.generator]
        for M in reaction_products:
            M.propagate_z(self.cfg.reactions,self.molecules)
        for M in self.molecules.values():
            logging.debug(f'Ring detector for {M.name}')
            M.label_ring_atoms(M.TopoCoord.ring_detector())
        for M in self.molecules.values():
            if len(M.reaction_bonds)>0:
                logging.debug(f'Crosslink bond(s) in {M.name}:')
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
                # logging.info(f'Generating {mname}')
                M.generate(available_molecules=self.molecules,**self.cfg.parameters)
                for ex in ['mol2','top','itp','gro']:
                    checkin(f'molecules/parameterized/{mname}.{ex}',overwrite=force_checkin)
                M.set_origin('newly parameterized')
            else:
                return
        else:
            logging.info(f'Fetching parameterized {mname}')
            for ex in ['mol2','top','itp','gro']:
                self.checkout(f'molecules/parameterized/{mname}.{ex}')
            M.load_top_gro(f'{mname}.top',f'{mname}.gro')
            M.set_sequence()
            M.set_reaction_bonds(self.molecules)
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
            # we have to transfer from the base reactants
            M.inherit_attribute_from_reactants('sea-idx',available_molecules=self.molecules)
            M.write_gro_attributes(['sea-idx'],f'{M.name}.sea')
        return True

    def initialize_global_topology(self,inpfnm='init'):
        ''' Create a full gromacs topology that includes all directives necessary 
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined. '''
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
            logging.info(f'initialize_topology merging {N} copies of {M.name} into global topology')
            self.TopoCoord.Topology.merge(t)
            already_merged.append(M.name)
        for othermol,M in self.molecules.items():
            if not othermol in already_merged:
                self.TopoCoord.Topology.merge_types(M.TopoCoord.Topology)
        logging.info(f'Extended topology has {self.TopoCoord.Topology.atomcount()} atoms.')
        self.TopoCoord.write_top(f'{inpfnm}.top')
        logging.info(f'Wrote {inpfnm}.top to {cwd}')

    def setup_liquid_simulation(self,inpfnm='init'):
        # go to the results path, make the directory 'init', cd into it
        cwd=pfs.go_to('systems/init')
        # fetch unreacted init.top amd all monomer gro's 
        # from parameterization directory
        # self.checkout(f'{inpfnm}.top',altpath=pfs.subpath('systems'))
        for n in self.cfg.molecules.keys():
            self.checkout(f'molecules/parameterized/{n}.gro',altpath=pfs.subpath('molecules'))
        self.checkout('mdp/em.mdp')
        self.checkout('mdp/npt-1.mdp')
        if 'initial_boxsize' in self.cfg.parameters:
            boxsize=self.cfg.parameters['initial_boxsize']
        elif 'initial_density' in self.cfg.parameters:
            mass_kg=self.TopoCoord.total_mass(units='SI')
            V0_m3=mass_kg/self.cfg.parameters['initial_density']
            L0_m=V0_m3**(1./3.)
            L0_nm=L0_m*1.e9
            logging.info(f'Initial density {self.cfg.parameters["initial_density"]} kg/m^3 and total mass {mass_kg:.3f} kg dictate an initial box side length of {L0_nm:.3f} nm')
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
        for m,M in m_togromacs.items():
            logging.info(f'Molecule to gromacs: {m} ({M.name})')
        for m,c in c_togromacs.items():
            logging.info(f'Composition to gromacs: {m} {c}')
        if not os.path.exists(f'{inpfnm}.gro') or not os.path.exists(f'{inpfnm}.top'): 
            msg=insert_molecules(m_togromacs,c_togromacs,boxsize,inpfnm)
            logging.info(f'Generated {inpfnm}.top and {inpfnm}.gro.')
        else:
            logging.info(f'Found {inpfnm}.gro.')
        self.TopoCoord.read_gro(f'{inpfnm}.gro')
        self.TopoCoord.atom_count()
        self.TopoCoord.inherit_attributes_from_molecules(['z','cycle-idx'],self.cfg.molecules)
        self.TopoCoord.write_gro_attributes(['z','cycle-idx'],f'{inpfnm}.grx')
        self.TopoCoord.make_ringlist()
        self.TopoCoord.make_resid_graph()

    def do_liquid_simulation(self,inpfnm='init',deffnm='npt-1'):
        cwd=pfs.go_to('systems/init')
        if os.path.exists(f'{deffnm}.gro') and self.cfg.parameters['restart']:
            logging.info(f'{deffnm}.gro exists in {os.getcwd()}; skipping initial NPT md.')
        else:
            msg=grompp_and_mdrun(gro=inpfnm,top=inpfnm,out=f'{inpfnm}-min-1',mdp='em')
            msg=grompp_and_mdrun(gro=f'{inpfnm}-min-1',top=inpfnm,out=deffnm,mdp=deffnm)
            logging.info(f'Generated configuration {deffnm}.gro\n')
        density_trace(deffnm)
        # update coordinates
        self.TopoCoord.copy_coords(TopoCoord(grofilename=f'{deffnm}.gro'))

    def SCUR(self):
        # Search - Connect - Update - Relax
        restart=self.cfg.parameters['restart']
        force_scur=self.cfg.parameters.get('force_scur',False)
        logging.info('*'*10+' SEARCH - CONNECT - UPDATE - RELAX  BEGINS '+'*'*10)
        max_nxlinkbonds=self.cfg.maxconv # only for a stoichiometric system
        logging.debug(f'SCUR: max_nxlinkbonds {max_nxlinkbonds}')
        if max_nxlinkbonds==0:
            logging.warning(f'Apparently there are no crosslink bonds to be made! (sum of z == 0)')
            return
        scur_finished=False
        scur_search_radius=self.cfg.parameters['SCUR_cutoff']
        desired_conversion=self.cfg.parameters['conversion']
        radial_increment=self.cfg.parameters.get('SCUR_radial_increment',0.5)
        maxiter=self.cfg.parameters.get('max_SCUR_iterations',20)
        iter=0
        curr_nxlinkbonds=0
        while not scur_finished:
            cwd=pfs.go_to(f'systems/iter{iter}')
            self.TopoCoord.show_z_report()
            num_newbonds=self.scur_iter_complete_check(iter)
            if num_newbonds>0:
                logging.info(f'SCUR iteration {iter+1}/{maxiter} already ran.')
            else:
                logging.info(f'SCUR iteration {iter+1}/{maxiter} begins in {cwd}')
                newbonds=self.scur_make_bonds(scur_search_radius,iter,force_scur=force_scur)
                if len(newbonds)>0:
                    top,gro,newbonds=self.update_topology_and_coordinates(newbonds,iter)
                    top,gro,grx=self.gromacs_stepwise_relaxation(newbonds,top,gro)
                    # self.Coordinates.show_z_report()
                num_newbonds=len(newbonds)
                logging.info(f'SCUR iteration {iter+1}/{maxiter} ends.')
                self.scur_iter_complete_register(iter,top,gro,grx,num_newbonds)

            curr_nxlinkbonds+=num_newbonds
            curr_conversion=curr_nxlinkbonds/max_nxlinkbonds
            scur_finished=curr_conversion>desired_conversion
            scur_finished=scur_finished or (iter+1)>=maxiter
            # logging.debug(f'iter {iter} maxiter {maxiter} iter>=maxiter {iter>=maxiter}')
            if not scur_finished:
                if num_newbonds==0:
                    logging.info(f'No new bonds in SCUR iteration {iter}')
                    logging.info(f'-> updating search radius to {scur_search_radius}')
                    scur_search_radius += radial_increment
                    # TODO: prevent search radius from getting too large for box
                    logging.info(f'-> updating search radius to {scur_search_radius}')
            logging.info(f'Current conversion: {curr_conversion}')
            logging.info(f'   SCUR finished? {scur_finished}')
            iter+=1
        logging.info(f'SCUR iterations finished.')
        # TODO: any post-cure reactions handled here

    def scur_iter_complete_register(self,iter,top,gro,grx,num_newbonds):
        with open('complete.yaml','w') as f:
            f.write(f'ITERATION: {iter+1}\n')
            f.write(f'TOPOLOGY: {top}\n')
            f.write(f'COORDINATES: {gro}\n')
            f.write(f'EXTRA_ATTRIBUTES: {grx}\n')
            f.write(f'NUM_NEWBONDS: {num_newbonds}\n')
            f.close()

    def scur_iter_complete_check(self,iter):
        num_newbonds=0
        if os.path.exists('complete.yaml'):
            with open('complete.yaml','r') as f:
                basedict=yaml.safe_load(f)
                top=os.path.basename(basedict['TOPOLOGY'])
                gro=os.path.basename(basedict['COORDINATES'])
                grx=os.path.basename(basedict['EXTRA_ATTRIBUTES'])
                self.set_system(top,gro,grx=grx)
                num_newbonds=basedict['NUM_NEWBONDS']
            logging.info(f'SCUR iteration {iter+1} marked complete:')
            logging.info(f'        topology    {top}')
            logging.info(f'        coordinates {gro}')
            logging.info(f'        extra_attr  {grx}')
            logging.info(f'        #newbonds   {num_newbonds}')
        return num_newbonds

    def set_system(self,top,gro,grx=None):
        logging.debug(f'RESETTING SYSTEM FROM {top} {gro} {grx}')
        self.TopoCoord=TopoCoord(top,gro)
        if (grx):
            self.TopoCoord.read_gro_attributes(grx)
            
    def scur_make_bonds(self,radius,iter,force_scur=False):
        if os.path.exists(f'keepbonds-{iter}.dat') and not force_scur:
            logging.debug(f'keepbonds-{iter}.dat found.  Reading.')
            with open (f'keepbonds-{iter}.dat','r') as f:
                lines=f.read().split('\n')
                keepbonds=[]
                for line in lines:
                    t=line.split()
                    b=(int(t[0]),int(t[1]))
                    keepbonds.append((b,t[2]))
                return keepbonds
        adf=self.TopoCoord.gro_DataFrame('atoms')
        self.TopoCoord.linkcell_initialize(radius)
        # self.Coordinates.linkcell.make_memberlists(self.Coordinates.A)
        raset=adf[adf['z']>0]  # this view will be used for downselecting to potential A-B partners
        logging.debug(f'make_bonds: there are {raset.shape[0]} reactive atoms')
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
                logging.debug(f'  -> bond\n    {R.reactants}\n    {bond}:\n    A {A}\n    B {B}')
                aname=A['atom']
                aresname_template=R.reactants[A['reactant']]
                aresid_template=A['resid']
                # get ACTUAL resname from the resid'th residue of reactant...
                aresname=self.molecules[aresname_template].get_resname(aresid_template)
                az=A['z']
                bname=B['atom']
                bresname_template=R.reactants[B['reactant']]
                bresid_template=B['resid']
                # get ACTUAL resname from the resid'th residue of reactant...
                bresname=self.molecules[bresname_template].get_resname(bresid_template)
                bz=B['z']
                Aset=raset[(raset['atomName']==aname)&(raset['resName']==aresname)&(raset['z']==az)]
                Bset=raset[(raset['atomName']==bname)&(raset['resName']==bresname)&(raset['z']==bz)]
                logging.debug(f'Aset.shape[0] {Aset.shape[0]}')
                logging.debug(f'Bset.shape[0] {Bset.shape[0]}')
                Pbonds=list(product(Aset['globalIdx'].to_list(),Bset['globalIdx'].to_list()))
                logging.debug(f'Examining {len(Pbonds)} {aresname}:{aname}({az})-{bresname}:{bname}({bz}) pairs')
                passbonds=[]
                if len(Pbonds)>0:
                    bondtestoutcomes={k:0 for k in BTRC}
                    logging.debug(f'Bond search will use {self.cfg.parameters["cpu"]} processes')
                    p = Pool(processes=self.cfg.parameters['cpu'])
                    Pbonds_split = np.array_split(Pbonds,self.cfg.parameters['cpu'])
                    results = p.map(partial(self.TopoCoord.bondtest_par,radius=radius), Pbonds_split)
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
        # keep shortest bonds up to point atoms are all used up
        keepbonds=[]
        for n in newbonds:
            b=n[0]
            if b[0] in atomset and b[1] in atomset:
                atomset.remove(b[0])
                atomset.remove(b[1])
                keepbonds.append((b,n[2]))
        logging.debug(f'*** accepted the {len(keepbonds)} shortest non-competing bonds')
        with open(f'keepbonds-{iter}.dat','w') as f:
            f.write('\n'.join(f'{i[0][0]} {i[0][1]} {i[1]}' for i in keepbonds))
        return keepbonds

    def update_topology_and_coordinates(self,keepbonds,iter):    
        ''' make the new bonds '''
        logging.debug(f'update_topology_and_coordinates begins.')
        if len(keepbonds)>0:
            pairs=[i[0] for i in keepbonds]
            # logging.debug(f'update_topo_coords: bondlist: {bondlist}')
            logging.debug(f'Making {len(pairs)} bonds.')
            idx_to_delete=self.TopoCoord.make_bonds(pairs)
            logging.debug(f'Deleting {len(idx_to_delete)} atoms.')
            idx_mapper=self.TopoCoord.delete_atoms(idx_to_delete) # will result in full reindexing
            reindexed_keepbonds=[((idx_mapper[i[0][0]],idx_mapper[i[0][1]]),i[1]) for i in keepbonds]
            pairs=[i[0] for i in reindexed_keepbonds]
            self.TopoCoord.decrement_z(pairs)
            self.TopoCoord.make_ringlist()
            self.TopoCoord.map_from_templates(reindexed_keepbonds,self.molecules)
            self.TopoCoord.adjust_charges(msg='You might want to increase the scope of template mapping for each new bond.')
            basefilename=f'scur-step-{iter}'
            self.TopoCoord.write_top_gro(basefilename+'.top',basefilename+'.gro')
            logging.debug(f'Wrote {basefilename}.top and {basefilename}.gro.')
            return (basefilename+'.top',basefilename+'.gro',reindexed_keepbonds)

    def gromacs_stepwise_relaxation(self,newbonds,fulltop,initcoords):
        self.checkout('mdp/em-inter-scur-relax-stage.mdp')
        self.checkout('mdp/nvt-inter-scur-relax-stage.mdp')
        self.checkout('mdp/npt-inter-scur-relax-stage.mdp')
        self.checkout('mdp/npt-inter-scur-iter.mdp')
        tmpTC=TopoCoord(fulltop,initcoords)
        pref,ext=os.path.splitext(fulltop)
        n_stages=self.cfg.parameters.get('max_bond_relaxation_stages',6)
        bonds=[b[0] for b in newbonds] # strip the product template names
        lengths=self.TopoCoord.return_bond_lengths(bonds)
        for i in range(n_stages):
            saveT=tmpTC.copy_bond_parameters(bonds)
            tmpTC.attenuate_bond_parameters(bonds,i,n_stages,lengths)
            stagepref=pref+f'-stage-{i}'
            tmpTC.write_top_gro(stagepref+'.top',stagepref+'.gro')
            msg=grompp_and_mdrun(gro=stagepref,top=stagepref,out=stagepref+'-min',mdp='em-inter-scur-relax-stage')
            msg=grompp_and_mdrun(gro=stagepref+'-min',top=stagepref,out=stagepref+'-nvt',mdp='nvt-inter-scur-relax-stage')
            msg=grompp_and_mdrun(gro=stagepref+'-min',top=stagepref,out=stagepref+'-npt',mdp='npt-inter-scur-relax-stage')
            sacmol=TopoCoord(grofilename=stagepref+'-npt.gro')
            tmpTC.copy_coords(sacmol)
            tmpTC.restore_bond_parameters(saveT)
        msg=grompp_and_mdrun(gro=stagepref+'-npt',top=pref,out=pref+'-post',mdp='npt-inter-scur-iter')
        sacmol=TopoCoord(grofilename=pref+'-post.gro')
        self.TopoCoord.copy_coords(sacmol)
        self.TopoCoord.write_gro_attributes(['z','cycle-idx'],pref+'-post.grx')
        return fulltop,pref+'-post.gro',pref+'-post.grx'

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
        self.initialize_global_topology()
        self.setup_liquid_simulation()
        self.do_liquid_simulation()
        self.SCUR()
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
