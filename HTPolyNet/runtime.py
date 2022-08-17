import logging
import os
import shutil
from termios import N_SLIP
import numpy as np
import random
from copy import deepcopy
from HTPolyNet.configuration import Configuration
from HTPolyNet.topology import select_topology_type_option
from HTPolyNet.topocoord import TopoCoord
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.gromacs import insert_molecules, mdp_modify, mdp_get
import HTPolyNet.checkpoint as cp
from HTPolyNet.plot import trace
from HTPolyNet.molecule import Molecule, MoleculeDict, is_reactant
from HTPolyNet.expandreactions import symmetry_expand_reactions, chain_expand_reactions
from HTPolyNet.curecontroller import CureController
from HTPolyNet.stringthings import my_logger

logger=logging.getLogger(__name__)

def logrotate(filename):
    if os.path.exists(filename):
        n=1
        while os.path.exists(f'#{n}#{filename}'):
            n+=1
        shutil.copyfile(filename,f'#{n}#{filename}')

_directives=['ps','nsteps','ncycles']
def _nonempty_directives(dirlist):
    amts=[]
    logger.debug(f'dirlist {dirlist}')
    for d in dirlist:
        ta=[d.get(x,0) for x in _directives]
        amts.append(any(ta))
        logger.debug(f'  dir {d} amts {amts}')
    logger.debug(f'returning {any(amts)}')
    return any(amts)

class Runtime:
    runtime_defaults={
        'gromacs': {
            'gmx': 'gmx',
            'gmx_options': '-quiet -nobackup',
            'mdrun': 'gmx mdrun'
        },
        'ambertools': {
            'charge_method': 'gas'
        },
        'densification': {
            'initial_density': 200.0,  # kg/m3
            'aspect_ratio': np.array([1.0,1.0,1.0]),
            'equilibration': [
                { 'ensemble': 'min' },
                { 'ensemble': 'nvt', 'temperature': 300, 'ps': 10 },
                { 'ensemble': 'npt', 'temperature': 300, 'pressure': 10, 'ps': 200 }
            ]
        },
        'precure': {  
            'preequilibration': {
                'ensemble': 'npt',
                'temperature': 300,        # K
                'pressure': 1,             # bar
                'ps': 100  # we should probably bring to a desired pressure after densification
            },
            'anneal': {
                'ncycles': 0,  # default none
                'initial_temperature': 300,
                'cycle_segments': [
                    { 'T': 300, 'ps': 0 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 300, 'ps': 20 },
                    { 'T': 300, 'ps': 20 }
                ]
            },
            'postequilibration': {
                'ensemble': 'npt',
                'temperature': 300,        # K
                'pressure': 1,             # bar
                'ps': 0 # default none
            }
        },
        'CURE': {},
        'postcure': {
            'anneal': {
                'ncycles': 0, # default none
                'initial_temperature': 300,
                'cycle_segments': [
                    { 'T': 300, 'ps': 0 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 600, 'ps': 20 },
                    { 'T': 300, 'ps': 20 },
                    { 'T': 300, 'ps': 20 }
                ]
            },
            'postequilibration': {
                'ensemble': 'npt',
                'temperature': 300,       # K
                'pressure': 1,            # bar
                'ps':  0 # default none
            }
        }
    }
    ''' Class for a single HTPolyNet runtime session '''
    def __init__(self,cfgfile='',restart=False):
        my_logger(software.to_string(),logger.info)
        self.cfgfile=cfgfile
        if cfgfile=='':
            logger.error('HTPolyNet requires a configuration file.\n')
            raise RuntimeError('HTPolyNet requires a configuration file')
        logger.info(f'Configuration: {cfgfile}')
        self.cfg=Configuration.read(os.path.join(pfs.root(),cfgfile))
        # add any missing defaults
        for k,default_values in self.runtime_defaults.items():
            if not k in self.cfg.parameters:
                self.cfg.parameters[k]=default_values
            if type(default_values)==dict:
                for kk,vv in default_values.items():
                    if not kk in self.cfg.parameters[k]:
                        self.cfg.parameters[k][kk]=vv
        software.set_gmx_preferences(self.cfg.parameters)
        self.TopoCoord=TopoCoord(system_name='htpolynet')
        self.cfg.parameters['restart']=restart
        if self.cfg.parameters['restart']:
            logger.info(f'Restarting in {pfs.proj()}')
        self.molecules:MoleculeDict={}
        cure_dict=self.cfg.parameters.get('CURE',{})
        if cure_dict:
            self.cc=CureController(cure_dict)
        self.ncpu=self.cfg.parameters.get('ncpu',os.cpu_count())

    def generate_molecules(self,force_parameterization=False,force_checkin=False):
        GAFF_dict=self.cfg.parameters.get('GAFF',{})
        self.molecules={}
        ''' configuration.parse() generated a list of Molecules implied by configuration; assume
            they are all unparameterized '''
        for mname,M in self.cfg.molecules.items():
            M.set_origin('unparameterized')
        cwd=pfs.go_to('molecules/parameterized')
        my_logger(f'Templates in {pfs.cwd()}',logger.info)
        ''' Each molecule implied by the cfg is 'generated' here, either by
            reading from the library or direct parameterization.  In some cases,
            the molecule is to be generated by a reaction; if so, it's
            `generator` attribute will be a Reaction instance '''
        ess='' if len(self.cfg.molecules)==1 else 's'
        logger.info(f'{len(self.cfg.molecules)} molecule{ess} explicit in {self.cfgfile}')
        # ml=list(self.cfg.molecules.keys())
        # my_logger(ml,logger.info)
        all_made=all([(m in self.molecules and M.get_origin()!='unparameterized') for m,M in self.cfg.molecules.items()])
        while not all_made:
            for mname,M in self.cfg.molecules.items():
                if mname not in self.molecules:
                    if self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin):
                        self.molecules[mname]=M
                        logger.debug(f'Generated {mname}')
            all_made=all([(m in self.molecules and M.get_origin()!='unparameterized') for m,M in self.cfg.molecules.items()])

        ''' Generate all reactions and products that result from invoking symmetry '''
        symmetry_relateds=self.cfg.parameters.get('symmetry_equivalent_atoms',{})
        if not symmetry_relateds:
            constituents=self.cfg.parameters.get('constituents',{})
            # if not constituents:
            #     raise Exception(f'Config file must have a "symmetry_equivalent_atoms" key if no "constituents" key is specified')
            for cname,crec in constituents.items():
                this_sr=crec.get('symmetry_equivalent_atoms',[])
                if len(this_sr)>0:
                    symmetry_relateds[cname]=this_sr
        if len(symmetry_relateds)>0:
            new_reactions,new_molecules=symmetry_expand_reactions(self.cfg.reactions,symmetry_relateds)
            ess='' if len(new_molecules)==1 else 's'
            logger.info(f'{len(new_molecules)} molecule{ess} implied by symmetry-equivalent atoms')
            ml=list(new_molecules.keys())
            logger.info(ml)
            self.cfg.reactions.extend(new_reactions)
            make_molecules={k:v for k,v in new_molecules.items() if k not in self.molecules}
            for mname,M in make_molecules.items():
                self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
                assert M.get_origin()!='unparameterized'
                self.molecules[mname]=M
                logger.debug(f'Generated {mname}')

        ''' Generate any required template products that result from reactions in which the bond generated creates
            dihedrals that span more than just the two monomers that are connected '''
        new_reactions,new_molecules=chain_expand_reactions(self.molecules)
        if len(new_molecules)>0:
            ess='' if len(new_molecules)==1 else 's'
            logger.info(f'{len(new_molecules)} molecule{ess} implied by chaining')
            ml=list(new_molecules.keys())
            logger.info(ml)
            self.cfg.reactions.extend(new_reactions)
            make_molecules={k:v for k,v in new_molecules.items() if k not in self.molecules}
            for mname,M in make_molecules.items():
                # logger.debug(f'Generating {mname}:')
                self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
                assert M.get_origin()!='unparameterized'
                self.molecules[mname]=M
                logger.debug(f'Generated {mname}')

        for M in self.molecules:
            self.molecules[M].is_reactant=is_reactant(M,self.cfg.reactions,stage='cure')

        resolve_type_discrepancies=GAFF_dict.get('resolve_type_discrepancies',[])
        if not resolve_type_discrepancies:
            resolve_type_discrepancies=self.cfg.parameters.get('resolve_type_discrepancies',[])
        if resolve_type_discrepancies:
            for resolve_directive in resolve_type_discrepancies:
                logger.info(f'Resolving type discrepancies for directive {resolve_directive}')
                self._type_consistency_check(typename=resolve_directive['typename'],funcidx=resolve_directive.get('funcidx',4),selection_rule=resolve_directive['rule'])

        ess='' if len(self.molecules)==1 else 's'
        logger.info(f'Generated {len(self.molecules)} molecule template{ess}')
        if self.cfg.initial_composition: 
            logger.info(f'Initial composition is {", ".join([(x["molecule"]+" "+str(x["count"])) for x in self.cfg.initial_composition])}')
            self.cfg.calculate_maximum_conversion()
            logger.info(f'100% conversion is {self.cfg.maxconv} bonds')

        # for mname,M in self.molecules.items():
        #     if M.generator:
        #         M.prepare_new_bonds(available_molecules=self.molecules)

        logger.debug(f'Reaction bond(s) in each molecular template:')
        for M in self.molecules.values():
            if len(M.reaction_bonds)>0:
                logger.debug(f'Bond {M.name}:')
                for b in M.reaction_bonds:
                    logger.debug(f'   {str(b)}')

        logger.debug(f'Bond template(s) in each molecular template:')
        for M in self.molecules.values():
            if len(M.bond_templates)>0:
                logger.debug(f'Template {M.name}:')
                for b in M.bond_templates:
                    logger.debug(f'   {str(b)}')

        for k,v in self.cfg.parameters.get('constituents',{}).items():
            relaxdict=v.get('relax',{})
            if relaxdict:
                self.molecules[k].relax(relaxdict)

    def do_initialization(self,inpfnm='init'):
        cwd=pfs.go_to('systems/init')
        my_logger(f'Initialization in {pfs.cwd()}',logger.info)
        self._initialize_topology(inpfnm)
        self._initialize_coordinates(inpfnm)

    def do_densification(self,deffnm='densified'):
        """do_liquid_simulation Manages execution of gmx mdrun to perform minimization
            and NPT MD simulation of the initial liquid system.  Final coordinates are
            loaded into the global TopoCoord.

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        :param deffnm: deffnm prefix fed to gmx mdrun, defaults to 'npt-1'
        :type deffnm: str, optional
        """
        if cp.passed('do_densification'): return
        cwd=pfs.go_to('systems/densification')
        my_logger(f'Densification in {pfs.cwd()}',logger.info)
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        densification_dict=self.cfg.parameters.get('densification',{})
        assert len(densification_dict)>0,'"densification" directives missing'
        equilibration=densification_dict.get('equilibration',{})
        assert len(equilibration)>0,'equilibration directives missing'
        TC=self.TopoCoord
        infiles=[TC.files[x] for x in ['gro','top','grx']]
        assert all([os.path.exists(x) for x in infiles]),f'One or more of {infiles} not found'
        self._do_equilibration_series(equilibration,deffnm=f'{deffnm}',plot_pfx='densification')
        gro_res=os.path.basename(TC.files["gro"])
        deffnm_res=gro_res.replace('.gro','')
        logger.info(f'Densified coordinates in {pfs.cwd()}/{gro_res}')
        cp.set(TC,'do_densification')

    def do_precure(self):
        if cp.passed('do_precure'): return
        self._do_pap(pfx='precure')
        cp.set(self.TopoCoord,'do_precure')

    def do_cure(self):
        if not hasattr(self,'cc'): return  # no cure controller
        if cp.passed('cure'): return
        cc=self.cc
        TC=self.TopoCoord
        RL=self.cfg.reactions
        MD=self.molecules
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        if cp.is_currentstepname('cure'):
            cc.iter=cp.last_substep()
        else:
            cc.reset()
        cc.setup(max_nxlinkbonds=self.cfg.maxconv,desired_nxlinkbonds=int(self.cfg.maxconv*cc.dicts['controls']['desired_conversion']),max_search_radius=min(TC.Coordinates.box.diagonal()/2))
        cure_finished=cc.is_cured()
        if cure_finished: return
        my_logger('Connect-Update-Relax-Equilibrate (CURE) begins',logger.info)
        logger.info(f'Attempting to form {cc.desired_nxlinkbonds} bonds')
        while not cure_finished:
            my_logger(f'Iteration {cc.iter} begins',logger.info,fill='~')
            reentry=pfs.go_to(f'systems/iter-{cc.iter}')
            if os.path.exists('cure_controller_state.yaml'):
                logger.debug(f'Reading new cure controller in {pfs.cwd()}')
                self.cc=CureController.from_yaml('cure_controller_state.yaml')
                cc=self.cc
                logger.info(f'Restarting at {cc.cum_nxlinkbonds} bonds')
            TC.grab_files() # copy files locally
            cc.do_bondsearch(TC,RL,MD,reentry=reentry)
            cc.do_preupdate_dragging(TC)
            cc.do_topology_update(TC,MD)
            cc.do_relax(TC)
            cc.do_equilibrate(TC,gromacs_dict)
            cp.subset(TC,'cure',cc.iter)
            logger.info(f'Iteration {cc.iter} current conversion {cc.curr_conversion():.3f} or {cc.cum_nxlinkbonds} bonds')
            cure_finished=cc.is_cured()
            if not cure_finished:
                cure_finished=cc.next_iter()
        my_logger(f'Capping begins',logger.info)
        cc.do_cap_bondsearch(TC,RL,MD)
        cc.do_topology_update(TC,MD)
        cc.do_relax(TC)
        cc.do_equilibrate(TC,gromacs_dict)
        my_logger('Connect-Update-Relax-Equilibrate (CURE) ends',logger.info)
        cp.set(TC,'cure')

    def do_postcure(self):
        if cp.passed('do_postcure'): return
        self._do_pap(pfx='postcure')
        cp.set(self.TopoCoord,'do_postcure')
        
    def save_data(self,result_name='final'):
        TC=self.TopoCoord
        pfs.go_to(f'systems/{result_name}-results')
        my_logger(f'Final data to {pfs.cwd()}',logger.info)
        TC.write_grx_attributes(f'{result_name}.grx')
        TC.write_gro(f'{result_name}.gro')
        TC.write_top(f'{result_name}.top')
        cp.set(TC,'final')

    def build(self,**kwargs):
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        checkpoint_file=kwargs.get('checkpoint_file','checkpoint_state.yaml')
        TC=self.TopoCoord

        pfs.go_proj()

        self.generate_molecules(
            force_parameterization=force_parameterization,  # force antechamber/GAFF parameterization
            force_checkin=force_checkin                     # force check-in to system libraries
        )

        pfs.go_proj()
        cp.setup(TC,filename=checkpoint_file)
        TC.load_files()
        self.do_initialization()
        self.do_densification()
        self.do_precure()
        self.do_cure()
        self.do_postcure()
        self.save_data()

    def _generate_molecule(self,M:Molecule,**kwargs):
        mname=M.name
        checkin=pfs.checkin
        # pfs.go_to(f'molecules/parameterized/work/{M.name}')
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        if force_parameterization or not M.previously_parameterized():
            logger.debug(f'Parameterization of {mname} requested -- can we generate {mname}?')
            generatable=(not M.generator) or (all([m in self.molecules for m in M.generator.reactants.values()]))
            if generatable:
                logger.debug(f'Generating {mname}')
                M.generate(available_molecules=self.molecules,**self.cfg.parameters)
                for ex in ['mol2','top','itp','gro','grx']:
                    checkin(f'molecules/parameterized/{mname}.{ex}',overwrite=force_checkin)
#                    checkin(f'molecules/parameterized/work/{mname}/{mname}.{ex}',overwrite=force_checkin)
#                    shutil.copy(f'molecules/parameterized/{mname}.{ex}','../../')
                M.set_origin('newly parameterized')
            else:
                logger.debug(f'...no, did not generate {mname}')
                logger.debug(f'not ({mname}.generator) {bool(not M.generator)}')
                if M.generator:
                    logger.debug(f'reactants {list(M.generator.reactants.values())}')
                return False
        else:
            logger.debug(f'Fetching parameterized {mname}')
            exts=pfs.fetch_molecule_files(mname)
            logger.debug(f'fetched {mname} exts {exts}')
            # for ex in ['mol2','top','itp','gro','grx']:
            #     pfs.checkout(f'molecules/parameterized/{mname}.{ex}')
            mol2=f'{mname}.mol2' if 'mol2' in exts else ''
            M.load_top_gro(f'{mname}.top',f'{mname}.gro',mol2filename='',wrap_coords=False)
            M.TopoCoord.read_gro_attributes(f'{mname}.grx')
            # logger.debug(f'{M.name} box {M.TopoCoord.Coordinates.box}')
            M.set_sequence()
            # if M.generator:
            #     M.prepare_new_bonds(available_molecules=self.molecules)
            M.set_origin('previously parameterized')

        M.generate_stereoisomers()
        M.generate_conformers(minimize=True)

        return True

    def _type_consistency_check(self,typename='dihedraltypes',funcidx=4,selection_rule='stiffest'):
        logger.debug(f'Consistency check of {typename} func {funcidx} on all {len(self.molecules)} molecules requested')
        mnames=list(self.molecules.keys())
        # checkin=pfs.checkin
        types_duplicated=[]
        for i in range(len(mnames)):
            logger.debug(f'{mnames[i]}...')
            for j in range(i+1,len(mnames)):
                mol1topo=self.molecules[mnames[i]].TopoCoord.Topology
                mol2topo=self.molecules[mnames[j]].TopoCoord.Topology
                this_duptypes=mol1topo.report_duplicate_types(mol2topo,typename=typename,funcidx=funcidx)
                logger.debug(f'...{mnames[j]} {len(this_duptypes)}')
                for d in this_duptypes:
                    if not d in types_duplicated:
                        types_duplicated.append(d)
        logger.debug(f'Duplicated {typename}: {types_duplicated}')
        options={}
        for t in types_duplicated:
            logger.debug(f'Duplicate {t}:')
            options[t]=[]
            for i in range(len(mnames)):
                moltopo=self.molecules[mnames[i]].TopoCoord.Topology
                this_type=moltopo.report_type(t,typename='dihedraltypes',funcidx=4)
                if len(this_type)>0 and not this_type in options[t]:
                    options[t].append(this_type)
                logger.debug(f'{mnames[i]} reports {this_type}')
            logger.debug(f'Conflicting options for this type: {options[t]}')
            selected_type=select_topology_type_option(options[t],typename,rule=selection_rule)
            logger.debug(f'Under selection rule "{selection_rule}", preferred type is {selected_type}')
            for i in range(len(mnames)):
                logger.debug(f'resetting {mnames[i]}')
                TC=self.molecules[mnames[i]].TopoCoord
                moltopo=TC.Topology
                moltopo.reset_type(typename,t,selected_type)

    def _initialize_topology(self,inpfnm='init'):
        """Create a full gromacs topology that includes all directives necessary
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined.

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        """
        TC=self.TopoCoord
        if cp.passed('initialize_topology'): return
        cwd=pfs.go_to('systems/init')
        if os.path.isfile(f'{inpfnm}.top'):
            logger.debug(f'{inpfnm}.top already exists in {cwd} but we will rebuild it anyway!')

        ''' for each monomer named in the cfg, either parameterize it or fetch its parameterization '''
        already_merged=[]
        for item in self.cfg.initial_composition:
            M=self.molecules[item['molecule']]
            N=item['count']
            t=deepcopy(M.TopoCoord.Topology)
            logger.debug(f'Merging {N} copies of {M.name}\'s topology into global topology')
            t.adjust_charges(atoms=t.D['atoms']['nr'].to_list(),desired_charge=0.0,overcharge_threshhold=0.1,msg='')
            t.rep_ex(N)
            TC.Topology.merge(t)
            already_merged.append(M.name)
        for othermol,M in self.molecules.items():
            if not othermol in already_merged:
                logger.debug(f'Merging types from {othermol}\'s topology into global topology')
                TC.Topology.merge_types(M.TopoCoord.Topology)
        logger.info(f'Topology "{inpfnm}.top" in {pfs.cwd()}')
        TC.write_top(f'{inpfnm}.top')
        cp.set(self.TopoCoord,stepname='initialize_topology')

    def _initialize_coordinates(self,inpfnm='init'):
        """Builds initial top and gro files for initial liquid simulation

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        """
        if cp.passed('initialize_coordinates'): return
        TC=self.TopoCoord
        cwd=pfs.go_to('systems/init')
        densification_dict=self.cfg.parameters.get('densification',{})
        # logger.debug(f'{densification_dict}')
        if not densification_dict:
            densification_dict['aspect_ratio']=self.cfg.parameters.get('aspect_ratio',np.array([1.,1.,1.]))
            if 'initial_density' in self.cfg.parameters:
                densification_dict['initial_density']=self.cfg.parameters['initial_density']
            if 'initial_boxsize' in self.cfg.parameters:
                densification_dict['initial_boxsize']=self.cfg.parameters['initial_boxsize']
        dspec=['initial_density' in densification_dict,'initial_boxsize' in densification_dict]
        # logger.debug(f'{dspec} {any(dspec)} {not all(dspec)}')
        assert any(dspec),'Neither "initial_boxsize" nor "initial_density" are specfied'
        assert not all(dspec),'Cannot specify both "initial_boxsize" and "initial_density"'
        if 'initial_boxsize' in densification_dict:
            boxsize=densification_dict['initial_boxsize']
        else:
            mass_kg=TC.total_mass(units='SI')
            V0_m3=mass_kg/densification_dict['initial_density']
            ar=densification_dict.get('aspect_ratio',np.array([1.,1.,1.]))
            assert ar[0]==1.,f'Error: parameter aspect_ratio must be a 3-element-list with first element 1'
            ar_yx=ar[1]*ar[2]
            L0_m=(V0_m3/ar_yx)**(1./3.)
            L0_nm=L0_m*1.e9
            boxsize=np.array([L0_nm,L0_nm,L0_nm])*ar
            logger.info(f'Initial density: {densification_dict["initial_density"]} kg/m^3')
            logger.info(f'Total mass: {mass_kg:.3e} kg')
            logger.info(f'Box aspect ratio: {" x ".join([str(x) for x in ar])}')
        logger.info(f'Initial box side lengths: {boxsize[0]:.3f} nm x {boxsize[1]:.3f} nm x {boxsize[2]:.3f} nm')

        clist=[cc for cc in self.cfg.initial_composition if 'count' in cc]
        c_togromacs={}
        for cc in clist:
            M=self.cfg.molecules[cc['molecule']]
            tc=cc['count']
            if len(M.conformers)>0:
                nconf=len(M.conformers)
                if tc<nconf:
                    c_togromacs.update({u:1 for u in random.sample(M.conformers,tc)})
                else:
                    q,r=divmod(tc,nconf)
                    logger.debug(f'divmod({nconf},{tc}) gives {q}, {r}')
                    c_togromacs.update({u:q for u in M.conformers})
                    k=0
                    while r:
                        c_togromacs[M.conformers[k]]+=1
                        k+=1
                        r-=1
                for gro in M.conformers:
                    pfs.checkout(f'molecules/parameterized/{gro}.gro',altpath=[pfs.subpath('molecules')])        
            else:
                ''' assuming racemic mixture of any stereoisomers '''
                total_isomers=len(M.stereoisomers)+1
                count_per_isomer=tc//total_isomers
                leftovers=tc%total_isomers
                c_togromacs[M.name]=count_per_isomer
                if leftovers>0:
                    c_togromacs[M.name]+=1
                    leftovers-=1
                pfs.checkout(f'molecules/parameterized/{M.name}.gro',altpath=[pfs.subpath('molecules')])
                for isomer in M.stereoisomers:
                    c_togromacs[isomer]=count_per_isomer
                    if leftovers>0:
                        c_togromacs[isomer]+=1
                        leftovers-=1
                    pfs.checkout(f'molecules/parameterized/{isomer}.gro',altpath=[pfs.subpath('molecules')])
        logger.debug(f'Sending to insert_molecules: {c_togromacs}')
        msg=insert_molecules(c_togromacs,boxsize,inpfnm,inputs_dir='../../molecules/parameterized',**self.cfg.parameters)
        TC.read_gro(f'{inpfnm}.gro')
        TC.atom_count()
        # TC.set_grx_attributes(['z','nreactions','reactantName','cycle','cycle_idx','chain','chain_idx'])
        TC.set_grx_attributes()
        TC.inherit_grx_attributes_from_molecules(self.cfg.molecules,self.cfg.initial_composition,
            globally_unique=[False,False,False,True,False,True,False,True,True],
            unset_defaults=[0,0,'UNSET',-1,-1,-1,-1,-1,'UNSET'])
        for list_name in ['cycle','chain']:
            TC.reset_idx_list_from_grx_attributes(list_name)
        chainlists=TC.idx_lists['chain']
        # logger.debug(f'virgin chains')
        # for i,c in enumerate(chainlists):
        #     logger.debug(f'  {i} {c}')
        TC.make_resid_graph()
        TC.write_grx_attributes(f'{inpfnm}.grx')
        logger.info(f'Coordinates "{inpfnm}.gro" in {pfs.cwd()}')
        logger.info(f'Extended attributes "{inpfnm}.grx" in {pfs.cwd()}')
        cp.set(TC,stepname='initialize_coordinates')

    def _do_equilibration(self,edict={},deffnm='equilibrate',plot_pfx=''):
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        TC=self.TopoCoord
        TC.equilibrate(deffnm,edict,gromacs_dict,plot_pfx)

    def _do_equilibration_series(self,eq_dict={},deffnm='equilibrate',plot_pfx=''):
        if not eq_dict: return
        for stg in eq_dict:
            self._do_equilibration(stg,deffnm,plot_pfx)
                
    def _do_pap(self,pfx=''):
        if not pfx: return
        assert pfx in ['precure','postcure']
        pap_dict=self.cfg.parameters.get(pfx,{})
        if not pap_dict: return
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        preequil=pap_dict.get('preequilibration',{})
        anneal=pap_dict.get('anneal',{})
        postequil=pap_dict.get('postequilibration',{})
        if not any([preequil,anneal,postequil]): return
        if not _nonempty_directives([x for x in [preequil,anneal,postequil] if x]): return
        cwd=pfs.go_to(f'systems/{pfx}')
        my_logger(f'{pfx.capitalize()} in {pfs.cwd()}',logger.info)
        if preequil: self._do_equilibration(preequil,deffnm='preequilibration',plot_pfx=f'{pfx}-preequilibration')
        if anneal and _nonempty_directives([anneal]): 
            self._do_anneal(anneal,deffnm='annealed')
            trace('Temperature',['annealed'],outfile=os.path.join(pfs.proj(),f'plots/{pfx}-anneal-T.png'))
        if postequil: self._do_equilibration(postequil,deffnm='postequilibration',plot_pfx=f'{pfx}-postequilibration')

    def _do_anneal(self,anneal_dict={},deffnm='anneal'):
        if not anneal_dict: return
        ncycles=anneal_dict.get('ncycles',0)
        if not ncycles: return
        mdp_pfx='nvt' # assume all annealing run at NVT
        TC=self.TopoCoord
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        pfs.checkout(f'mdp/{mdp_pfx}.mdp')
        timestep=float(mdp_get(f'{mdp_pfx}.mdp','dt'))
        cycle_segments=anneal_dict.get('cycle_segments',[])
        temps=[str(r['T']) for r in cycle_segments]
        durations=[r.get('ps',0) for r in cycle_segments]
        if not any(durations):
            durations=['{:.2f}'.format(r.get('nsteps',)*timestep) for r in cycle_segments]
            if not any(durations):
                raise Exception(f'No "ps" or "nsteps" found.')
        cycle_duration=sum(durations)
        total_duration=cycle_duration*ncycles
        nsteps=int(total_duration/timestep)
        cum_time=durations.copy()
        for i in range(1,len(cum_time)):
            cum_time[i]+=cum_time[i-1]
        logger.info(f'Annealing: {len(durations)} points for {ncycles} cycles over {total_duration} ps')
        mod_dict={
            'ref_t':anneal_dict.get('initial_temperature',300.0),
            'gen-temp':anneal_dict.get('initial_temperature',300.0),
            'gen-vel':'yes',
            'annealing-npoints':len(cycle_segments),
            'annealing-temp':' '.join(temps),
            'annealing-time':' '.join([f'{x:.2f}' for x in cum_time]),
            'annealing':'periodic' if ncycles>1 else 'single',
            'nsteps':nsteps
            }
        mdp_modify(f'{mdp_pfx}.mdp',mod_dict)
        msg=TC.grompp_and_mdrun(out=deffnm,mdp=mdp_pfx,quiet=False,**gromacs_dict)
        logger.info(f'Annealed coordinates in {deffnm}.gro')

