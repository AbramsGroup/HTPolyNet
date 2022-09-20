import logging
import os
import shutil
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
from HTPolyNet.molecule import Molecule, MoleculeDict
from HTPolyNet.reaction import is_reactant
from HTPolyNet.expandreactions import chain_expand_reactions
from HTPolyNet.reaction import reaction_stage
from HTPolyNet.curecontroller import CureController, CureState
from HTPolyNet.stringthings import my_logger

logger=logging.getLogger(__name__)

def logrotate(filename):
    """logrotate renames any existing log file with name 'filename' using a pattern that
    prevents overwriting any old logs

    :param filename: name of log file
    :type filename: str
    """
    if os.path.exists(filename):
        n=1
        while os.path.exists(f'#{n}#{filename}'):
            n+=1
        shutil.copyfile(filename,f'#{n}#{filename}')

_directives=['ps','nsteps','ncycles']
def _nonempty_directives(dirlist):
    """_nonempty_directives checks the list of directives for existence of at least
    one of the strings in the _directives list

    :param dirlist: list of dictionaries
    :type dirlist: lsit
    :return: True if one of the strings in _directives is present in one of the dictionaries
    :rtype: Boolean
    """
    amts=[]
    for d in dirlist:
        ta=[d.get(x,0) for x in _directives]
        amts.append(any(ta))
    return any(amts)

class Runtime:
    """ Runtime class manages all aspects of a system build using HTPolyNet.
    """
    default_edict={ 'ensemble': 'npt', 'temperature': 300, 'pressure': 10, 'ps': 200, 'nsteps': -2, 'repeat': 0 }
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

    def __init__(self,cfgfile='',restart=False):
        my_logger(software.to_string(),logger.info)
        self.cfgfile=cfgfile
        if cfgfile=='':
            logger.error('HTPolyNet requires a configuration file.')
            raise RuntimeError('HTPolyNet requires a configuration file')
        logger.info(f'Configuration: {cfgfile}')
        self.cfg=Configuration.read(os.path.join(pfs.root(),cfgfile))
        """ Fill in any default values """
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
            logger.debug('Setting up cure controller')
            self.cc=CureController(cure_dict)
        self.ncpu=self.cfg.parameters.get('ncpu',os.cpu_count())

    def generate_molecules(self,force_parameterization=False,force_checkin=False):
        """generate_molecules manages creation and parameterization of all monomers
        and oligomer templates

        :param force_parameterization: forces AmberTools to run parameterizations, defaults to False
        :type force_parameterization: bool, optional
        :param force_checkin: forces HTPolyNet to overwrite molecule library, defaults to False
        :type force_checkin: bool, optional
        """
        GAFF_dict=self.cfg.parameters.get('GAFF',{})
        self.molecules={}
        ''' configuration.parse() generated a list of Molecules implied by configuration; assume
            they are all unparameterized '''
        for mname,M in self.cfg.molecules.items():
            M.set_origin('unparameterized')
        pfs.go_to('molecules/parameterized')
        my_logger(f'Templates in {pfs.cwd()}',logger.info)
        ''' Each molecule implied by the cfg is 'generated' here, either by
            reading from the library or direct parameterization.  In some cases,
            the molecule is to be generated by a reaction; if so, it's
            `generator` attribute will be a Reaction instance '''
        ess='' if len(self.cfg.molecules)==1 else 's'
        logger.info(f'{len(self.cfg.molecules)} molecule{ess} detected in {self.cfgfile}')
        replns=[f'{msg:>30s}: {count:<5d}' for msg,count in self.cfg.molecule_report.items()]
        for ln in replns:
            logger.info(ln)
        ml=list(self.cfg.molecules.keys())
        logger.debug(f'Generating: {ml}')
        for mname,M in self.cfg.molecules.items():
            self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
            self.molecules[mname]=M
        ''' Generate any required template products that result from reactions in which the bond generated creates
            dihedrals that span more than just the two monomers that are connected '''
        new_reactions,new_molecules=chain_expand_reactions(self.molecules)
        if len(new_molecules)>0:
            ess='' if len(new_molecules)==1 else 's'
            logger.info(f'{len(new_molecules)} molecule{ess} implied by chaining')
            ml=list(new_molecules.keys())
            # logger.info(ml)
            self.cfg.reactions.extend(new_reactions)
            make_molecules={k:v for k,v in new_molecules.items() if k not in self.molecules}
            for mname,M in make_molecules.items():
                # logger.debug(f'Generating {mname}:')
                self._generate_molecule(M,force_parameterization=force_parameterization,force_checkin=force_checkin)
                assert M.get_origin()!='unparameterized'
                self.molecules[mname]=M
                logger.debug(f'Generated {mname}')

        for M in self.molecules:
            self.molecules[M].is_reactant=is_reactant(M,self.cfg.reactions,stage=reaction_stage.cure)

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
            cmp_msg=", ".join([(x["molecule"]+" "+str(x["count"])) for x in self.cfg.initial_composition if x["count"]>0])
            logger.info(f'Initial composition is {cmp_msg}')
            self.cfg.calculate_maximum_conversion()
            logger.info(f'100% conversion is {self.cfg.maxconv} bonds')

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

    @cp.enableCheckpoint
    def do_initialization(self,inpfnm='init'):
        """do_initialization manages creation of the global system topology file
        and the execution of 'gmx insert-molecules' to create the initial
        simulation box

        :param inpfnm: basename for output files, defaults to 'init'
        :type inpfnm: str, optional
        :return: dictionary of Gromacs files
        :rtype: dict
        """
        my_logger(f'Initialization in {pfs.cwd()}',logger.info)
        top=self._initialize_topology(inpfnm)
        gro,grx=self._initialize_coordinates(inpfnm)
        return {'top':top,'gro':gro,'grx':grx}

    @cp.enableCheckpoint
    def do_densification(self,deffnm='densified'):
        """do_densification manages execution of gmx mdrun to perform minimization
            and NPT MD simulation of the initial liquid system.  Final coordinates are
            loaded into the global TopoCoord.

        :param inpfnm: input file name prefix, defaults to 'init'
        :type inpfnm: str, optional
        :param deffnm: deffnm prefix fed to gmx mdrun, defaults to 'npt-1'
        :type deffnm: str, optional
        :return: dictionary of Gromacs file names
        :rtype: dict
        """
        my_logger(f'Densification in {pfs.cwd()}',logger.info)
        densification_dict=self.cfg.parameters.get('densification',{})
        assert len(densification_dict)>0,'"densification" directives missing'
        equilibration=densification_dict.get('equilibration',[])
        assert len(equilibration)>0,'equilibration directives missing'
        TC=self.TopoCoord
        infiles=[TC.files[x] for x in ['gro','top','grx']]
        assert all([os.path.exists(x) for x in infiles]),f'One or more of {infiles} not found'
        self._do_equilibration_series(equilibration,deffnm=f'{deffnm}',plot_pfx='densification')
        logger.info(f'Densified coordinates in {pfs.cwd()}/{os.path.basename(TC.files["gro"])}')
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def do_precure(self):
        """do_precure manages execution of mdrun to perform any pre-cure MD simulation(s).

        :return: dictionary of Gromacs file names
        :rtype: dict
        """
        self._do_pap(pfx='precure')
        return {c:os.path.basename(x) for c,x in self.TopoCoord.files.items() if c!='mol2'}

    def do_cure(self):
        """do_cure manages CURE algorithm execution.

        :return: only returns None if CURE cannot execute
        :rtype: None
        """
        if not hasattr(self,'cc'): 
            logger.debug(f'no cure controller')
            return  # no cure controller
        cc=self.cc
        TC=self.TopoCoord
        RL=self.cfg.reactions
        MD=self.molecules
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        pfs.go_proj()
        ''' read cure state or initialize new cure '''
        if os.path.exists('systems/cure_state.yaml'):
            cc.state=CureState.from_yaml('systems/cure_state.yaml')
            my_logger(f'Connect-Update-Relax-Equilibrate (CURE) resumes',logger.info)
            my_logger(f'at iteration {cc.state.iter} and {cc.state.cum_nxlinkbonds} bonds',logger.info)
        else:
            cc.setup(max_nxlinkbonds=self.cfg.maxconv,desired_nxlinkbonds=int(self.cfg.maxconv*cc.dicts['controls']['desired_conversion']),max_search_radius=float(min(TC.Coordinates.box.diagonal()/2)))
            cc.state.iter=1
            my_logger('Connect-Update-Relax-Equilibrate (CURE) begins',logger.info)
        cure_finished=cc.is_cured()
        if cure_finished: 
            logger.debug('cure finished even before loop')
            return
        ''' perform CURE iterations '''
        logger.info(f'Attempting to form {cc.state.desired_nxlinkbonds} bonds')
        while not cure_finished:
            pfs.go_to(f'systems/iter-{cc.state.iter}')
            cc.do_iter(TC,RL,MD,gromacs_dict=gromacs_dict)
            cure_finished=cc.is_cured()
            if not cure_finished:
                cure_finished=cc.next_iter()
        ''' perform capping if necessary '''
        my_logger(f'Capping begins',logger.info)
        pfs.go_to(f'systems/capping')
        cc.do_capping(TC,RL,MD,gromacs_dict=gromacs_dict)
        my_logger('Connect-Update-Relax-Equilibrate (CURE) ends',logger.info)

    @cp.enableCheckpoint
    def do_postcure(self):
        """do_postcure manages execution of mdrun to perform any post-cure MD simulation(s).

        :return: dictionary of Gromacs file names
        :rtype: dict
        """
        self._do_pap(pfx='postcure')
        return {c:os.path.basename(x) for c,x in self.TopoCoord.files.items() if c!='mol2'}

    @cp.enableCheckpoint
    def save_data(self,result_name='final'):
        """save_data writes 'gro', 'top', and 'grx' files for system

        :param result_name: output file basename, defaults to 'final'
        :type result_name: str, optional
        :return: dictionary of Gromacs file basenames
        :rtype: dict
        """
        TC=self.TopoCoord
        my_logger(f'Final data to {pfs.cwd()}',logger.info)
        TC.write_grx_attributes(f'{result_name}.grx')
        TC.write_gro(f'{result_name}.gro')
        TC.write_top(f'{result_name}.top')
        return {c:os.path.basename(x) for c,x in TC.files.items() if c!='mol2'}

    def do_workflow(self,**kwargs):
        """do_workflow manages runtime for one entire system-build workflow
        """
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        TC=self.TopoCoord
        pfs.go_proj()
        self.generate_molecules(
            force_parameterization=force_parameterization,  # force antechamber/GAFF parameterization
            force_checkin=force_checkin                     # force check-in to system libraries
        )
        pfs.go_proj()
        last_data=cp.read_checkpoint()
        logger.debug(f'Checkpoint last_data {last_data}')
        TC.load_files(last_data)
        pfs.go_to(f'systems/init')
        self.do_initialization()
        pfs.go_to(f'systems/densification')
        self.do_densification()
        pfs.go_to(f'systems/precure')
        self.do_precure()
        self.do_cure()
        pfs.go_to(f'systems/postcure')
        self.do_postcure()
        pfs.go_to(f'systems/final-results')
        self.save_data()
        pfs.go_proj()

    def _generate_molecule(self,M:Molecule,**kwargs):
        """_generate_molecule generates and parameterizes a single molecule based on the partially 
        complete instance in parameter M

        :param M: partially-complete molecule instance (populated by config reader)
        :type M: Molecule
        """
        if M.origin!='unparameterized' and M.origin!='symmetry_product': return
        mname=M.name
        checkin=pfs.checkin
        force_parameterization=kwargs.get('force_parameterization',False)
        force_checkin=kwargs.get('force_checkin',False)
        ''' Either perform a fresh parameterization or fetch parameterization files '''
        if force_parameterization or not M.previously_parameterized():
            logger.debug(f'Parameterization of {mname} requested -- can we generate {mname}?')
            generatable=(not M.generator) or (all([m in self.molecules for m in M.generator.reactants.values()]))
            if generatable:
                logger.debug(f'Generating {mname}')
                M.generate(available_molecules=self.molecules,**self.cfg.parameters)
                for ex in ['mol2','top','itp','gro','grx']:
                    checkin(f'molecules/parameterized/{mname}.{ex}',overwrite=force_checkin)
                M.set_origin('newly parameterized')
            else:
                logger.debug(f'Error: could not generate {mname}')
                logger.debug(f'not ({mname}.generator) = {bool(not M.generator)}')
                if M.generator:
                    logger.debug(f'reactants {list(M.generator.reactants.values())}')
                return
        else:
            logger.debug(f'Fetching parameterized {mname}')
            exts=pfs.fetch_molecule_files(mname)
            logger.debug(f'fetched {mname} exts {exts}')
            M.load_top_gro(f'{mname}.top',f'{mname}.gro',mol2filename='',wrap_coords=False)
            M.TopoCoord.read_gro_attributes(f'{mname}.grx')
            M.set_sequence_from_coordinates()
            if M.generator:
                M.prepare_new_bonds(available_molecules=self.molecules)
            M.set_origin('previously parameterized')

        ''' Generate any stereoisomers and/or conformers '''
        M.generate_stereoisomers()
        M.generate_conformers()

    def _type_consistency_check(self,typename='dihedraltypes',funcidx=4,selection_rule='stiffest'):
        """_type_consistency_check checks current topology for any inconsistencies and corrects them

        :param typename: type of Gromacs interaction to check, defaults to 'dihedraltypes'
        :type typename: str, optional
        :param funcidx: function type of interaction, defaults to 4
        :type funcidx: int, optional
        :param selection_rule: one-word designation for correction mechanism, defaults to 'stiffest'
        :type selection_rule: str, optional
        """
        logger.debug(f'Consistency check of {typename} func {funcidx} on all {len(self.molecules)} molecules requested')
        mnames=list(self.molecules.keys())
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
        """_initialize_topology creates a full gromacs topology that includes all directives necessary
            for an initial liquid simulation.  This will NOT use any #include's;
            all types will be explicitly in-lined.

        :param inpfnm: input file basename, defaults to 'init'
        :type inpfnm: str, optional
        :return: name of top file
        :rtype: str
        """
        # if cp.passed('initialize_topology'): return
        if os.path.isfile(f'{inpfnm}.top'):
            logger.debug(f'{inpfnm}.top already exists in {pfs.cwd()} but we will rebuild it anyway!')

        ''' for each monomer named in the cfg, either parameterize it or fetch its parameterization '''
        TC=self.TopoCoord
        already_merged=[]
        for item in self.cfg.initial_composition:
            M=self.molecules[item['molecule']]
            N=item['count']
            if not N: continue
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
        return f'{inpfnm}.top'

    def _initialize_coordinates(self,inpfnm='init'):
        """_initialize_coordinates builds initial top and gro files for initial liquid simulation

        :param inpfnm: input file basename, defaults to 'init'
        :type inpfnm: str, optional
        :return: 'gro' and 'grx' file names
        :rtype: tuple
        """
        densification_dict=self.cfg.parameters.get('densification',{})
        # logger.debug(f'{densification_dict}')
        TC=self.TopoCoord
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
            if not tc: continue
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
        TC.set_grx_attributes()
        TC.inherit_grx_attributes_from_molecules(self.cfg.molecules,self.cfg.initial_composition)
        for list_name in ['cycle','chain']:
            TC.reset_idx_list_from_grx_attributes(list_name)
        TC.make_resid_graph()
        TC.write_grx_attributes(f'{inpfnm}.grx')
        logger.info(f'Coordinates "{inpfnm}.gro" in {pfs.cwd()}')
        logger.info(f'Extended attributes "{inpfnm}.grx" in {pfs.cwd()}')
        return f'{inpfnm}.gro',f'{inpfnm}.grx'

    def _do_equilibration(self,edict={},deffnm='equilibrate',plot_pfx=''):
        """_do_equilibration manages execution of TopoCoord.equilibrate using
        the directive passed in parameter 'edict'

        :param edict: equilibration directives, defaults to {}
        :type edict: dict, optional
        :param deffnm: mdrun default file basename, defaults to 'equilibrate'
        :type deffnm: str, optional
        :param plot_pfx: basename for any output plots, defaults to ''
        :type plot_pfx: str, optional
        """
        gromacs_dict=self.cfg.parameters.get('gromacs',{})
        TC=self.TopoCoord
        for k in edict.keys():
            if not k in self.default_edict:
                logger.debug(f'ignoring unknown equilibration directive "{k}"')
        edr_list=TC.equilibrate(deffnm,edict,gromacs_dict)
        ens=edict['ensemble']
        if ens=='npt' and plot_pfx!='':
            trace('Density',edr_list,outfile=os.path.join(pfs.proj(),f'plots/{plot_pfx}-density.png'))

    def _do_equilibration_series(self,eq_stgs=[],deffnm='equilibrate',plot_pfx=''):
        """_do_equilibration_series manages a series of equilibrations

        :param eq_stgs: list of stage designations, defaults to []
        :type eq_stgs: list, optional
        :param deffnm: mdrun default file basename, defaults to 'equilibrate'
        :type deffnm: str, optional
        :param plot_pfx: basename for any output plots, defaults to ''
        :type plot_pfx: str, optional
        """
        if not eq_stgs: return
        for stg in eq_stgs:
            self._do_equilibration(stg,deffnm,plot_pfx)
                
    def _do_pap(self,pfx='precure'):
        """_do_pap manages a preanneal-anneal-postanneal series of MD simulations

        :param pfx: either 'precure' or 'postcure', defaults to 'precure'
        :type pfx: str, optional
        """
        if not pfx: return
        assert pfx in ['precure','postcure']
        pap_dict=self.cfg.parameters.get(pfx,{})
        if not pap_dict: return
        preequil=pap_dict.get('preequilibration',{})
        anneal=pap_dict.get('anneal',{})
        postequil=pap_dict.get('postequilibration',{})
        if not any([preequil,anneal,postequil]): return
        if not _nonempty_directives([x for x in [preequil,anneal,postequil] if x]): return
        my_logger(f'{pfx.capitalize()} in {pfs.cwd()}',logger.info)
        if preequil: self._do_equilibration(preequil,deffnm='preequilibration',plot_pfx=f'{pfx}-preequilibration')
        if anneal and _nonempty_directives([anneal]): 
            self._do_anneal(anneal,deffnm='annealed')
            trace('Temperature',['annealed'],outfile=os.path.join(pfs.proj(),f'plots/{pfx}-anneal-T.png'))
        if postequil: self._do_equilibration(postequil,deffnm='postequilibration',plot_pfx=f'{pfx}-postequilibration')

    def _do_anneal(self,anneal_dict={},deffnm='anneal'):
        """_do_anneal manages execution of an annealing MD simulation

        :param anneal_dict: annealing directives, defaults to {}
        :type anneal_dict: dict, optional
        :param deffnm: mdrun default file basename, defaults to 'anneal'
        :type deffnm: str, optional
        :raises Exception: If no run duration is indicated in the any of the annealing segments
        """
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

