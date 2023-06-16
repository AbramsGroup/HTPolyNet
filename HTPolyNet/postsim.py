"""

.. module:: postsim
   :synopsis: handles the postsim subcommand
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import shutil
import numpy as np
import os
import json
import yaml
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.gromacs import mdp_get, mdp_modify, gmx_energy_trace
import HTPolyNet.software as software
from HTPolyNet.configuration import Configuration
from HTPolyNet.plot import scatter

logger=logging.getLogger(__name__)

class PostSimMD:
    """ Generic class for handling post-cure md simulations; this one just does simple NPT MD equilibration;
    Classes that inherit from this class should define their own default_params and build_npt
    """
    default_params={
        'subdir': 'postsim/equilibrate',
        'input_top': 'systems/final-results/final.top',
        'input_gro': 'systems/final-results/final.gro',
        'input_grx': 'systems/final-results/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'ps': 1000,
        'T': 300,
        'P':1,
        'output_deffnm': 'equilibrate',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png')
    }
    def __init__(self,indict,strict=True):
        self.params={}
        for p,v in self.default_params.items():
            self.params[p]=indict.get(p,v)
        for p,v in indict.items():
            if not p in self.default_params:
                if strict:
                    logger.info(f'Ignoring directive \'{p}\' in yaml input file')
                else:
                    self.params[p]=v
                    
    def do(self,mdp_pfx='npt',**gromacs_dict):
        """do handles executing the postsim MD simulation

        :param mdp_pfx: filename prefix for output files, defaults to 'npt'
        :type mdp_pfx: str, optional
        """
        p=self.params
        logger.info(f'do {p}')
        # if a gromacs dict is passed in, assume this overrides the one read in from the file
        if gromacs_dict:
            software.set_gmx_preferences(gromacs_dict)
        else:
            software.set_gmx_preferences(p['gromacs'])
        logger.info(f'going to {p["subdir"]}')
        pfs.go_to(p['subdir'])
        for input_file in ['input_top','input_gro','input_grx']:
            srcnm=os.path.join(pfs.proj(),p[input_file])
            shutil.copy(srcnm,'.')
        local_top=os.path.basename(p['input_top'])
        local_gro=os.path.basename(p['input_gro'])
        local_grx=os.path.basename(p['input_grx'])
        TC=TopoCoord(topfilename=local_top,grofilename=local_gro,grxfilename=local_grx)
        logger.info(f'{TC.Coordinates.A.shape[0]} atoms {TC.total_mass(units="gromacs"):.2f} amu')
        box=TC.Coordinates.box
        pfs.checkout('mdp/npt.mdp')
        os.rename('npt.mdp',f'{mdp_pfx}.mdp')
        self.build_mdp(f'{mdp_pfx}.mdp',box=box)
        msg=TC.grompp_and_mdrun(out=p['output_deffnm'],mdp=mdp_pfx,quiet=False,mylogger=logger.info,**gromacs_dict)
        df=gmx_energy_trace(p['output_deffnm'],p['traces'])
        use_scatter=p['scatter']
        temp_x=None
        temp_y=None
        scat_file=use_scatter[2]
        for sznm in ['Box-X','Box-Y','Box-Z']:
            if sznm in p['traces']:
                L0=box[0][0] if sznm=='Box-X' else box[1][1] if sznm=='Box-Y' else box[2][2]
                temp_x=f'{sznm}-strain'
                df[temp_x]=df[sznm]/L0-1.0
        for sznm in ['Pres-XX','Pres-YY','Pres-ZZ']:
            if sznm in p['traces']:
                df[f'{sznm}-stress']=df[sznm]*(-1)
                temp_y=[f'{sznm}-stress']
        if temp_x and temp_y:
            use_scatter=(temp_x,temp_y,scat_file)
        df.to_csv(f'{p["output_deffnm"]}.csv',header=True,index=False)
        scatter(df,*use_scatter)
        logger.info(f'Final coordinates in {p["output_deffnm"]}.gro')
        logger.info(f'Traces saved in {p["output_deffnm"]}.csv')

    def build_mdp(self,mdpname,**kwargs):
        """build_mdp builds the GROMACS mdp file required for an NPT equilibration

        :param mdpname: name of mdp file
        :type mdpname: str
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        duration=params['ps']
        nsteps=int(duration/timestep)
        mod_dict={
            'ref_t':params['T'],
            'ref_p':params['P'],
            'gen-temp':params['T'],
            'gen-vel':'yes',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)
        
class PostSimAnneal(PostSimMD):
    """ a class to handle temperature annealing MD simulation 
    """
    default_params={
        'subdir':'postsim/anneal',
        'input_top': 'systems/final-results/final.top',
        'input_gro': 'systems/final-results/final.gro',
        'input_grx': 'systems/final-results/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun'
        },
        'output_deffnm':'anneal',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png'),
        'T0': 300,
        'T1': 600,
        'ncycles': 1,
        'T0_to_T1_ps': 1000,
        'T1_ps': 1000,
        'T1_to_T0_ps': 1000,
        'T0_ps': 1000,
        'P':1
    }
    def build_mdp(self,mdpname,**kwargs):
        """build_mdp builds the GROMACS mdp file required for an annealing MD simulation

        :param mdpname: name of mdp file
        :type mdpname: str
        """        
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        timeints=[0.0,params['T0_to_T1_ps'],params['T1_ps'],params['T1_to_T0_ps'],params['T0_ps']]
        timepoints=[0.0]
        temppoints=[params['T0'],params['T1'],params['T1'],params['T0'],params['T0']]
        for i in range(1,len(timeints)):
            timepoints.append(timepoints[-1]+timeints[i])
        duration=timepoints[-1]
        nsteps=int(duration/timestep)*params['ncycles']
        mod_dict={
            'ref_t':params['T0'],
            'ref_p':params['P'],
            'gen-temp':params['T0'],
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'periodic' if params['ncycles']>1 else 'single',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)

class PostSimLadder(PostSimMD):
    """ a class to handle a temperature-ladder MD simulation
    """
    default_params={
        'subdir':'postsim/ladder',
        'input_top': 'systems/final-results/final.top',
        'input_gro': 'systems/final-results/final.gro',
        'input_grx': 'systems/final-results/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'output_deffnm':'ladder',
        'traces': ['Temperature','Density','Volume'],
        'scatter': ('time(ps)',['Density'],'rho_v_ns.png'),
        'Tlo':300.0,
        'Thi':600.0,
        'deltaT':5,
        'ps_per_run':1000,
        'ps_per_rise':1000,
        'warmup_ps':5000,
        'P':1
    }

    def build_mdp(self,mdpname,**kwargs):
        """build_mdp builds the GROMACS mdp file required for a temperature-ladder MD simulation

        :param mdpname: name of mdp file
        :type mdpname: str
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        Tdomain=params['Thi']-params['Tlo']
        nSteps=int(Tdomain/np.abs(params['deltaT']))+1
        T0=params['Tlo'] if params['deltaT']>0.0 else params['Thi']
        T1=params['Thi'] if params['deltaT']>0.0 else params['Tlo']
        Tladder=np.linspace(T0,T1,nSteps)
        timepoints=[0.0,params['warmup_ps']]
        temppoints=[T0,T0]
        for i in range(nSteps):
            cT=Tladder[i]
            timepoints.append(timepoints[-1]+params['ps_per_rise'])
            temppoints.append(cT)
            timepoints.append(timepoints[-1]+params['ps_per_run'])
            temppoints.append(cT)
        duration=timepoints[-1]
        nMDsteps=int(duration/timestep)
        mod_dict={
            'ref_t':T0,
            'ref_p':params['P'],
            'gen-temp':T0,
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'single',
            'nsteps':nMDsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(mdpname,mod_dict)

class PostSimDeform(PostSimMD):
    """ a class to handle a uniaxial deformation MD simulation
    """
    default_params={
        'subdir':'postsim/deform-x',
        'input_top': 'systems/final-results/final.top',
        'input_gro': 'postsim/equilibrate/equilibrate.gro',
        'input_grx': 'systems/final-results/final.grx',
        'gromacs' : {
            'gmx': 'gmx',
            'mdrun': 'gmx mdrun',
            'options': '-quiet -nobackup',
            'mdrun_single_molecule': 'gmx mdrun mdrun'
        },
        'output_deffnm':'deform-x',
        'traces': ['Box-X','Pres-XX'],
        'scatter': ('Box-X',['Pres-XX'],'tension_v_xlength.png'),
        'direction':'x',
        'T':300.0,
        'P':1.0,
        'ps':1000,
        'edot': 0.001 # strain rate in ps^-1
    }

    def build_mdp(self,mdpname,**kwargs):
        """build_mdp builds the GROMACS mdp file required for a temperature-ladder MD simulation

        :param mdpname: name of mdp file
        :type mdpname: str
        """
        params=self.params
        timestep=float(mdp_get(mdpname,'dt'))
        duration=params['ps']
        nsteps=int(duration/timestep)
        box=kwargs.get('box',np.array([[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]))
        edot=params.get('edot',0.0)
        direction=params.get('direction','x')
        mod_dict={
            'ref_t':params['T'],
            'ref_p':params['P'],
            'gen-temp':params['T'],
            'gen-vel':'yes',
            'tcoupl':'v-rescale',
            'nsteps': nsteps,
            'rlist': 1.2,
            'rcoulomb': 1.2,
            'rvdw': 1.2,
            'tau_t':0.5,
            'tau_p':1.0,
            'refcoord_scaling': 'com',
            'pcoupltype': 'anisotropic'
            }
        if direction=='x':
            strain_vel=box[0][0]*edot
            mod_dict['ref_p']='0.0 1.0 1.0 0 0 0'
            mod_dict['compressibility']='0.0 4.5e-5 4.5e-5 0 0 0'
            mod_dict['deform']=f'{strain_vel:.3e} 0 0 0 0 0'
            params['output_deffnm'] = 'deform-x'
            params['traces']=['Box-X','Pres-XX']
            params['scatter']=('Box-X',['Pres-XX'],'tension_v_xlength.png')
        elif direction=='y':
            strain_vel=box[1][1]*edot
            mod_dict['ref_p']='1.0 0.0 1.0 0 0 0'
            mod_dict['compressibility']='4.5e-5 0.0 4.5e-5 0 0 0'
            mod_dict['deform']=f'0 {strain_vel:.3e} 0 0 0 0'
            params['output_deffnm'] = 'deform-y'
            params['traces']=['Box-Y','Pres-YY']
            params['scatter']=('Box-Y',['Pres-YY'],'tension_v_ylength.png')
        elif direction=='z':
            strain_vel=box[2][2]*edot
            mod_dict['ref_p']='1.0 1.0 0.0 0 0 0'
            mod_dict['compressibility']='4.5e-5 4.5e-5 0.0 0 0 0'
            mod_dict['deform']=f'0 0 {strain_vel:.3e} 0 0 0'
            params['output_deffnm'] = 'deform-z'
            params['traces']=['Box-Z','Pres-ZZ']
            params['scatter']=('Box-Z',['Pres-ZZ'],'tension_v_zlength.png')
        else:
            logger.error(f'Bad direction for uniaxial strain {direction}')

        mdp_modify(mdpname,mod_dict)

class PostsimConfiguration:
    """ handles reading and parsing a postcure simulation input config file.
        Config file format
        
        - { key1: {<paramdict>}}
        - { key2: {<paramdict>}} 
        
        ...

        The config file is a list of single-element dictionaries, whose single keyword
        indicates the type of MD simulation to be run; simulations are run in the order
        they appear in the config file.

        Currently allowed simulation types:

        - 'equilibrate': simple NPT equilibration;
        - 'anneal': simple simulated annealing;
        - 'ladder': temperature ladder;
        - 'deform: constant strain-rate uniaxial deformation;
        
        """
    default_classes={'equilibrate':PostSimMD,'anneal':PostSimAnneal,'ladder':PostSimLadder,'deform':PostSimDeform}
    def __init__(self):
        self.cfgFile=''
        self.baselist=[]
        self.stagelist=[]

    @classmethod
    def read(cls,filename,parse=True,**kwargs):
        """read generates a new PostsimConfiguration object by reading in the JSON or YAML file indicated by filename

        :param filename: name of file from which to read new PostsimConfiguration object
        :type filename: str
        :param parse: if True, parse the input configuration file, defaults to True
        :type parse: bool, optional
        :raises Exception: if extension of filename is not '.json' or '.yaml' or '.yml'
        :return: a new PostsimConfiguration object
        :rtype: PostsimConfiguration
        """
        basename,extension=os.path.splitext(filename)
        if extension=='.json':
            return cls._read_json(filename,parse,**kwargs)
        elif extension=='.yaml' or extension=='.yml':
            return cls._read_yaml(filename,parse,**kwargs)
        else:
            raise Exception(f'Unknown config file extension {extension}')

    @classmethod
    def _read_json(cls,filename,parse=True,**kwargs):
        """_read_json create a new PostsimConfiguration object by reading from JSON input

        :param filename: name of JSON file
        :type filename: str
        :param parse: if True, parse the JSON data, defaults to True
        :type parse: bool, optional
        :return: a new PostsimConfiguration object
        :rtype: PostsimConfiguration
        """
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.baselist=json.load(f)
            assert type(inst.baselist)==list,f'Poorly formatted {filename}'
        if parse: inst.parse(**kwargs)
        return inst

    @classmethod
    def _read_yaml(cls,filename,parse=True,**kwargs):
        """_read_yaml create a new PostsimConfiguration object by reading from YAML input

        :param filename: name of YAML file
        :type filename: str
        :param parse: if True, parse the YAML data, defaults to True
        :type parse: bool, optional
        :return: a new PostsimConfiguration object
        :rtype: PostsimConfiguration
        """
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.baselist=yaml.safe_load(f)
            assert type(inst.baselist)==list,f'Poorly formatted {filename}'
        if parse: inst.parse(**kwargs)
        return inst

    def parse(self,**kwargs):
        """parse parses a PostsimConfiguration file to build the list of stages to run
        """
        for p in self.baselist:
            assert len(p)==1,f'Poorly formatted {self.cfgFile}; each stanza may have only one keyword'
            simtype=list(p.keys())[0]
            assert simtype in self.default_classes,f'Simulation type "{simtype}" in {self.cfgFile} not understood.'
            logger.info(f'passing in {p[simtype]}')
            self.stagelist.append(self.default_classes[simtype](p[simtype]))

def postsim(args):
    """postsim handles the postsim subcommand for managing post-cure production MD simulations

    :param args: command-line arguments
    :type args: argparse.Namespace
    """
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(format='%(levelname)s> %(message)s',level=loglevel_numeric)
    ess='y' if len(args.proj)==0 else 'ies'
    ogromacs={}
    if args.ocfg:
        ocfg=Configuration.read(args.ocfg,parse=False)
        ogromacs=ocfg.basedict.get('gromacs',{})
    cfg=PostsimConfiguration.read(args.cfg)
    logger.debug(f'{cfg.baselist}')
    logger.info(f'Project director{ess}: {args.proj}')
    software.sw_setup()
    logger.debug(f'ogromacs {ogromacs}')
    for d in args.proj:
        pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots','postsim'],verbose=True,projdir=d,reProject=False,userlibrary=args.lib)
        pfs.go_to('postsim')
        for stage in cfg.stagelist:
            stage.do(mdp_pfx='local',**ogromacs)
        pfs.go_root()

