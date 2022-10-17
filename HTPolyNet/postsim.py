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

class Tanneal:
    """ a class to handle temperature annealing MD simulation 
    """
    default_params={
        'subdir':'postsim/anneal',
        'input_deffnm':'systems/final-results/final',
        'output_deffnm':'anneal',
        'T0': 300,
        'T1': 600,
        'ncycles': 1,
        'T0_to_T1_ps': 1000,
        'T1_ps': 1000,
        'T1_to_T0_ps': 1000,
        'T0_ps': 1000,
        'gromacs': {'gmx':'gmx'}
    }
    def __init__(self,indict):
        self.params={}
        for p,v in self.default_params.items():
            self.params[p]=indict.get(p,v)
    def __str__(self):
        p=self.params
        restr=f'Tanneal: {p["T0"]} to {p["T1"]}'
        return restr
    
    def do(self,mdp_pfx='npt',**gromacs_dict):
        """do handles executing the temperature-anneal on the passed-in system

        :param pfs: the project file system
        :type pfs: ProjectFileSystem
        :param mdp_pfx: filename prefix for output files, defaults to 'npt'
        :type mdp_pfx: str, optional
        """
        p=self.params
        software.set_gmx_preferences(p.get('gromacs',gromacs_dict))
        logger.info(f'going to {p["subdir"]}')
        pfs.go_to(p['subdir'])
        for sfx in ['.gro','.top','.grx']:
            srcnm=os.path.join(pfs.proj(),p['input_deffnm']+sfx)
            shutil.copy(srcnm,'.')
            logger.info(f'Copied {srcnm} to ./')
        local_deffnm=os.path.basename(p['input_deffnm'])
        TC=TopoCoord(topfilename=f'{local_deffnm}.top',grofilename=f'{local_deffnm}.gro',grxfilename=f'{local_deffnm}.grx')
        logger.info(f'{TC.Coordinates.A.shape[0]} atoms {TC.total_mass(units="gromacs"):.2f} amu')
        pfs.checkout('mdp/npt.mdp')
        os.rename('npt.mdp',f'{mdp_pfx}.mdp')
        timestep=float(mdp_get(f'{mdp_pfx}.mdp','dt'))
        timeints=[0.0,p['T0_to_T1_ps'],p['T1_ps'],p['T1_to_T0_ps'],p['T0_ps']]
        timepoints=[0.0]
        temppoints=[p['T0'],p['T1'],p['T1'],p['T0'],p['T0']]
        for i in range(1,len(timeints)):
            timepoints.append(timepoints[-1]+timeints[i])
        duration=timepoints[-1]
        nsteps=int(duration/timestep)*p['ncycles']
        mod_dict={
            'ref_t':p['T0'],
            'gen-temp':p['T0'],
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'periodic' if p['ncycles']>1 else 'single',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(f'{mdp_pfx}.mdp',mod_dict)
        msg=TC.grompp_and_mdrun(out=p['output_deffnm'],mdp=mdp_pfx,quiet=False,mylogger=logger.info,**gromacs_dict)
        df=gmx_energy_trace(p['output_deffnm'],['Temperature','Density','Volume'])
        scatter(df,'Temperature',['Density'],'rho_v_T.png')
        df.to_csv(f'{p["output_deffnm"]}.csv',header=True,index=False)
        logger.info(f'Final coordinates in {p["output_deffnm"]}.gro')
        logger.info(f'Traces saved in {p["output_deffnm"]}.csv')

class Tladder:
    """ a class to handle a temperature-ladder MD simulation
    """
    default_params={
        'subdir':'postsim/ladder',
        'input_deffnm':'systems/final-results/final',
        'output_deffnm':'ladder',
        'Tlo':300.0,
        'Thi':600.0,
        'Ntemps':31,
        'ps_per_run':1000,
        'ps_per_rise':1000,
        'warmup_ps':5000,
        'gromacs': {'gmx':'gmx'}
    }
    def __init__(self,indict):
        self.params={}
        for p,v in self.default_params.items():
            self.params[p]=indict.get(p,v)

    def __str__(self):
        p=self.params
        restr=f'Tladder: {p["Tlo"]}K -> {p["Tlo"]} K in {p["Ntemps"]} rungs at {p["ps_per_run"]} ps per run and {p["ps_per_rise"]} ps per rise; warmup for {p["warmup_ps"]} ps'
        return restr

    def do(self,mdp_pfx='npt',**gromacs_dict):
        """do handles executing the temperature-ladder

        :param mdp_pfx: filename prefix for output files, defaults to 'npt'
        :type mdp_pfx: str, optional
        """
        p=self.params
        software.set_gmx_preferences(p.get('gromacs',gromacs_dict))
        pfs.go_to(p['subdir'])
        for sfx in ['.gro','.top','.grx']:
            srcnm=os.path.join(pfs.proj(),p['input_deffnm']+sfx)
            shutil.copy(srcnm,'.')
        local_deffnm=os.path.basename(p['input_deffnm'])
        TC=TopoCoord(topfilename=f'{local_deffnm}.top',grofilename=f'{local_deffnm}.gro',grxfilename=f'{local_deffnm}.grx')
        logger.info(f'{TC.Coordinates.A.shape[0]} atoms {TC.total_mass(units="gromacs"):.2f} amu')
        pfs.checkout(f'mdp/npt.mdp')
        os.rename('npt.mdp',f'{mdp_pfx}.mdp')
        timestep=float(mdp_get(f'{mdp_pfx}.mdp','dt'))
        Tladder=np.linspace(p['Tlo'],p['Thi'],p['Ntemps'])
        timepoints=[0.0,p['warmup_ps']]
        temppoints=[p['Tlo'],p['Tlo']]
        for i in range(p['Ntemps']):
            cT=Tladder[i]
            timepoints.append(timepoints[-1]+p['ps_per_rise'])
            temppoints.append(cT)
            timepoints.append(timepoints[-1]+p['ps_per_run'])
            temppoints.append(cT)
        duration=timepoints[-1]
        nsteps=int(duration/timestep)
        mod_dict={
            'ref_t':p['Tlo'],
            'gen-temp':p['Thi'],
            'gen-vel':'yes',
            'annealing-npoints':len(timepoints),
            'annealing-temp':' '.join([str(x) for x in temppoints]),
            'annealing-time':' '.join([str(x) for x in timepoints]),
            'annealing':'single',
            'nsteps':nsteps,
            'tcoupl':'v-rescale','tau_t':0.5
            }
        mdp_modify(f'{mdp_pfx}.mdp',mod_dict)
        msg=TC.grompp_and_mdrun(out=p['output_deffnm'],mdp=mdp_pfx,quiet=False,mylogger=logger.info,**gromacs_dict)
        df=gmx_energy_trace(p['output_deffnm'],['Temperature','Density','Volume'])
        scatter(df,'Temperature',['Density'],'rho_v_T.png')
        df.to_csv(f'{p["output_deffnm"]}.csv',header=True,index=False)
        logger.info(f'Final coordinates in {p["output_deffnm"]}.gro')
        logger.info(f'Traces saved in {p["output_deffnm"]}.csv')

class PostsimConfiguration:
    """ handles reading and parsing a postcure simulation input config file.
        Config file format:
        
        - { key1: {<paramdict>}}
        - { key2: {<paramdict>}} 
        ...

        The config file is a list of single-element dictionaries, whose single keyword
        indicates the type of MD simulation to be run; simulations are run in the order
        they appear in the config file.

        Currently allowed simulation types:

        - 'anneal': simple simulated annealing;
        - 'ladder': temperature ladder
        - 'deform': deform along one axis (not yet implemented)
        
        """
    default_classes={'anneal':Tanneal,'ladder':Tladder,'deform':None}
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
        """_read_json create a new Configuration object by reading from JSON input

        :param filename: name of JSON file
        :type filename: str
        :param parse: if True, parse the JSON data, defaults to True
        :type parse: bool, optional
        :return: a new Configuration object
        :rtype: Configuration
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
        """_read_yaml create a new Configuration object by reading from YAML input

        :param filename: name of YAML file
        :type filename: str
        :param parse: if True, parse the YAML data, defaults to True
        :type parse: bool, optional
        :return: a new Configuration object
        :rtype: Configuration
        """
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.baselist=yaml.safe_load(f)
            assert type(inst.baselist)==list,f'Poorly formatted {filename}'
        if parse: inst.parse(**kwargs)
        return inst

    def parse(self,**kwargs):
        for p in self.baselist:
            assert len(p)==1,f'Poorly formatted {self.cfgFile}; each stanza may have only one keyword'
            simtype=list(p.keys())[0]
            assert simtype in self.default_classes,f'Simulation type "{simtype}" in {self.cfgFile} not understood.'
            self.stagelist.append(self.default_classes[simtype](p[simtype]))


def postsim(args):
    """postsim handles the postsim subcommand for managing post-cure production MD simulations

    :param args: command-line arguments
    :type args: argparse.Namespace
    """
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(format='%(levelname)s> %(message)s',level=loglevel_numeric)
    ess='y' if len(args.proj)==0 else 'ies'
    ocfg=Configuration.read(args.ocfg,parse=False)
    cfg=PostsimConfiguration.read(args.cfg)
    logger.info(f'{cfg.baselist}')
    logger.info(f'Project director{ess}: {args.proj}')
    software.sw_setup()
    ogromacs=ocfg.basedict.get('gromacs',{})
    for d in args.proj:
        pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots','postsim'],verbose=True,projdir=d,reProject=False,userlibrary=args.lib)
        pfs.go_to('postsim')
        for stage in cfg.stagelist:
            stage.do(mdp_pfx='local',**ogromacs)

