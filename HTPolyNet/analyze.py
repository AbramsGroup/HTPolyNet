"""

.. module:: analyze
   :synopsis: handles the analyze subcommand
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import json
import yaml
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.gromacs import gmx_command
import HTPolyNet.software as software
from HTPolyNet.configuration import Configuration
from pathlib import Path

logger=logging.getLogger(__name__)

class Analyze:
    allowed_keys=['gromacs','command','subdir','options','links','outfile','console-input','matchlines']
    required_keys=['command','subdir']
    default_params={
        'gromacs' : {
            'gmx': 'gmx'
        },
    }
    def __init__(self,indict,strict=True):
        self.params={}
        for p,v in self.default_params.items():
            self.params[p]=indict.get(p,v)
        for p,v in indict.items():
            if not p in self.allowed_keys:
                if strict:
                    logger.info(f'Ignoring directive \'{p}\' in yaml input file')
            if p in self.default_params:
                logger.info(f'Overwriting default {p} value')
            self.params[p]=v
        self.console_output=None

    def do(self,**gromacs_dict):
        """do handles executing the analysis

        """
        p=self.params
        print(p)
        for rk in self.required_keys:
            assert rk in p, f'Error: no {rk} value found'
        # logger.info(f'do {p}')
        # if a gromacs dict is passed in, assume this overrides the one read in from the file
        if gromacs_dict:
            software.set_gmx_preferences(gromacs_dict)
        else:
            software.set_gmx_preferences(p['gromacs'])
        logger.info(f'going to {p["subdir"]}')
        pfs.go_to(p['subdir'])
        # make symlinks to requested files
        symlinks=p.get('links',[])
        for input_file in symlinks:
            srcnm=os.path.join(pfs.proj(),input_file)
            bsnm=os.path.basename(srcnm)
            chk=Path(bsnm)
            if not chk.is_symlink():
                os.symlink(srcnm,bsnm)
            else:
                logger.info(f'Symlink {bsnm} already exists.')
        cfile=''
        if 'console-input' in p:
            cfile='console-in.txt'
            ci=p['console-input']
            with open(cfile,'w') as f:
                for ch in ci:
                    f.write(ch+'\n')
        self.console_output=gmx_command(p['command'],p.get('options',{}),console_in=cfile)
        logger.info(f'Command {p["command"]} completed.')
    
    def parse_console_output(self):
        if not self.console_output:
            logger.info(f'No console output')
            return
        p=self.params
        # either we are grepping out lines from console output or putting it all out there
        if not 'outfile' in p:
            logger.info(f'Here is the console output')
            logger.info(self.console_output)
        else:
            if not 'matchlines' in p:
                with open(p['outfile'],'w') as f:
                    f.write(self.console_output)
            else:
                svlns=[]
                console_lines=self.console_output.split('\n')
                for cl in console_lines:
                    for ml in p['matchlines']:
                        if ml in cl:
                            svlns.append(cl)
                with open(p['outfile'],'w') as f:
                    for s in svlns:
                        f.write(s+'\n')
            logger.info(f'Created {p["outfile"]} in {p["subdir"]}')

class AnalyzeDensity(Analyze):
    """ Analyze class for handling trajectory density profile calculation
    """
    default_params={
        'subdir': 'analyze/density',
        'links': ['postsim/equilibrate/equilibrate.tpr','postsim/equilibrate/equilibrate.trr'],
        'gromacs' : {
            'gmx': 'gmx'
        },
        'command': 'density',
        'options': {
            's':'equilibrate.tpr',
            'f':'equilibrate.trr',
            'o':'density.xvg',
            'xvg': 'none',
            'b': 0,
            'd': 'Z',
            'sl': 50
        },
        'console-input': ['0']
    }

class AnalyzeFFV(Analyze):
    default_params={
        'subdir': 'analyze/freevolume',
        'links': ['postsim/equilibrate/equilibrate.tpr','postsim/equilibrate/equilibrate.trr'],
        'gromacs' : {
            'gmx': 'gmx'
        },
        'command': 'freevolume',
        'options': {
            's':'equilibrate.tpr',
            'f':'equilibrate.trr',
            'o':'ffv.xvg',
            'xvg': 'none',
            'b': 0.0
        },
        'outfile': 'ffv.dat',
        'matchlines': ['Free volume','Total volume','Number of molecules','Average molar mass','Density','Molecular volume Vm assuming homogeneity:','Molecular van der Waals volume assuming homogeneity:','Fractional free volume']
    }


class AnalyzeConfiguration:
    """ handles reading and parsing an analysis input config file.
        Config file format
        
        - { key1: {<paramdict>}}
        - { key2: {<paramdict>}} 
        
        ...

        The config file is a list of single-element dictionaries, whose single keyword
        indicates the type of analysis to be run; analyses are run in the order
        they appear in the config file.
        
        """
    default_class=Analyze
    predefined_classes={'density':AnalyzeDensity,'freevolume':AnalyzeFFV}
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
            content=p['analysis']
            assert len(p)==1,f'Poorly formatted {self.cfgFile}; each stanza may have only one keyword \'analysis\''
            analysistype=content['command']#list(p.keys())[0]
            if analysistype in self.predefined_classes:
                self.stagelist.append(self.predefined_classes[analysistype](content))
            else:
                self.stagelist.append(self.default_class(content))

def analyze(args):
    """postsim handles the analyze subcommand for managing gromacs-based trajectory analyses

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
    cfg=AnalyzeConfiguration.read(args.cfg)
    logger.debug(f'{cfg.baselist}')
    logger.info(f'Project director{ess}: {args.proj}')
    software.sw_setup()
    logger.debug(f'ogromacs {ogromacs}')
    for d in args.proj:
        pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots','postsim','analyze'],verbose=True,projdir=d,reProject=False,userlibrary=args.lib)
        pfs.go_to('analyze')
        for stage in cfg.stagelist:
            stage.do(**ogromacs)
            stage.parse_console_output()
        pfs.go_root()
