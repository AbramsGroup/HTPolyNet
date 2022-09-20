"""

.. module:: configuration
   :synopsis: Manages reading and parsing of YAML configuration files
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import json
import yaml
import os
import logging
from collections import namedtuple
from HTPolyNet.molecule import Molecule, MoleculeDict, generate_stereo_reactions, generate_symmetry_reactions
from HTPolyNet.reaction import Reaction, ReactionList, parse_reaction_list, extract_molecule_reactions, reaction_stage

logger=logging.getLogger(__name__)

class Configuration:
    default_directives = {
        'Title': 'No title provided',
        'constituents': {},
        'reactions':[],
        'densification':{},
        'precure':{},
        'postcure':{},
        'CURE':{},
        'ambertools':{},
        'gromacs':{}
    }
    def __init__(self):
        self.cfgFile = ''
        self.Title = ''
        ''' List of (Molecule, count) '''
        self.constituents = {}
        ''' Dictionary of name:Molecule '''
        self.molecules:MoleculeDict = {}
        ''' List of Reaction instances '''
        self.reactions:ReactionList = []
        ''' all other parameters in cfg file '''
        self.parameters = {}

        self.molecule_report={}

        self.initial_composition = []
        ''' raw dict read from JSON/YAML '''
        self.basedict = {}
        self.maxconv=0.0

    @classmethod
    def read(cls,filename,parse=True,**kwargs):
        """read generates a new Configuration object by reading in the JSON or YAML file indicated by filename

        :param filename: name of file from which to read new Configuration object
        :type filename: str
        :param parse: if True, parse the input configuration file, defaults to True
        :type parse: bool, optional
        :raises Exception: if extension of filename is not '.json' or '.yaml' or '.yml'
        :return: a new Configuration object
        :rtype: Configuration
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
            inst.basedict=json.load(f)
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
            inst.basedict=yaml.safe_load(f)
        if parse: inst.parse(**kwargs)
        return inst
    
    def NewMolecule(self,mol_name,molrec={}):
        """NewMolecule generate and return a new Molecule object with name mol_name and populated via directives in molrec

        :param mol_name: name of new molecule
        :type mol_name: str
        :param molrec: dictionary of new molecule directives, defaults to {}
        :type molrec: dict, optional
        :return: the new Molecule object
        :rtype: Molecule
        """
        return Molecule.New(mol_name,molrec)
        
    def parse(self,**kwargs):
        """parse self.basedict to set Title, initial_composition, and lists of
           reactions and molecules.
        """
        plot_reaction_network=kwargs.get('plot_reaction_network',True)
        assert not self.basedict=={},f'Reading error for config'
        for d in self.basedict.keys():
            if not d in self.default_directives:
                logging.debug(f'Ignoring unknown directive "{d}" in {self.cfgFile}')

        self.Title=self.basedict.get('Title',self.default_directives['Title'])
        self.parameters=self.basedict
        if not 'ncpu' in self.parameters:
            self.parameters['ncpu']=os.cpu_count()
    
        self.constituents=self.basedict.get('constituents',{})
        rlist=self.basedict.get('reactions',[])

        base_reaction_list=[Reaction(r) for r in rlist]
        self.reactions=parse_reaction_list(base_reaction_list)
        mol_reac_detected=extract_molecule_reactions(self.reactions,plot=plot_reaction_network)
        self.molecule_report['explicit']=len(mol_reac_detected)
        self.molecule_report['implied by stereochemistry']=0
        self.molecule_report['implied by symmetry']=0
        for mname,gen in mol_reac_detected:
            self.molecules[mname]=Molecule.New(mname,gen,self.constituents.get(mname,{}))
            for si,S in self.molecules[mname].stereoisomers.items():
                self.molecules[si]=S
                self.molecule_report['implied by stereochemistry']+=1
        
        for mname,M in self.molecules.items():
            M.set_sequence_from_moldict(self.molecules)
            logging.debug(f'{mname} seq: {M.sequence}')

        self.molecule_report['implied by stereochemistry']+=generate_stereo_reactions(self.reactions,self.molecules)
        self.molecule_report['implied by symmetry']+=generate_symmetry_reactions(self.reactions,self.molecules)
        for R in self.reactions:
            for rnum,rname in R.reactants.items():
                zrecs=[]
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        zrecs.append(cprec)
                self.molecules[rname].update_zrecs(zrecs,self.molecules)

        for r in self.reactions:
            logger.debug(f'{str(r)}')

        for m,M in self.molecules.items():
            R=M.generator
            if not R==None:
                logger.debug(f'{m}: generator: {R.name}')
            else:
                logger.debug(f'{m}: generator: None')
            logger.debug(f'{m}: zrecs: {M.zrecs}')
        
        self.initial_composition=[]
        for molecule,mrec in self.constituents.items():
            self.initial_composition.append({'molecule':molecule,'count':mrec.get('count',0)})
        
        for molec in self.initial_composition:
            mname=molec['molecule']
            if not mname in self.molecules:
                self.molecules[mname]=Molecule.New(mname,None,self.constituents.get(mname,{}))
                self.molecules[mname].set_sequence_from_moldict(self.molecules)
                logging.debug(f'{mname} seq: {self.molecules[mname].sequence}')


    def calculate_maximum_conversion(self):
        """calculate_maximum_conversion calculates the maximum number of polymerization bonds that can form based on specified system composition and reactions
        """
        Atom=namedtuple('Atom',['name','resid','reactantKey','reactantName','z'])
        Bond=namedtuple('Bond',['ai','aj'])
        N={}
        # may have composite molecules
        for item in self.initial_composition:
            molecule_name=item['molecule']
            molecule_count=item.get('count',0)
            if molecule_count:
                for res in self.molecules[molecule_name].sequence:
                    if not res in N:
                        N[res]=0
                    N[res]+=molecule_count
        logger.debug(f'Base residue counts: {N}')
        Bonds=[]
        Atoms=[]
        for R in [x for x in self.reactions if x.stage==reaction_stage.cure]:
            logger.debug(str(R))
            for b in R.bonds:
                A,B=b['atoms']
                a,b=R.atoms[A],R.atoms[B]
                aan,ban=a['atom'],b['atom']
                ari,bri=a['resid'],b['resid']
                arnum,brnum=a['reactant'],b['reactant']
                arn,brn=R.reactants[arnum],R.reactants[brnum]
                if arnum==brnum:  continue # this is an intramolecular reaction
                az,bz=a['z'],b['z']
                ia=Atom(aan,ari,arnum,arn,az)
                ib=Atom(ban,bri,brnum,brn,bz)
                logger.debug(f'ia {ia} ib {ib} arn {arn} brn {brn}')
                b=Bond(ia,ib)
                if ia not in Atoms and arn in N:
                    Atoms.append(ia)
                if ib not in Atoms and brn in N:
                    Atoms.append(ib)
                if b not in Bonds and arn in N and brn in N:
                    Bonds.append(b)
        logger.debug(f'atomset: {Atoms}')
        Z=[]
        for a in Atoms:
            Z.append(a.z*N[a.reactantName])
        logger.debug(f'Z: {Z}')
        logger.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            az=Z[Atoms.index(B.ai)]
            bz=Z[Atoms.index(B.aj)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(B.ai)]-=MaxB[-1]
            Z[Atoms.index(B.aj)]-=MaxB[-1]
        logger.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        self.maxconv=sum(MaxB)
