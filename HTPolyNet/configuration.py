import json
import yaml
import os
from copy import deepcopy
import logging
from collections import namedtuple
from HTPolyNet.molecule import Molecule, MoleculeDict, generate_stereo_reactions
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

        self.initial_composition = []
        ''' raw dict read from JSON/YAML '''
        self.basedict = {}
        self.maxconv=0.0

    @classmethod
    def read(cls,filename):
        extension=filename.split('.')[-1]
        if extension=='json':
            return cls.read_json(filename)
        elif extension=='yaml' or extension=='yml':
            return cls.read_yaml(filename)
        else:
            raise Exception(f'Unknown config file extension {extension}')

    @classmethod
    def read_json(cls,filename):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=json.load(f)
        inst.parse()
        return inst

    @classmethod
    def read_yaml(cls,filename):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=yaml.safe_load(f)
        inst.parse()
        return inst
    
    def NewMolecule(self,mol_name,molrec={}):
        return Molecule.New(mol_name,molrec)
        
    def parse(self):
        """parse self.basedict to set Title, initial_composition, and lists of
           reactions and molecules.
        """
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
        mol_reac_detected=extract_molecule_reactions(self.reactions)
        for mname,gen in mol_reac_detected:
            self.molecules[mname]=Molecule.New(mname,gen,self.constituents.get(mname,{}))
            for si,S in self.molecules[mname].stereoisomers.items():
                self.molecules[si]=S
        for mname,M in self.molecules.items():
            M.set_sequence_from_moldict(self.molecules)
            logging.debug(f'{mname} seq: {M.sequence}')
        for R in self.reactions:
            for rnum,rname in R.reactants.items():
                zrecs=[]
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        zrecs.append(cprec)
                self.molecules[rname].update_zrecs(zrecs,self.molecules)

        generate_stereo_reactions(self.reactions,self.molecules)
        # self.reactions.extend(new_reactions)

        for r in self.reactions:
            logger.debug(f'{str(r)}')

        for m,M in self.molecules.items():
            R=M.generator
            if not R==None:
                logger.debug(f'{m}: {R.name}')
            else:
                logger.debug(f'{m}: None')
            logger.debug(f'zrecs: {M.zrecs}')
        
        self.initial_composition=[]
        for molecule,mrec in self.constituents.items():
            self.initial_composition.append({'molecule':molecule,'count':mrec.get('count',0)})


    def calculate_maximum_conversion(self):
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
        Bonds=[]
        Atoms=[]
        for R in [x for x in self.reactions if x.stage==reaction_stage.cure]:
            for b in R.bonds:
                A,B=b['atoms']
                a,b=R.atoms[A],R.atoms[B]
                aan,ban=a['atom'],b['atom']
                ari,bri=a['resid'],b['resid']
                arnum,brnum=a['reactant'],b['reactant']
                arn,brn=R.reactants[arnum],R.reactants[brnum]
                if arnum==brnum:  continue # this is an intermolecular reaction
                az,bz=a['z'],b['z']
                ia=Atom(aan,ari,arnum,arn,az)
                ib=Atom(ban,bri,brnum,brn,bz)
                b=Bond(ia,ib)
                if ia not in Atoms and arn in N:
                    Atoms.append(ia)
                if ib not in Atoms and brn in N:
                    Atoms.append(ib)
                if b not in Bonds and arn in N and brn in N:
                    Bonds.append(b)
        # logger.debug(f'atomset: {Atoms}')
        Z=[]
        for a in Atoms:
            Z.append(a.z*N[a.reactantName])
            # Z.append(a[4]*N[a[3]])
        # logger.debug(f'Z: {Z}')
        # logger.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            # a,b=B
            az=Z[Atoms.index(B.ai)]
            bz=Z[Atoms.index(B.aj)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(B.ai)]-=MaxB[-1]
            Z[Atoms.index(B.aj)]-=MaxB[-1]
        # logger.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        self.maxconv=sum(MaxB)
        # return sum(MaxB)
