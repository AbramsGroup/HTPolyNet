from imp import init_builtin
import json
import yaml
import os
from copy import deepcopy
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
        extension=filename.split('.')[-1]
        if extension=='json':
            return cls._read_json(filename,parse,**kwargs)
        elif extension=='yaml' or extension=='yml':
            return cls._read_yaml(filename,parse,**kwargs)
        else:
            raise Exception(f'Unknown config file extension {extension}')

    @classmethod
    def _read_json(cls,filename,parse=True,**kwargs):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=json.load(f)
        if parse: inst.parse(**kwargs)
        return inst

    @classmethod
    def _read_yaml(cls,filename,parse=True,**kwargs):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=yaml.safe_load(f)
        if parse: inst.parse(**kwargs)
        return inst
    
    def NewMolecule(self,mol_name,molrec={}):
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
                # logger.debug(f'z: update {rname} num {rnum} in rxn {R.name}')
                zrecs=[]
                for atnum,atrec in R.atoms.items():
                    # logger.debug(f'query {atrec} for "reactant" {rnum}')
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        zrecs.append(cprec)
                self.molecules[rname].update_zrecs(zrecs,self.molecules)
        # self.reactions.extend(new_reactions)

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
                # TODO: fix this -- z's defined by residue not reactant!!
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
            # Z.append(a[4]*N[a[3]])
        logger.debug(f'Z: {Z}')
        logger.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            # a,b=B
            az=Z[Atoms.index(B.ai)]
            bz=Z[Atoms.index(B.aj)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(B.ai)]-=MaxB[-1]
            Z[Atoms.index(B.aj)]-=MaxB[-1]
        logger.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        self.maxconv=sum(MaxB)
        # return sum(MaxB)
