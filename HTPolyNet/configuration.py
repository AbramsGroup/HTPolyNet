import json
import yaml
import os
from copy import deepcopy
import logging
from collections import namedtuple
from HTPolyNet.molecule import Molecule, MoleculeDict, Reaction, ReactionList

logger=logging.getLogger(__name__)

def _determine_sequence(m,moldict:MoleculeDict,atoms):
    if not moldict[m].generator:
        return [m],[atoms]
    thisseq=[]
    thisatm=[]
    for rid,mname in moldict[m].generator.reactants.items():
        atoms=[a for a in moldict[m].generator.atoms.values() if a['reactant']==rid]
        newseq,newatm=_determine_sequence(mname,moldict,atoms)
        thisatm.extend(newatm)
        thisseq.extend(newseq)
    return thisseq,thisatm

class Configuration:
    def __init__(self):
        self.cfgFile = ''
        self.Title = ''
        ''' List of (Molecule, count) '''
        self.initial_composition = []
        self.constituents = {}
        ''' Dictionary of name:Molecule '''
        self.molecules:MoleculeDict = {}
        ''' List of Reaction instances '''
        self.reactions:ReactionList = []
        ''' all other parameters in cfg file '''
        self.parameters = {}
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
    
    def NewMolecule(self,mol_name):
        M=Molecule(mol_name)
        if not mol_name in self.constituents: return M 
        molrec=self.constituents[mol_name]
        symmetry_relateds=molrec.get('symmetry_equivalent_atoms',[])
        if not symmetry_relateds:
            sea=self.parameters.get('symmetry_equivalent_atoms',{})
            symmetry_relateds=sea.get(mol_name,[])
        M.symmetry_relateds=symmetry_relateds

        stereocenters=molrec.get('stereocenters',[])
        if not stereocenters:
            sc=self.parameters.get('stereocenters',{})
            stereocenters=sc.get(mol_name,[])
        M.stereocenters=stereocenters
        extra_stereocenters=[]
        for stc in M.stereocenters:
            for sc in M.symmetry_relateds:
                if stc in sc:
                    sc_copy=sc.copy()
                    sc_copy.remove(stc)
                    if not sc_copy in extra_stereocenters:
                        extra_stereocenters.extend(sc_copy)
        M.stereocenters.extend(extra_stereocenters)
        M.conformers_dict=molrec.get('conformers',{})
        # M.conformer_generator=molrec.get('conformer_generator','gromacs')
        return M

    def parse(self):
        """parse self.basedict to set Title, initial_composition, and lists of
           reactions and molecules.
        """
        self.Title=self.basedict.get('Title','No Title Provided')
        self.parameters=self.basedict
        if not 'ncpu' in self.parameters:
            self.parameters['ncpu']=os.cpu_count()
        '''
        add a molecule for each unique molecule specified in initial_composition
        '''
        self.constituents=self.basedict.get('constituents',{})
        self.initial_composition=self.basedict.get('initial_composition',[])
        # assert self.constituents or self.initial_composition,'Must specify initial composition in config file'
        if not self.constituents:
            for crec in self.initial_composition:
                self.constituents[crec['molecule']]['count']=crec['count']
        if not self.initial_composition:
            for molecule,mrec in self.constituents.items():
                self.initial_composition.append({'molecule':molecule,'count':mrec['count']})
        # logger.debug(f'{self.constituents}')
        # logger.debug(f'{self.initial_composition}')
        for item in self.initial_composition:
            m=item['molecule']
            if m not in self.molecules:
                # logger.debug(f'new {m}')
                self.molecules[m]=self.NewMolecule(m)
        '''
        add additional molecules that appear as either reactants or products of reactions
        '''
        rlist=self.basedict.get('reactions',[])
        self.reactions=[Reaction(r) for r in rlist]
        '''
        check for an processive reactions
        '''
        processives=[]
        for i,R in enumerate(self.reactions):
            if R.procession:
                nReactions=[]
                logger.debug(f'Reaction {R.name} is processive: {R.procession}')
                final_product_name=R.product
                R.product=f'{final_product_name}_I0'
                original_name=R.name
                R.name=f'{R.name}-i0'
                base_reactant_key=R.procession['increment_resid']
                for c in range(R.procession['count']):
                    nR=deepcopy(R)
                    nR.name=f'{original_name}-i{c+1}'
                    nR.product=f'{final_product_name}_I{c+1}' if c<(R.procession['count']-1) else final_product_name
                    nR.reactants[base_reactant_key]=f'{final_product_name}_I{c}'
                    for ak,arec in nR.atoms.items():
                        if arec['reactant']==base_reactant_key: arec['resid']+=(c+1)
                    nR.procession={}
                    nReactions.append(nR)
                processives.append((i,nReactions))
        lasti=0
        for rec in processives:
            atidx,sublist=rec
            atidx+=lasti
            self.reactions=self.reactions[0:atidx+1]+sublist+self.reactions[atidx+1:]
            lasti+=len(sublist)

        for r in self.reactions:
            logger.debug(f'{str(r)}')

        for R in self.reactions:
            '''
            add product of this reaction or update generator if product is already in the list of molecules
            '''
            r=R.product
            if not r in self.molecules:
                self.molecules[r]=self.NewMolecule(r)
            self.molecules[r].generator=R

        for R in self.reactions:
            '''
            add every reactant in this reaction to the list of molecules
            '''
            for rnum,rname in R.reactants.items():
                zrecs=[]
                # logger.debug(f'{R.name} rname {rname}')
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        # this atom is in this reactant
                        zrecs.append(cprec)
                # logger.debug(f'zrecs {zrecs}')
                if not rname in self.molecules:
                    self.molecules[rname]=self.NewMolecule(rname)
                ''' provide molecule with records of atoms that have z values '''
                self.molecules[rname].update_zrecs(zrecs)

    def calculate_maximum_conversion(self):
        Atom=namedtuple('Atom',['name','resid','reactantKey','reactantName','z'])
        Bond=namedtuple('Bond',['ai','aj'])
        N={}
        for item in self.initial_composition:
            N[item['molecule']]=item['count']
        Bonds=[]
        Atoms=[]
        for R in [x for x in self.reactions if x.stage=='cure']:
            for b in R.bonds:
                A,B=b['atoms']
                a,b=R.atoms[A],R.atoms[B]
                aan,ban=a['atom'],b['atom']
                ari,bri=a['resid'],b['resid']
                arnum,brnum=a['reactant'],b['reactant']
                arn,brn=R.reactants[arnum],R.reactants[brnum]
                if arnum!=brnum:  # this is an intermolecular reaction
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
