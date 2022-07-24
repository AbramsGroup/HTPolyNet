"""

read cfg file
@author: huang, abrams

"""
import json
import yaml
import os
import logging
import numpy as np
from copy import deepcopy
from itertools import product

from HTPolyNet.molecule import Molecule, Reaction

logger=logging.getLogger(__name__)

def _determine_sequence(m,moldict,atoms):
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
        ''' Dictionary of name:Molecule '''
        self.molecules = {}
        ''' List of Reaction instances '''
        self.reactions = []
        ''' all other parameters in cfg file '''
        self.parameters = {}

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
        if 'symmetry_equivalent_atoms' in self.parameters:
            if mol_name in self.parameters['symmetry_equivalent_atoms']:
                M.symmetry_relateds=self.parameters['symmetry_equivalent_atoms'][mol_name]
        if 'stereocenters' in self.parameters:
            if mol_name in self.parameters['stereocenters']:
                M.stereocenters=self.parameters['stereocenters'][mol_name]
                extra_stereocenters=[]
                for stc in M.stereocenters:
                    for sc in M.symmetry_relateds:
                        if stc in sc:
                            sc_copy=sc.copy()
                            sc_copy.remove(stc)
                            extra_stereocenters.extend(sc_copy)
                M.stereocenters.extend(extra_stereocenters)
        return M

    def parse(self):
        """parse parse self.basedict to set Title, initial_composition, and lists of
           reactions and molecules.
           
        """
        self.Title=self.basedict.get('Title','No Title Provided')
        self.parameters=self.basedict
        if not 'ncpu' in self.parameters:
            self.parameters['ncpu']=os.cpu_count()
        '''
        add a molecule for each unique molecule specified in initial_composition
        '''
        self.initial_composition=self.basedict.get('initial_composition',[])
        for item in self.initial_composition:
            m=item['molecule']
            if m not in self.molecules:
                logger.debug(f'new {m}')
                self.molecules[m]=self.NewMolecule(m)
        '''
        add additional molecules that appear as either reactants or products of reactions
        '''
        rlist=self.basedict.get('reactions',[])
        self.reactions=[Reaction(r) for r in rlist]
        for R in self.reactions:
            '''
            add every reactant in this reaction to the list of molecules
            '''
            for rnum,rname in R.reactants.items():
                zrecs=[]
                logger.debug(f'{R.name} rname {rname}')
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        # this atom is in this reactant
                        zrecs.append(cprec)
                if not rname in self.molecules:
                    self.molecules[rname]=self.NewMolecule(rname)
                ''' provide molecule with records of atoms that have z values '''
                self.molecules[rname].update_zrecs(zrecs)
        for R in self.reactions:
            '''
            add product of this reaction or update generator if product is already in the list of molecules
            '''
            r=R.product
            if not r in self.molecules:
                self.molecules[r]=self.NewMolecule(r)
            self.molecules[r].generator=R

    def get_reaction(self,product_name):
        if '-' in product_name:
            base_product_name,sym=product_name.split('-')
        else:
            base_product_name,sym=product_name,'0'
        for r in self.reactions:
            if r.product==base_product_name:
                return r
        return None

    def calculate_maximum_conversion(self):
        N={}
        for item in self.initial_composition:
            N[item['molecule']]=item['count']
        Bonds=[]
        Atoms=[]
        for R in self.reactions:
            if R.stage=='cure':
                for b in R.bonds:
                    A,B=b['atoms']
                    a,b=R.atoms[A],R.atoms[B]
                    aan,ban=a['atom'],b['atom']
                    ari,bri=a['resid'],b['resid']
                    arnum,brnum=a['reactant'],b['reactant']
                    arn,brn=R.reactants[arnum],R.reactants[brnum]
                    if arnum!=brnum:
                        az,bz=a['z'],b['z']
                        ia=(aan,ari,arnum,arn,az)
                        ib=(ban,bri,brnum,brn,bz)
                        b=(ia,ib)
                        if ia not in Atoms and arn in N:
                            Atoms.append(ia)
                        if ib not in Atoms and brn in N:
                            Atoms.append(ib)
                        if b not in Bonds and arn in N and brn in N:
                            Bonds.append(b)
        # logger.debug(f'atomset: {Atoms}')
        Z=[]
        for a in Atoms:
            Z.append(a[4]*N[a[3]])
        # logger.debug(f'Z: {Z}')
        # logger.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            a,b=B
            az=Z[Atoms.index(a)]
            bz=Z[Atoms.index(b)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(a)]-=MaxB[-1]
            Z[Atoms.index(b)]-=MaxB[-1]
        # logger.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        return sum(MaxB)
