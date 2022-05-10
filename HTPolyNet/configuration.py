# -*- coding: utf-8 -*-
"""

read cfg file
@author: huang, abrams

"""
import json
from threading import currentThread
import yaml
import os
import logging
import numpy as np
from copy import deepcopy
from itertools import product

from HTPolyNet.molecule import Molecule, Reaction

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
        ''' List of molecules for which symmetry-equivalence is requested '''
        self.use_sea=[]
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
    
    def parse(self):
        privileged_keys=['reactions','Title','use_sea','initial_composition']
        self.Title=self.basedict.get('Title','No Title Provided')
        ''' reactions must declare molecules '''
        self.reactions=[Reaction(r) for r in self.basedict['reactions']]
        for R in self.reactions:
            for rkey,r in R.reactants.items():
                ''' reactants do not get assigned generators if they are *only* reactants '''
                if not r in self.molecules:
                    self.molecules[r]=Molecule(r)
                else:
                    logging.debug(f'Reactant {r} in reaction {R.name} is already on the global Molecules list')
                    if self.molecules[r].generator:
                        logging.debug(f'{r} is a product of reaction {self.molecules[r].generator.name}')
            if not R.product in self.molecules:
                self.molecules[R.product]=Molecule(R.product,generator=R)
        self.initial_composition=self.basedict.get('initial_composition',[])
        ''' any molecules lists in the initial composition
        may not have been declared in reactions '''
        for item in self.initial_composition:
            m=item['molecule']
            if m not in self.molecules:
                self.molecules[m]=Molecule(m)

        for mname,M in self.molecules.items():
            M.sequence,M.reactive_atoms_seq=_determine_sequence(mname,self.molecules,[])
            logging.debug(f'Sequence of {mname}: {M.sequence}')
            logging.debug(f'Reactive atoms per sequence: {M.reactive_atoms_seq}')

        self.use_sea=self.basedict.get('use_sea',[])
        for m in self.use_sea:
            if not m in self.molecules:
                logging.error(f'Configuration {self.cfgFile} references undeclared molecule {m} in use_sea')
                logging.error(f'Molecules must be declared in reactions or initial_composition')
                raise Exception('Configuration error')

        self.privileged_items={}
        for k in privileged_keys:
            self.privileged_items[k]=self.basedict[k]
            del self.basedict[k]

        self.parameters=self.basedict
        if not 'cpu' in self.parameters:
            self.parameters['cpu']=os.cpu_count()
        return self

    def symmetry_expand_reactions(self,unique_molecules):
        logging.debug('symmetry_expand_reactions')
        extra_reactions=[]
        current_molecules=unique_molecules
        extra_molecules={}
        trydict={}
        sclass={}
        for R in self.reactions:
            logging.debug(f'Reaction {R.name}')
            pro_seq,pro_ra=_determine_sequence(R.product,current_molecules,[])
            logging.debug(f'{R.product} {pro_seq} {pro_ra}')
            sseq=[]
            for mname,matrectlist in zip(pro_seq,pro_ra):
                # how many se versions of this reactant?
                # -> product of number of members of symmetry class of each reactive atom in this
                # reactant
                residue=current_molecules[mname]
                sclass[mname]={}
                sp=[]
                for atrec in matrectlist:
                    atomName=atrec['atom']
                    clu=residue.atoms_w_same_attribute_as(find_dict={'atomName':atomName},    
                                                same_attribute='sea-idx',
                                                return_attribute='atomName')
                    logging.debug(f'atoms symmetry-equivalent to {mname} {atomName}: {clu}')
                    sclass[mname][atomName]=list(clu)
                    sp.append(list(clu))
                logging.debug(f'sp {sp}')
                if len(R.reactants)>1:
                    sseq.append(list(product(*sp)))
                else:
                    sseq.append(list(zip(*sp)))
            logging.debug(f'sclass {sclass}')
            logging.debug(f'Reaction reactant sealist: {sseq}')                      
            P=list(product(*[v for v in sseq]))
            logging.debug(f'P {P}')
            idx=1
            trydict[P[0]]=R.product
            for p in P[1:]:
                logging.debug(f'{R.name}-{idx} {p}({len(p)}) -> {R.product}-{idx}')
                trydict[p]=f'{R.product}-{idx}'
                newR=deepcopy(R)
                newR.name+=f'-{idx}'
                newR.product+=f'-{idx}'
                ip=0
                pdiv={}
                for nri,nR in R.reactants.items():
                    nres_nR=len(current_molecules[nR].sequence)
                    subp=p[ip:ip+nres_nR]
                    ip+=nres_nR
                    nresname=nR
                    if subp in trydict:
                        nresname=trydict[subp]
                    logging.debug(f'  {nR} {nresname} {subp}')
                    newR.reactants[nri]=nresname
                    pdiv[nri]=subp
                for aL,atomrec in R.atoms.items():
                    atomName=atomrec['atom']
                    z=atomrec['z']
                    resid=atomrec['resid']
                    reactant_idx=atomrec['reactant']
                    reactant=R.reactants[reactant_idx]
                    resname=current_molecules[reactant].sequence[resid-1]
                    logging.debug(f'need to alter atom: {atomName} resid: {resid} reactant: {reactant_idx} ({reactant})')
                    logging.debug(f'have {pdiv[reactant_idx]} {len(pdiv[reactant_idx])}')
                    new_atomName=atomName
                    for q in pdiv[reactant_idx]:
                        logging.debug(f'{q}')
                        for qq in q:
                            logging.debug(f'is it {qq}?')
                            isit=qq in sclass[resname][atomName]
                            logging.debug(f'{isit}')
                            if isit:
                                new_atomName=qq
                    newR.atoms[aL]={'atom':new_atomName,'resid':resid,'reactant':reactant_idx,'z':z}
                    # pdiv[] here has one element for each reactive atom in reactant, but
                extra_reactions.append(newR)
                newP=Molecule(name=newR.product,generator=newR)
                extra_molecules[newR.product]=newP
                idx+=1
        logging.debug(f'trydict {trydict}')
        for R in self.reactions:
            logging.debug(R)
        logging.debug(f'extra reactions:')
        for R in extra_reactions:
            logging.debug(R)
        return extra_molecules

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
        logging.debug(f'CMC: extracting atoms from {len(self.reactions)} reactions')
        for R in self.reactions:
            for b in R.bonds:
                A,B=b['atoms']
                a,b=R.atoms[A],R.atoms[B]
                aan,ban=a['atom'],b['atom']
                ari,bri=a['resid'],b['resid']
                arn,brn=R.reactants[a['reactant']],R.reactants[b['reactant']]
                if arn!=brn:
                    az,bz=a['z'],b['z']
                    ia=(aan,ari,arn,az)
                    ib=(ban,bri,brn,bz)
                    b=(ia,ib)
                    if ia not in Atoms and arn in N:
                        Atoms.append(ia)
                    if ib not in Atoms and brn in N:
                        Atoms.append(ib)
                    if b not in Bonds and arn in N and brn in N:
                        Bonds.append(b)
        # logging.debug(f'atomset: {Atoms}')
        Z=[]
        for a in Atoms:
            Z.append(a[3]*N[a[2]])
        # logging.debug(f'Z: {Z}')
        # logging.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            a,b=B
            az=Z[Atoms.index(a)]
            bz=Z[Atoms.index(b)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(a)]-=MaxB[-1]
            Z[Atoms.index(b)]-=MaxB[-1]
        # logging.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        return sum(MaxB)
