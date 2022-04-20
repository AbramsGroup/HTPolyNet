from itertools import combinations_with_replacement, product
import os
from copy import deepcopy
import pandas as pd
import numpy as np
import logging

from HTPolyNet.topology import Topology
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.ambertools import GAFFParameterize
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.gromacs import grompp_and_mdrun, analyze_sea

# class Atom:
#     def __init__(self,datadict):
#         self.monomer=datadict['monomer']
#         self.name=datadict['atom']
#         # z: maximum number of bonds this atom will form
#         # in a crosslinking reaction
#         self.z=int(datadict.get('z',0))
#     def __str__(self):
#         return f'{self.monomer}_{self.name}({self.z})'
#     def to_yaml(self):
#         return r'{'+f'monomer: {self.monomer}, atom: {self.name}, z: {self.z}'+r'}'

class Reaction:
    ''' reactions may only be initialized after monomers '''
    def __init__(self,jsondict):
        self.jsondict=jsondict
        ''' attributes expected in cfg file '''
        self.name=jsondict.get('name','unnamed reaction')
        self.atoms=jsondict.get('atoms','empty reaction')
        self.bonds=jsondict.get('bonds','empty reaction')
        self.reactants=jsondict.get('reactants','empty reaction')
        self.product=jsondict.get('product','empty reaction')
        self.restrictions=jsondict.get('restrictions',{})
        self.stage=jsondict.get('stage','cure')

    def __str__(self):
        return f'Reaction "{self.name}"'

class Molecule:
    def __init__(self,name='',generator=None):
        self.name=name
        self.Topology=None
        self.Coords={}
        self.generator=generator
        self.cstale=''  # ['','gro','mol2']

    def __str__(self):
        restr=f'{self.name} '
        if not self.Topology:
            restr+=f'(empty topology) '
        if len(self.Coords)==0:
            restr+=f'(empty coords) '
        return restr+'\n'

    def parameterize(self,outname='',**kwargs):
        pfs.checkout(f'{self.name}.mol2')
        assert os.path.exists(f'{self.name}.mol2'),f'Cannot parameterize molecule {self.name} without {self.name}.mol2'
        if outname=='':
            outname=f'{self.name}-p'
        GAFFParameterize(self.name,outname,**kwargs)
        self.read_topology(f'{outname}.top')
        self.read_coords(f'{outname}.gro')
        self.read_coords(f'{outname}.mol2')
        assert self.cstale=='',f'Error: {self.cstale} coords are stale'

    def calculate_sea(self):
        ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
        n=self.name
        boxsize=np.array(self.Coords['gro'].maxspan())+2*np.ones(3)
        pfs.checkout('nvt-sea.mdp')
        for ex in ['top','itp','gro']:
            pfs.checkout(f'{n}-p.{ex}')
        logging.info(f'Hot md running...output to {n}-p-sea')
        grompp_and_mdrun(gro=f'{n}-p',top=f'{n}-p',
                        mdp='nvt-sea',out=f'{n}-p-sea',boxSize=boxsize)
        self.Coords['gro'].set_attribute('sea-idx',analyze_sea(f'{n}-p-sea'))
        self.Coords['gro'].write_sea(f'{n}-p.sea')

    def read_sea(self,filename):
        assert os.path.exists(filename),f'SEA file {filename} not found.'
        self.Coords['gro'].read_sea(filename)

    def generate(self,outname='',available_molecules={},**kwargs):
        logging.info(f'Generating Molecule {self.name}')
        if outname=='':
            outname=f'{self.name}-p'
        if self.generator:
            R=self.generator
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            can_react=all([a in available_molecules for a in R.reactants.values()])
            if not can_react:
                raise Exception(f'Cannot generate {self.name} because reactants not available')
            # add z attribute to molecules in input list based on bonds
            for b in R.bonds:
                atoms=[R.atoms[i] for i in b.atoms]
                mols=[available_molecules[R.reactants[a.reactant]] for a in atoms]
                for a,m in zip(atoms,mols):
                    if not 'z' in m.Coords['gro']:
                        logging.info(f'Initializing z attribute of all atoms in {m.name}')
                        m.Coords['gro']['z']=np.zeros(m.Coords['gro'].shape[0])
                    idx=m.Coords['gro'].get_idx({'atomName':a.atom,'resNum':a.resid})-1
                    logging.info(f'Modifing z attribute of molecule {m.name} atom {idx+1}')
                    m.Coords['gro']['z'].iloc[idx]=a.z

            self.reactants={}
            for n,r in R.reactants.items():
                self.reactants[n]=deepcopy(available_molecules[r])
            base=None
            for b in R.bonds:
                # every bond names exactly two atoms, A and B
                # here we associate the identifiers A and B with
                # their entry in the atoms dictionary for this reaction
                A,B=[R.atoms[i] for i in b.atoms]
                mA,mB=[self.reactants[a.reactant] for a in [A,B]]
                if not base:
                    base=mA          
                Aidx=m.Coords['gro'].get_idx({'atomName':A.atom,'resNum':A.resid})
                Bidx=m.Coords['gro'].get_idx({'atomName':B.atom,'resNum':B.resid})
                mA.react_with(mB,Aidx,Bidx)
            self.merge(base)
            self.toggle('gro')
            self.sync_coords()
            self.write_mol2()
        else:
            logging.info(f'Using {self.name}.mol2 as a source.')
        self.parameterize(outname,**kwargs)

    def toggle(self,t):  # toggle stalecoords after update of type t
        def other(t):
            return 'gro' if t=='mol2' else 'mol2'
        if self.cstale==t: # indicates we just updated stale coords
            self.cstale='' # none are stale
        else:
            if self.cstale=='': # none are stale, but we just updated t ...
                self.cstale=other(t) #... so other is now stale
            else: # other is stale (cstale is not t nor '') ...
                pass # ... do nothing
    
    def sync_coords(self):
        ''' keep the gro and mol2 format coordinates in sync '''
        if self.cstale=='gro':
            self.Coords['gro'].copy_coords(self.Coords['mol2'])
            self.toggle('gro')
        elif self.cstale=='mol2':
            self.Coords['mol2'].copy_coords(self.Coords['gro'])
            self.toggle('mol2')
        else:
            pass

    def read_topology(self,filename):
        assert os.path.exists(filename),f'Topology file {filename} not found.'
        if self.Topology:
            logging.warning(f'Overwriting topology of monomer {self.name} from file {filename}')
        self.Topology=Topology.read_gro(filename)

    def read_coords(self,filename):
        assert os.path.exists(filename),f'Coordinate file {filename} not found.'
        if self.Coords:
            logging.warning(f'Overwriting coordinates of monomer {self.name} from file {filename}')
        basename,ext=os.path.splitext(filename)
        if ext=='.mol2':
            self.Coords['mol2']=Coordinates.read_mol2(filename)
            self.toggle('mol2')
        elif ext=='.gro':
            self.Coords['gro']=Coordinates.read_gro(filename)
            self.toggle('gro')
        else:
            raise Exception(f'Coordinate filename extension {ext} is not recognized.')

    def merge(self,other):
        if not self.Topology:
            self.Topology=other.Topology
        else:
            self.Topology.merge(other.Topology)
        for typ in ['gro','mol2']:
            if not typ in self.Coords:
                if typ in other.Coords:
                    self.Coords[typ]=other.Coords[typ]
            else:
                if typ in other.Coords:
                    self.Coords[typ].merge(other.Coords[typ])

    def update_topology(self,t):
        self.Topology.merge(t)

    def update_coords(self,c):
        self.Coords.copy_coords(c)

    def react_with(self,other,myidx,otidx,**kwargs):
        myC=self.Coords['gro']
        myA=myC['atoms']
        myT=self.Topology
        otC=other.Coords['gro']
        otA=otC['atoms']
        otT=other.Topology

        mypartners=myT.bondlist.partners_of(myidx)
        otpartners=otT.bondlist.partners_of(otidx)
        myHpartners=[k for k,v in zip(mypartners,[myA[myA['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')]
        otHpartners=[k for k,v in zip(otpartners,[otA[otA['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')]
        assert len(myHpartners)>0,f'Error: atom {myidx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {otidx} does not have a deletable H atom!'

        Ri=myC.get_R(myidx)
        Rj=otC.get_R(otidx)
        if self!=other: # not an intramolecular reaction
            # TODO: optimize position/orientation of other
            pass

        # TODO: get the H on myidx that is closest to otidx
        # and the H on otidx that is closest to myidx
        # these are the H's to delete

        pass



def get_conn(mol):
    conn=[]
    for mrname,mr in mol.reactive_atoms.items():
        ''' check to see if this atom's symmetry partners are already on the list '''
        donotadd=False
        adf=mol.Coords['active'].D['atoms']
        if 'sea-idx' in adf:
            isea=adf[adf['atomName']==mrname]['sea-idx'].values[0]
            logging.info(f'get_conn: looking for other atoms with sea-idx = {isea}')
            for s in adf[adf['sea-idx']==isea]['atomName']:
                logging.info(f'...found {s}')
                if s in conn:
                    donotadd=True
        if not donotadd:
            conn.append(mrname)
    return conn

def get_oligomers(m,n):
    ''' m and n are Monomer instances (defined in here) '''
    if not isinstance(m,Monomer) or not isinstance(n,Monomer):
        raise Exception(f'react_mol2 needs Monomers not {type(m)} and {type(n)}')
    if not 'active' in m.Topology or not 'active' in n.Topology:
        raise Exception('react_mol2 needs Monomers with active topologies')
    if not 'active' in m.Coords or not 'active' in n.Coords:
        raise Exception('react_mol2 needs Monomers with active coordinates')

    oligdict=_oligomerize(m,n)
    # check to see if we need to consider n as a basemol
    # we do if any reactive atom on n has z>1
    nbase=any([n.reactive_atoms[i].z>1 for i in get_conn(n)])
    if nbase:
        oligdict.update(_oligomerize(n,m))
    return oligdict

def _oligomerize(basemol,othermol):
    oligdict={}
    logging.info(f'react_mol2: basemol {basemol.name} and othermol {othermol.name}')
    # list of all asymmetric reactive atoms on base
    basera=get_conn(basemol)
    # number of connections for each of those reactive atoms
    baseconn=[basemol.reactive_atoms[i].z for i in basera]
    # list of all asymmetric reactive atoms on other
    otherra=get_conn(othermol)
    # options for a connection on basemol are empty ('') or any one of 
    # asymmetric atoms on other
    otherra=['']+otherra
    # enumerate all symmetry-unique configurations of oligomers
    # on *each* reactive atom on basemol *independently*, for all
    # connections on each reactive atom
    basearr=[]
    for z in baseconn:
        # this reactive atom has z connections, each of which can be
        # 'occupied' by one member of otherra in all possible combinations
        basearr.append(list(combinations_with_replacement(otherra,z)))
    # for n,c in zip(basera,basearr):
    #     print(f'{n} can have following connection configurations:',c)
    # enumerate all unique configurations of connections on all reactive
    # atoms.  This is relevant if there is more than one asymmetric 
    # reactive atom on basemol.
    o=product(*basearr)
    # the first one is one where all connections on all reactive atoms
    # of basemol are empty, so skip it
    next(o)
    for oligo in o:
        # logging.debug('making oligo',oligo)
        # make working copy of basemol coordinates
        wc=deepcopy(basemol.Coords['active'])
        stoich=[(basemol,1)]
        oname=basemol.name
        for c,b in zip(oligo,basera):
            # print(f'establishing connection(s) to atom {b} of {basemol.name}:')
            nconn=len([x for x in c if x!=''])
            if nconn>0:
                oname+=f'@{b}-'+','.join([f'{othermol.name}#{a}' for a in c if a!=''])
            oc=0
            for a in c:
                # print(f'  to atom {a} of {othermol.name}')
                if a != '':
                    owc=deepcopy(othermol.Coords['active'])
                    oc+=1
                    # bring copy of coords from othermol into 
                    # working copy of basemol
                    wc.bond_to(owc,acc=b,don=a)
            if oc>0:
                stoich.append((othermol,oc))
        # print('-> prefix',oname)
        # wc.write_mol2(f'OLIG-{oname}.mol2')
        oligdict[oname]=Oligomer(coord=wc,stoich=stoich)

    return oligdict