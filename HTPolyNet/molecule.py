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

class Reaction:
    ''' reactions may only be initialized after monomers '''
    def __init__(self,jsondict):
        self.jsondict=jsondict
        ''' attributes expected in cfg file '''
        self.name=jsondict.get('name','')
        self.atoms=jsondict.get('atoms',{})
        self.bonds=jsondict.get('bonds',{})
        self.reactants=jsondict.get('reactants',{})
        self.product=jsondict.get('product','')
        self.restrictions=jsondict.get('restrictions',{})
        self.stage=jsondict.get('stage','')

    def __str__(self):
        return f'Reaction "{self.name}"'

class Molecule:
    def __init__(self,name='',generator=None):
        self.name=name
        self.Topology=None
        self.Coords={} # ['gro', 'mol2]
        self.generator=generator
        self.cstale=''  # ['','gro','mol2']

    def __str__(self):
        restr=f'{self.name} '
        if not self.Topology:
            restr+=f'(empty topology) '
        if len(self.Coords)==0:
            restr+=f'(empty coords) '
        return restr+'\n'

    def num_atoms(self):
        N=[]
        for typ in ['gro','mol2']:
            n=self.Coords[typ].N
            if not n in N:
                N.append(n)
        assert len(N)==1,'gro and mol2 coordinates are not in sync'
        return N[0]
        
    def previously_parameterized(self):
        rval=True
        for ext in ['mol2','top','itp','gro']:
            rval=rval and pfs.exists(os.path.join('molecules/parameterized',f'{self.name}.{ext}'))
        return rval

    def parameterize(self,outname='',**kwargs):
        assert os.path.exists(f'{self.name}.mol2'),f'Cannot parameterize molecule {self.name} without {self.name}.mol2'
        if outname=='':
            outname=f'{self.name}'
        GAFFParameterize(self.name,outname,**kwargs)
        self.read_topology(f'{outname}.top')
        self.read_coords(f'{outname}.gro')
        self.read_coords(f'{outname}.mol2')
        assert self.cstale=='',f'Error: {self.cstale} coords are stale'

    def calculate_sea(self):
        ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
        n=self.name
        boxsize=np.array(self.Coords['gro'].maxspan())+2*np.ones(3)
        pfs.checkout('mdp/nvt-sea.mdp')
        for ex in ['top','itp','gro']:
            pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        logging.info(f'Hot md running...output to {n}-sea')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='nvt-sea',out=f'{n}-sea',boxSize=boxsize)
        self.Coords['gro'].set_atomset_attribute(attribute='sea-idx',srs=analyze_sea(f'{n}-sea'))
        self.Coords['gro'].write_atomset_attributes(attributes=['sea-idx'],filename=f'{n}.sea')

    def minimize(self,outname='',**kwargs):
        if outname=='':
            outname=f'{self.name}'
        n=self.name
        boxsize=np.array(self.Coords['gro'].maxspan())+2*np.ones(3)
        pfs.checkout('mdp/em-single-molecule.mdp')
        if 'checkout_required' in kwargs:
            for ex in ['top','itp','gro']:
                pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        logging.info(f'Hot md running...output to {n}-sea')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='em-single-molecule',out=f'{outname}',boxSize=boxsize)
        self.read_coords(f'{n}.gro')
        self.toggle('gro')
        self.sync_coords()

    def has_atom_attribute(self,attribute):
        hasit=True
        for typ in ['gro','mol2']:
            hasit = hasit and self.Coords[typ].has_atom_attribute(attribute)
        return hasit

    def set_atomset_attribute(self,name,srs):
        for typ in ['gro','mol2']:
            self.Coords[typ].set_atomset_attribute(name,srs)

    def set_atom_attribute(self,name,value,attributes):
        for typ in ['gro','mol2']:
            self.Coords[typ].set_atom_attribute(name,value,attributes)

    def generate(self,outname='',available_molecules={},**kwargs):
        logging.info(f'Generating Molecule {self.name}')
        if outname=='':
            outname=f'{self.name}'
        if self.generator:
            R=self.generator
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            can_react=all([a in available_molecules for a in R.reactants.values()])
            if not can_react:
                raise Exception(f'Cannot generate {self.name} because reactants not available')
            # add z attribute to molecules in input list based on bonds
            print(R)
            for b in R.bonds:
                atoms=[R.atoms[i] for i in b['atoms']]
                mols=[available_molecules[R.reactants[a['reactant']]] for a in atoms]
                for a,m in zip(atoms,mols):
                    if not m.has_atom_attribute('z'):
#                    if not 'z' in m.Coords['mol2'].D['atoms']:
                        logging.info(f'Initializing z attribute of all atoms in molecule {m.name}')
                        m.set_atomset_attribute('z',np.zeros(m.num_atoms()))
                    logging.info(f'Modifing z attribute of molecule {m.name} resnum {a["resid"]} atom {a["atom"]}')
                    m.set_atom_attribute('z',a['z'],{'atomName':a['atom'],'resNum':a['resid']})

            # local copies of all reactant molecules
            reactants={}
            for n,r in R.reactants.items():
                reactants[n]=deepcopy(available_molecules[r])
            bases=[]
            # TODO: enforce any symmetries
            for b in R.bonds:
                # every bond names exactly two atoms, A and B
                # here we associate the identifiers A and B with
                # their entry in the atoms dictionary for this reaction
                A,B=[R.atoms[i] for i in b['atoms']]
                mA,mB=[reactants[a['reactant']] for a in [A,B]]
                Aidx=mA.Coords['gro'].get_idx({'atomName':A['atom'],'resNum':A['resid']})
                Bidx=mB.Coords['gro'].get_idx({'atomName':B['atom'],'resNum':B['resid']})
                mA.new_bond(mB,Aidx,Bidx)
                if len(bases)==0 or mA not in bases:
                    bases.append(mA)
            assert len(bases)==1,f'Error: Reaction {R.name} results in more than one molecular fragment product'
            base=bases[0]
            # at this point, base.Coord['mol'] contains our coordinates, and base.Coord['gro'] still
            # reflects the original base, so we need to fully generate a new base.Coord['gro'] from 
            # the base.Coord['mol2']
            base.regenerate_coordinates('gro')
            self.merge_coordinates(base)
            self.Coords['mol2'].write_mol2(filename=f'{self.name}.mol2')
        else:
            logging.info(f'Using existing molecules/inputs/{self.name}.mol2 as a source.')
            pfs.checkout(f'molecules/inputs/{self.name}.mol2')

        self.parameterize(outname,**kwargs)
        self.minimize(outname,**kwargs)

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

    def regenerate_coordinates(self,typ):
        assert typ in ['gro','mol2']
        if typ=='gro':  # we have to generate gro from mol2
            self.Coords['gro'].D['atoms']=pd.DataFrame()
            df=self.Coords['gro'].D['atoms']
            for c in Coordinates.gro_colnames[:-3]:
                df[c]=self.Coords['mol2'].D['atoms'][c]
            for c in ['posX','posY','posZ']:
                df[c]/=10.0 # nm!
        else: # we have to generate mol2 from gro
            self.Coords['mol2'].D['atoms']=pd.DataFrame()
            df=self.Coords['mol2'].D['atoms']
            for c in Coordinates.gro_colnames[:-3]:
                df[c]=self.Coords['gro'].D['atoms'][c]
            for c in ['posX','posY','posZ']:
                df[c]*=10.0 # A!            

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
        self.merge_topologies(other)
        self.merge_coordinates(other)
    
    def merge_topologies(self,other):
        if not self.Topology:
            self.Topology=other.Topology
        else:
            self.Topology.merge(other.Topology)

    def merge_coordinates(self,other):
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

    def new_bond(self,other,myidx,otidx,**kwargs):
        myC=self.Coords['mol2']
        myA=myC.D['atoms']
        otC=other.Coords['mol2']
        otA=otC.D['atoms']
        myz=myC.get_atom_attribute('z',{'globalIdx':myidx})
        otz=otC.get_atom_attribute('z',{'globalIdx':otidx})
        assert myz>0 and otz>0,f'No bond permissible: z({myidx}):{myz}, z({otidx}):{otz}'

        mypartners=myC.bondlist.partners_of(myidx)
        otpartners=otC.bondlist.partners_of(otidx)
        myHpartners=[k for k,v in zip(mypartners,[myA[myA['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')]
        otHpartners=[k for k,v in zip(otpartners,[otA[otA['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')]
        assert len(myHpartners)>0,f'Error: atom {myidx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {otidx} does not have a deletable H atom!'
            
        Ri=myC.get_R(myidx)
        Rj=otC.get_R(otidx)
        minHH=(1.e9,-1,-1)
        if self!=other:
            overall_maximum=(-1.e9,-1,-1)
        ''' Identify the H atom on each atom in be pair to
            be bonded that result in the "optimum" choice
            for deletion or, if this is intermolecular,
            the "optimum" choice for the location/orientation
            of the other. '''
        totc={}
        for myH in myHpartners:
            totc[myH]={}
            Rh=myC.get_R(myH)
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            for otH in otHpartners:
                totc[myH][otH]=deepcopy(otC)
                Rk=totc[myH][otH].get_R(otH)
                logging.debug(f'   otH {otH} Rk {Rk}')
                Rkj=Rk-Rj
                Rkj*=1.0/np.linalg.norm(Rkj)
                Rhk=Rh-Rk
                rhk=np.linalg.norm(Rhk)
                if self!=other: # this is intermolecular
                    cp=np.cross(Rkj,Rih)
                    c=np.dot(Rkj,Rih)
                    v=np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
                    v2=np.dot(v,v)
                    I=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                    # R is the rotation matrix that will rotate donb to align with accb
                    R=I+v+v2/(1.+c)
                    # rotate translate all donor atoms!
                    totc[myH][otH].rotate(R)
                    Rk=totc[myH][otH].get_R(otH)
                    # overlap the other H atom with self's 
                    # reactive atom by translation
                    Rik=Ri-Rk
                    totc[myH][otH].translate(Rik)
                    minD=myC.minimum_distance(totc[myH][otH],self_excludes=[myH],other_excludes=[otH])
                    if minD>overall_maximum[0]:
                        overall_maximum=(minD,myH,otH)
                else:
                    if rhk<minHH[0]:
                        minHH=(rhk,myH,otH)
        if self!=other:
            minD,myH,otH=overall_maximum
            shifts=myC.merge(totc[myH][otH])
            idxshift=shifts[0]
            otidx+=idxshift
            otH+=idxshift
        else:
            mhh,myH,otH=minHH
            logging.info(f'Executing intramolecular reaction  in {self.name}')
            logging.info(f'Two Hs closest to each other are {myH} and {otH}')
    
        myC.add_bonds(pairs=[(myidx,otidx)])
        myC.set_atom_attribute('z',myz-1,{'globalIdx':myidx})
        otC.set_atom_attribute('z',otz-1,{'globalIdx':otidx})
        myC.delete_atoms(idx=[myH,otH])
