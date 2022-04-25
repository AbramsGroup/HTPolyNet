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
        self.Topology=Topology()
        self.Coords=Coordinates()
        self.generator=generator

    def __str__(self):
        restr=f'{self.name} '
        if not self.Topology:
            restr+=f'(empty topology) '
        if len(self.Coords)==0:
            restr+=f'(empty coords) '
        return restr+'\n'

    def num_atoms(self):
        if hasattr(self.Coords,'N'):
            return self.Coords.N
        else:
            return self.Coords.shape[0]

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
        # self.read_coords(f'{outname}.gro')
        self.read_coords(f'{outname}.mol2')
        #assert self.cstale=='',f'Error: {self.cstale} coords are stale'

    def calculate_sea(self):
        ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
        n=self.name
        boxsize=np.array(self.Coords.maxspan())+2*np.ones(3)
        pfs.checkout('mdp/nvt-sea.mdp')
        for ex in ['top','itp','gro']:
            pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        logging.info(f'Hot md running...output to {n}-sea')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='nvt-sea',out=f'{n}-sea',boxSize=boxsize)
        self.set_atomset_attribute('sea-idx',analyze_sea(f'{n}-sea'))
        self.write_atomset_attributes(['sea-idx'],f'{n}.sea')

    def minimize(self,outname='',**kwargs):
        if outname=='':
            outname=f'{self.name}'
        n=self.name
        boxsize=np.array(self.Coords.maxspan())+2*np.ones(3)
        pfs.checkout('mdp/em-single-molecule.mdp')
        if 'checkout_required' in kwargs:
            for ex in ['top','itp','gro']:
                pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        logging.info(f'Hot md running...output to {n}-sea')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='em-single-molecule',out=f'{outname}',boxSize=boxsize)
        self.read_coords(f'{n}.gro')

    def has_atom_attributes(self,attributes):
        return self.Coords.has_atom_attributes(attributes)

    def set_atomset_attribute(self,attribute,srs):
        self.Coords.set_atomset_attribute(attribute,srs)

    def write_atomset_attributes(self,name,filename):
        if self.has_atom_attributes(name):
            self.Coords.write_atomset_attributes(name,filename)

    def set_atom_attribute(self,name,value,attributes):
        self.Coords.set_atom_attribute(name,value,attributes)

    def read_atomset_attributes(self,filename,attributes=[]):
        ''' filename is expected to be a column-oriented file with a 'globalIdx' header and zero or more headed columns explicitly named in attributes[] '''
        self.Coords.read_atomset_attributes(filename,attributes=attributes)

    def generate(self,outname='',available_molecules={},**kwargs):
        logging.info(f'Generating Molecule {self.name}')
        if outname=='':
            outname=f'{self.name}'
        ''' wipe topology and coordinate '''
        if not self.Topology.empty:
            logging.info(f'Called generate() on {self.name} with non-empty topology')
            self.Topology=Topology()
        if not self.Coords.empty:
            logging.info(f'Called generate() on {self.name} with non-empty coordinates')
            self.Coords=Coordinates()

        if self.generator:
            R=self.generator
            assert type(R)==Reaction,'Can only handle Reaction-type generators at the moment'
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            can_react=all([a in available_molecules for a in R.reactants.values()])
            if not can_react:
                raise Exception(f'Cannot generate {self.name} because reactants not available')
            # add z attribute to molecules in input list based on bonds
#            print(R)
            for b in R.bonds:
                atoms=[R.atoms[i] for i in b['atoms']]
                mols=[available_molecules[R.reactants[a['reactant']]] for a in atoms]
                for a,m in zip(atoms,mols):
                    if not m.has_atom_attributes(['z']):
#                    if not 'z' in m.Coords['mol2'].D['atoms']:
                        logging.info(f'Initializing z attribute of all atoms in molecule {m.name}')
                        m.set_atomset_attribute('z',np.zeros(m.num_atoms()))
                    logging.info(f'Setting z of {m.name}-R{a["resid"]}-{a["atom"]} to {a["z"]}')
                    m.set_atom_attribute('z',a['z'],{'atomName':a['atom'],'resNum':a['resid']})
                    logging.info(f'Are we setting z values of symmetry-equivalent atoms?')
                    if m.Coords.has_atom_attributes(['sea-idx']):
                        Aidx=m.Coords.get_atom_attribute('globalIdx',{'atomName':a['atom'],'resNum':a['resid']})
                        Asea=m.Coords.get_atom_attribute('sea-idx',{'globalIdx':Aidx})
                        Aclu=m.Coords.get_atoms_w_attribute('globalIdx',{'sea-idx':Asea})
                        logging.info(f'Yes.  Atom {Aidx} belongs to sea-class {Asea} whose total member idx list is {Aclu}')
                        for aa in Aclu:
                            if aa!=Aidx:
                                an=m.Coords.get_atom_attribute('atomName',{'globalIdx':aa})
                                rn=m.Coords.get_atom_attribute('resNum',{'globalIdx':aa})
                                logging.info(f'Setting z of {m.name}-R{rn}-{an} to {a["z"]}')
                                m.set_atom_attribute('z',a['z'],{'globalIdx':aa})
            # local copies of all reactant molecules
            reactants={}
            for n,r in R.reactants.items():
                reactants[n]=deepcopy(available_molecules[r])
            bases=[]
            for b in R.bonds:
                # every bond names exactly two atoms, A and B
                # here we associate the identifiers A and B with
                # their entry in the atoms dictionary for this reaction
                A,B=[R.atoms[i] for i in b['atoms']]
                mA,mB=[reactants[a['reactant']] for a in [A,B]]
                aan=A['atom']
                ban=B['atom']
                arn=A['resid']
                brn=B['resid']
                Aidx=mA.Coords.get_idx({'atomName':A['atom'],'resNum':A['resid']})
                Bidx=mB.Coords.get_idx({'atomName':B['atom'],'resNum':B['resid']})
                logging.info(f'Explicitly, I think {mA.name}-R{arn}-{aan}({Aidx}) should bond to {mB.name}-R{brn}-{ban}({Bidx})')
                logging.info(f'Creating new bond {Aidx}-{Bidx} in molecule {mA.name}')
                Aidx,Bidx=mA.new_bond(mB,Aidx,Bidx)
                aan=mA.Coords.get_atom_attribute('atomName',{'globalIdx':Aidx})
                arn=mA.Coords.get_atom_attribute('resNum',{'globalIdx':Aidx})
                ban=mA.Coords.get_atom_attribute('atomName',{'globalIdx':Bidx})
                brn=mA.Coords.get_atom_attribute('resNum',{'globalIdx':Bidx})
                logging.info(f'Due to post-bonding reindexing, the two bonded atoms are now {mA.name}-R{arn}-{aan}({Aidx}) and {mB.name}-R{brn}-{ban}({Bidx})')
                logging.info(f'Are we automatically making symmetry-equivalent bonds in molecule {mA.name}?')
                if mA==mB and mA.Coords.has_atom_attributes(['sea-idx']):
                    Asea=mA.Coords.get_atom_attribute('sea-idx',{'globalIdx':Aidx})
                    Bsea=mA.Coords.get_atom_attribute('sea-idx',{'globalIdx':Bidx})
                    Aclu=mA.Coords.get_atoms_w_attribute('globalIdx',{'sea-idx':Asea})
                    Bclu=mA.Coords.get_atoms_w_attribute('globalIdx',{'sea-idx':Bsea})
                    Aclu=np.delete(Aclu,np.where(Aclu==Aidx))
                    Bclu=np.delete(Bclu,np.where(Bclu==Bidx))
                    logging.info(f'Yes.  Atom {Aidx} belongs to sea-cluster {Asea} whose other member(s) is/are {Aclu}')
                    logging.info(f'      Atom {Bidx} belongs to sea-cluster {Bsea} whose other member(s) is/are {Bclu}')
                    if mA==mB:
                        for aa,bb in zip(Aclu,Bclu):
                            aan=mA.Coords.get_atom_attribute('atomName',{'globalIdx':aa})
                            arn=mA.Coords.get_atom_attribute('resNum',{'globalIdx':aa})
                            ban=mA.Coords.get_atom_attribute('atomName',{'globalIdx':bb})
                            brn=mA.Coords.get_atom_attribute('resNum',{'globalIdx':bb})
                            logging.info(f'By symmetry, I think {mA.name}-R{arn}-{aan}({aa}) should bond to {mB.name}-R{brn}-{ban}({bb})')
                            new_idx=mA.new_bond(mB,aa,bb,return_idx_of=[aa,bb]+list(Aclu)+list(Bclu))
                            aa,bb=new_idx[0:2]
                            Aclu=new_idx[2:2+len(Aclu)]
                            Bclu=new_idx[2+len(Aclu):2+len(Aclu)+len(Bclu)]
                else:
                    logging.info('No.')
                if len(bases)==0 or mA not in bases:
                    bases.append(mA)
            assert len(bases)==1,f'Error: Reaction {R.name} results in more than one molecular fragment product'
            base=bases[0]
            self.merge(base)
            self.write_mol2(filename=f'{self.name}.mol2')
        else:
            logging.info(f'Using existing molecules/inputs/{self.name}.mol2 as a source.')
            pfs.checkout(f'molecules/inputs/{self.name}.mol2')

        self.parameterize(outname,**kwargs)
        self.minimize(outname,**kwargs)  

    def inherit_sea_from_reactants(self,molecules,sea_list):
        R=self.generator
        if R and len(sea_list)>0:
            self.Coords.A['sea-idx']=np.arange(self.Coords.A.shape[0])
            logging.debug(f'Inheriting sea: {self.name} using reaction {R.name} and sea_list {sea_list}')
            logging.debug(f'Coords prior to inheriting sea:\n'+self.Coords.A.to_string())

            donors=[molecules[x] for x in R.reactants.values() if x in sea_list]
            for D in donors:
                logging.debug(f'Donor {D.name}')
                dfsea=D.Coords.A[['resNum','resName','atomName','sea-idx']]
                logging.debug(f'  coords subset:\n'+dfsea.to_string())
                for i,a in self.Coords.A.iterrows():
                    if a['resName']==D.name:
                        logging.debug(f'Inheriting from donor {D.name} to {a["resName"]}')
                        sealist=dfsea[(dfsea['resName']==a['resName'])&(dfsea['atomName']==a['atomName'])]['sea-idx'].values
                        logging.debug(f"atom {a['resName']} {a['atomName']} sea: {sealist}")
                        if len(sealist)==0:
                            logging.debug(f"atom {a['resName']} {a['atomName']} not found in sea")
                        logging.debug(f"setting sea of {a['resName']} {a['atomName']} to {sealist[0]}")
                        self.Coords.set_atom_attribute('sea-idx',sealist[0],{'resName':a['resName'],'atomName':a['atomName']})
            logging.debug(f'Molecule {self.name} after inheriting sea:\n'+self.Coords.A.to_string())

    def read_topology(self,filename):
        assert os.path.exists(filename),f'Topology file {filename} not found.'
        if not self.Topology.empty:
            logging.warning(f'Overwriting topology of monomer {self.name} from file {filename}')
        sv=pd.DataFrame
        if 'mol2_bonds' in self.Topology.D:
            sv=self.Topology.D['mol2_bonds'].copy()
        self.Topology=Topology.read_gro(filename)
        if not sv.empty:
            self.Topology.D['mol2_bonds']=sv
            self.Topology.bond_source_check()

    def read_coords(self,filename):
        assert os.path.exists(filename),f'Coordinate file {filename} not found.'
        if not self.Coords.empty:
            logging.warning(f'Overwriting coordinates of monomer {self.name} from file {filename}')
        basename,ext=os.path.splitext(filename)
        if ext=='.mol2':
            self.read_mol2(filename)
        elif ext=='.gro':
            self.read_gro(filename)
#            print(filename,self.Coords.A.to_string())
        else:
            raise Exception(f'Coordinate filename extension {ext} is not recognized.')
        if not self.Coords.mol2_bonds.empty: # we just read in some bonding info
            self.Topology.D['mol2_bonds']=self.Coords.mol2_bonds.copy()
            if 'bonds' in self.Topology.D: # this molecule has bonds already defined in its topology
                self.Topology.bond_source_check()

    def read_mol2(self,filename):
        self.Coords=Coordinates.read_mol2(filename)

    def read_gro(self,filename):
        self.Coords=Coordinates.read_gro(filename)

    def write_mol2(self,filename,molname=''):
        if molname=='':
            molname=self.name
        if 'mol2_bonds' in self.Topology.D:
            self.Coords.write_mol2(filename,molname=molname,bondsDF=self.Topology.D['mol2_bonds'])
        else:
            self.Coords.write_mol2(filename,molname=molname)

    def merge(self,other):
        self.merge_topologies(other)
        self.merge_coordinates(other)
    
    def merge_topologies(self,other):
        if not self.Topology:
            self.Topology=other.Topology
        else:
            self.Topology.merge(other.Topology)

    def merge_coordinates(self,other):
        self.Coords.merge(other.Coords)

    def update_topology(self,t):
        self.Topology.merge(t)

    def update_coords(self,c):
        self.Coords.copy_coords(c)

    def add_bonds(self,pairs=[]):
        self.Topology.add_bonds(pairs)

    def delete_atoms(self,idx=[],return_idx_of=[]):
        new_idx=self.Topology.delete_atoms(idx,return_idx_of=return_idx_of)
        self.Coords.delete_atoms(idx)
        return new_idx

    def new_bond(self,other,myidx,otidx,**kwargs):
        myC=self.Coords
        myA=myC.A
        myT=self.Topology
        otC=other.Coords
        otA=otC.A
        otT=other.Topology
        return_idx_of=kwargs.get('return_idx_of',[myidx,otidx])
        logging.info(f'new_bond: will request updated global indices of {return_idx_of}')
        assert not myT.empty,f'Error: empty topology -- cannot make a new bond'
        assert not otT.empty,f'Error: empty topology -- cannot make a new bond'
        myz=myC.get_atom_attribute('z',{'globalIdx':myidx})
        otz=otC.get_atom_attribute('z',{'globalIdx':otidx})
        assert myz>0 and otz>0,f'No bond permissible: z({myidx}):{myz}, z({otidx}):{otz}'

        mypartners=myT.bondlist.partners_of(myidx)
        otpartners=otT.bondlist.partners_of(otidx)
        logging.info(f'Partners of {myidx} {mypartners}')
        logging.info(f'Partners of {otidx} {otpartners}')
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
                if self!=other: # this is intermolecular--most likely building a molecule
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
    
        if self!=other:
            self.Topology.merge(other.Topology)
        self.add_bonds(pairs=[(myidx,otidx)])
        myC.set_atom_attribute('z',myz-1,{'globalIdx':myidx})
        otC.set_atom_attribute('z',otz-1,{'globalIdx':otidx})
        new_idx=self.delete_atoms(idx=[myH,otH],return_idx_of=return_idx_of)
        self.Topology.bond_source_check()
        return new_idx
