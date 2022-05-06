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

def get_base_reactants(mname,mdict):
    mol=mdict[mname]
    if not mol.generator:
        return [mol]
    mine=[]
    for r in mol.generator.reactants.values():
        mine.extend(get_base_reactants(r,mdict))
    return mine

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
        self.sym=0 # symmetry code: 0 means this reaction's molecules
                   # are the ones that are parameterized

    def __str__(self):
        return f'Reaction "{self.name}"'

class Molecule:
    def __init__(self,name='',generator=None):
        self.name=name
        self.Topology=Topology()
        self.Coords=Coordinates()
        self.generator=generator
        self.sequence=[]
        self.origin=None

    def __str__(self):
        restr=f'{self.name} '
        if not self.Topology:
            restr+=f'(empty topology) '
        if len(self.Coords)==0:
            restr+=f'(empty coords) '
        return restr+'\n'

    def set_origin(self,value):
        self.origin=value

    def get_origin(self):
        return self.origin

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
        sea_srs=analyze_sea(f'{n}-sea')
        self.set_atomset_attribute('sea-idx',sea_srs)
        # now, let's look at the topology parameters and see if this symmetry holds in them

    def analyze_sea_topology(self):
        tadf=self.Topology.D['atoms']
        cadf=self.Coords.A
        maxsea=cadf['sea-idx'].max()
        for i in range(maxsea+1):
            sea_indexes=cadf[cadf['sea-idx']==i]['globalIdx'].to_list()
            sea_cls=tadf[tadf['nr'].isin(sea_indexes)]
            # logging.debug(f'{self.name} symmetry class {i}:\n{sea_cls.to_string()}')
            for attr in ['type', 'residue', 'charge', 'mass']:
                values=sea_cls[attr].values
                flg=values[0]
                chk=all(values==flg)
                if not chk:
                    logging.debug(f'Error: atoms in symmetry class {i} have different values of {attr}')

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

    def propagate_z(self,reactions,mdict):
        ''' Assign all z values in Coords.A dataframes for each molecule '''
        logging.debug(f'propagate_z called for {self.name}')
        self.set_atomset_attribute('z',[0]*self.Coords.A.shape[0])
        if self.generator:  # this molecule is a product
            R=self.generator
            logging.debug(f'{self.name} was generated by {R.name}')
            resids=list(set(self.Coords.A['resNum']))
            resnames=[]
            logging.debug(f'{self.name} has resids {resids}')
            for resid in resids:
                resnames.append(list(set(self.Coords.A[self.Coords.A['resNum']==resid]['resName']))[0])
            logging.debug(f'{self.name} has resnames {resnames}')
            for resid,resname in zip(resids,resnames):
                logging.debug(f'taking z from {resname}')
                mol=mdict[resname]
                adf=mol.Coords.A
                zs=adf[adf['z']>0]
                for i,r in zs.iterrows():
                    atomName=r['atomName']
                    resNum=r['resNum']
                    resName=r['resName']
                    z=r['z']
                    self.Coords.set_atom_attribute('z',z,{'atomName':atomName,'resName':resName})
            ''' correct from self's bondlist '''
            bondlist=self.Topology.bondlist
            adf=self.Coords.A
            reactive_atoms_idx=list(adf[adf['z']>0]['globalIdx'])
            # logging.debug(f'reactive_atoms_idx {reactive_atoms_idx}')
            for i in range(len(reactive_atoms_idx)):
                ix=reactive_atoms_idx[i]
                for j in range(i,len(reactive_atoms_idx)):
                    jx=reactive_atoms_idx[j]
                    iz=self.Coords.get_atom_attribute('z',{'globalIdx':ix})
                    jz=self.Coords.get_atom_attribute('z',{'globalIdx':jx})
                    if bondlist.are_bonded(ix,jx):
                        self.Coords.set_atom_attribute('z',iz-1,{'globalIdx':ix})
                        self.Coords.set_atom_attribute('z',jz-1,{'globalIdx':jx})
            # logging.debug(f'after propagate_z product {self.name}:\n{self.Coords.A.to_string()}')
        else:  # this molecule is not a product -- was generated using an input mol2
            ''' find any reaction in which this is a reactant as long as it is not a product '''
            atoms=[]
            for R in reactions:
                if self.name in R.reactants.values():
                    ridx=[k for k,v in R.reactants.items() if v==self.name][0]
                    for a in R.atoms.values(): # find all reactive atoms in this reactant
                        aa=deepcopy(a)
                        if aa['reactant']==ridx:  # this is an atom in this molecule
                            aa['reactant']=R.reactants[ridx]
                            if not aa in atoms:
                                atoms.append(aa)
            logging.debug(f'atoms in reactions for which {self.name} is a reactant:')
            logging.debug(f'{atoms}')
            for a in atoms:
                z=a['z']  # this is the z value listed in the cfg file
                atomName=a['atom']
                resNum=a['resid']
                resName=a['reactant']
                if mdict[resName].generator:
                    logging.debug(f'{self.name} is not a precursor -- it is generated by {mdict[resName].generator.name}')
                assert resName==self.name
                Aidx=self.Coords.get_atom_attribute('globalIdx',{'atomName':atomName,'resName':resName,'resNum':resNum})
                self.Coords.set_atom_attribute('z',z,{'globalIdx':Aidx})
        # logging.debug(f'after propagate_z on {self.name}:\n{self.Coords.A.to_string()}')

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
            assert type(R)==Reaction,'HTPolyNet only recognizes Reaction-type generators at the moment'
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2.')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            can_react=all([a in available_molecules for a in R.reactants.values()])
            if not can_react:
                raise Exception(f'Cannot generate {self.name} because required reactants have not been generated')

            # local copies of all reactant molecules
            reactants={}
            for n,r in R.reactants.items():
                reactants[n]=deepcopy(available_molecules[r])
                ts=reactants[n].sequence
                logging.debug(f'{n} sequence: {reactants[n].sequence}')
                for rn in ts:
                    self.sequence.append(rn)
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
                ban=mB.Coords.get_atom_attribute('atomName',{'globalIdx':Bidx})
                brn=mB.Coords.get_atom_attribute('resNum',{'globalIdx':Bidx})
                logging.info(f'Due to post-bonding reindexing, the two bonded atoms are now {mA.name}-R{arn}-{aan}({Aidx}) and {mB.name}-R{brn}-{ban}({Bidx})')
                if len(bases)==0 or mA not in bases:
                    bases.append(mA)
            assert len(bases)==1,f'Error: Reaction {R.name} results in more than one molecular fragment product'
            base=bases[0]
            self.merge(base)
            self.write_mol2(filename=f'{self.name}.mol2')
        else:
            logging.info(f'Using input molecules/inputs/{self.name}.mol2 as a generator.')
            pfs.checkout(f'molecules/inputs/{self.name}.mol2')
            self.sequence.append(self.name)

        self.parameterize(outname,**kwargs)
        # identify new transferable topological units (angles, dihedrals)
        self.minimize(outname,**kwargs)

    def label_ring_atoms(self,cycles):
        adf=self.Coords.A
        self.set_atomset_attribute('cycle-idx',np.zeros(adf.shape[0]).astype(int))
        cidx=1
        for l,cl in cycles.items():
            for c in cl:
                for idx in c:
                    self.set_atom_attribute('cycle-idx',cidx,{'globalIdx':idx})
                cidx+=1
        # logging.debug(f'label_ring_atoms for {self.name}:\n{adf.to_string()}')

    def get_resname(self,internal_resid):
        logging.debug(f'{self.name} sequence: {self.sequence}')
        return self.sequence[internal_resid-1]

    def inherit_sea_from_reactants(self,molecules,sea_list):
        if self.name in sea_list:
            logging.debug(f'No need to inherit sea for {self.name}')
            return
        adf=self.Coords.A
        # logging.debug(f'Inherit sea for {self.name}; Atoms data frame initially:\n{adf.to_string()}')
        self.set_atomset_attribute('sea-idx',np.arange(adf.shape[0]))
        donors=get_base_reactants(self.name,molecules)
        logging.debug(f'Inherit sea: {self.name} base reactants {", ".join([d.name for d in donors])}')
        seaidx_shift=0
        for D in donors:
            if not D.name in sea_list:
                continue
            logging.debug(f' from Donor {D.name} seaidx_shift {seaidx_shift}:')
            # logging.debug(f'Here is the whole {D.name}:\n{D.Coords.A.to_string()}')
            recv=adf[adf['resName']==D.name]
            resids=list(set(recv['resNum']))
            logging.debug(f'  resids in {self.name} of {D.name} are {resids}')
            for rid in resids:
                logging.debug(f'    resid {rid}')
                recvr=recv[recv['resNum']==rid]
                for i,r in recvr.iterrows():
                    aidx=r['globalIdx']
                    donor_seaidx=D.Coords.get_atom_attribute('sea-idx',{'atomName':r['atomName']})
                    self.Coords.set_atom_attribute('sea-idx',donor_seaidx+seaidx_shift,{'globalIdx':aidx})
                logging.debug(f'     adding {max(recvr["sea-idx"])} to seaidx_shift {seaidx_shift}')
                seaidx_shift+=max(D.Coords.A['sea-idx'])
        # logging.debug(f'Molecule {self.name} after inheriting sea:\n'+self.Coords.A.to_string())

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
        self.Topology.add_bonds(pairs,enumerate_others=False)

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
        # myz=myC.get_atom_attribute('z',{'globalIdx':myidx})
        # otz=otC.get_atom_attribute('z',{'globalIdx':otidx})
        # assert myz>0 and otz>0,f'No bond permissible: z({myidx}):{myz}, z({otidx}):{otz}'

        mypartners=myT.bondlist.partners_of(myidx)
        otpartners=otT.bondlist.partners_of(otidx)
        logging.info(f'Partners of {myidx} {mypartners}')
        logging.info(f'Partners of {otidx} {otpartners}')
        myHpartners={k:v for k,v in zip(mypartners,[myA[myA['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners={k:v for k,v in zip(otpartners,[otA[otA['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        assert len(myHpartners)>0,f'Error: atom {myidx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {otidx} does not have a deletable H atom!'
            
        Ri=myC.get_R(myidx)
        Rj=otC.get_R(otidx)
        minHH=(1.e9,(-1,'x'),(-1,'x'))
        if self!=other:
            overall_maximum=(-1.e9,-1,-1)
        ''' Identify the H atom on each atom in be pair to
            be bonded that result in the "optimum" choice
            for deletion or, if this is intermolecular,
            the "optimum" choice for the location/orientation
            of the other. '''
        totc={}
        for myH,myHnm in myHpartners.items():  # keys are globalIdx's, values are names
            totc[myH]={}
            Rh=myC.get_R(myH)
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            for otH,otHnm in otHpartners.items():
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
        idxshift=0
        if self!=other:
            minD,myH,otH=overall_maximum
            myHnm=myHpartners[myH]
            otHnm=otHpartners[otH]
            shifts=myC.merge(totc[myH][otH])
            idxshift=shifts[0]
            otidx+=idxshift
        else:
            mhh,myH,otH=minHH
            myHnm=myHpartners[myH]
            otHnm=otHpartners[otH]
            logging.info(f'Executing intramolecular reaction  in {self.name}')
            logging.info(f'Two Hs closest to each other are {myH}({myHnm}) and {otH}({otHnm})')

        if self!=other:
            self.Topology.merge(other.Topology)

        otH+=idxshift
        otHpartners={(k+idxshift):v for k,v in otHpartners.items()}
    
        # change H atom names so that it looks like the highest-name-value one was deleted
        logging.debug(f'myHpartners: {myHpartners}')
        logging.debug(f'otHpartners: {otHpartners}')
        myavails=list(sorted(myHpartners.values()))[:-1]
        otavails=list(sorted(otHpartners.values()))[:-1]
        del myHpartners[myH]
        del otHpartners[otH]
        T=self.Topology.D['atoms']
        C=self.Coords.A
        for h in myHpartners:
            myHpartners[h]=myavails.pop(0)
            T.iloc[h-1,T.columns=='atom']=myHpartners[h]
            C.iloc[h-1,C.columns=='atomName']=myHpartners[h]
        for h in otHpartners:
            otHpartners[h]=otavails.pop(0)
            T.iloc[h-1,T.columns=='atom']=otHpartners[h]
            C.iloc[h-1,C.columns=='atomName']=otHpartners[h]
        # this makes sure that it always looks like the same atom was deleted
        self.add_bonds(pairs=[(myidx,otidx)])
        # myC.set_atom_attribute('z',myz-1,{'globalIdx':myidx})
        # otC.set_atom_attribute('z',otz-1,{'globalIdx':otidx})
        new_idx=self.delete_atoms(idx=[myH,otH],return_idx_of=return_idx_of)


        self.Topology.bond_source_check()
        return new_idx
