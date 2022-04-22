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
        self.Coords['gro'].set_attribute('sea-idx',analyze_sea(f'{n}-sea'))
        self.Coords['gro'].write_sea(f'{n}.sea')

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
        # TODO: sync mol2 coords

    def read_sea(self,filename):
        assert os.path.exists(filename),f'SEA file {filename} not found.'
        self.Coords['gro'].read_sea(filename)

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
                    if not 'z' in m.Coords['gro'].D['atoms']:
                        logging.info(f'Initializing z attribute of all atoms in {m.name}')
                        m.Coords['gro'].set_attribute('z',np.zeros(m.Coords['gro'].D['atoms'].shape[0]))
                    idx=m.Coords['gro'].get_idx({'atomName':a['atom'],'resNum':a['resid']})-1
                    logging.info(f'Modifing z attribute of molecule {m.name} atom {idx+1}')
                    m.Coords['gro'].set_atom_attribute('z',a['z'],{'globalIdx':idx})#   D['atoms']['z'].iloc[idx]=a.z

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

        logging.debug(f'new_bond {myidx}:{self.Coords["mol2"].get_atom_attribute("atomName",{"globalIdx":myidx})} - {otidx}:{other.Coords["mol2"].get_atom_attribute("atomName",{"globalIdx":otidx})}')

        mypartners=myC.bondlist.partners_of(myidx)
        otpartners=otC.bondlist.partners_of(otidx)
        myHpartners=[k for k,v in zip(mypartners,[myA[myA['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')]
        otHpartners=[k for k,v in zip(otpartners,[otA[otA['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')]
        assert len(myHpartners)>0,f'Error: atom {myidx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {otidx} does not have a deletable H atom!'
        logging.debug(f'Hs on {myidx}: '+','.join([f'{h}:{self.Coords["mol2"].get_atom_attribute("atomName",{"globalIdx":h})}' for h in myHpartners]))
        logging.debug(f'Hs on {otidx}: '+','.join([f'{h}:{other.Coords["mol2"].get_atom_attribute("atomName",{"globalIdx":h})}' for h in otHpartners]))
            
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
            logging.debug(f'myH {myH} Rh {Rh}')
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            logging.debug(f'    Rih {Rih}')
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
                    logging.debug(f'   R {R}')
                    # rotate translate all donor atoms!
                    totc[myH][otH].rotate(R)
                    #Rj=totc[otH].get_R(otidx)
                    Rk=totc[myH][otH].get_R(otH)
                    #Rkj=Rk-Rj
                    #Rkj*=1.0/np.linalg.norm(Rkj)
                    #Rhk=Rh-Rk
                    #rhk=np.linalg.norm(Rhk)
                    # overlap the other H atom with self's 
                    # reactive atom by translation
                    Rik=Ri-Rk
                    totc[myH][otH].translate(Rik)
                    totc[myH][otH].write_mol2(f'tmp-{myH}-{otH}.mol2')
                    #Rj=totc[otH].get_R(otidx)
                    #Rk=totc[otH].get_R(otH)
                    #Rkj=Rk-Rj
                    #Rkj*=1.0/np.linalg.norm(Rkj)
                    #Rhk=Rh-Rk
                    #rhk=np.linalg.norm(Rhk)
                    minD=myC.minimum_distance(totc[myH][otH],self_excludes=[myH],other_excludes=[otH])
                    if minD>overall_maximum[0]:
                        overall_maximum=(minD,myH,otH)
                else:
                    if rhk<minHH[0]:
                        minHH=(rhk,myH,otH)
        if self!=other:
            minD,myH,otH=overall_maximum
            shifts=myC.merge(totc[myH][otH])
            myC.write_mol2('tmp.mol2')
            idxshift=shifts[0]
            otidx+=idxshift
            otH+=idxshift
        else:
            mhh,myH,otH=minHH
            logging.info(f'Executing intramolecular reaction  in {self.name}')
            logging.info(f'Two Hs closest to each other are {myH} and {otH}')
    
        myC.add_bonds(pairs=[(myidx,otidx)])
        myC.delete_atoms(idx=[myH,otH])


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