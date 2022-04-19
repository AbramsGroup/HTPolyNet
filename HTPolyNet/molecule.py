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

class Atom:
    def __init__(self,datadict):
        self.monomer=datadict['monomer']
        self.name=datadict['atom']
        # z: maximum number of bonds this atom will form
        # in a crosslinking reaction
        self.z=int(datadict.get('z',0))
    def __str__(self):
        return f'{self.monomer}_{self.name}({self.z})'
    def to_yaml(self):
        return r'{'+f'monomer: {self.monomer}, atom: {self.name}, z: {self.z}'+r'}'


# class CappingBond:
#     boc={1:'-',2:'=',3:'≡'}
#     def __init__(self,jsondict):
#         self.pairnames=jsondict["pair"]
#         self.bondorder=jsondict.get("order",1)
#         self.deletes=jsondict.get("deletes",[])
#     def to_yaml(self):
#         return r'{'+f'pair: {self.pairnames}, order: {self.bondorder}, deletes: {self.deletes}'+r'}'
#     def __str__(self):
#         s=self.pairnames[0]+CappingBond.boc[self.bondorder]+self.pairnames[1]
#         if len(self.deletes)>0:
#             s+=' D['+','.join(self.deletes)+']'
#         return s


class Monomer:
    def __init__(self,name='',use_sea=False):
        self.name=name
        self.use_sea=use_sea
        self.reactive_atoms=[]
        self.Topology=None
        self.Coords={}  # 'gro' and 'mol2' are the two keys
        self.cstale=''  # ['','gro','mol2']

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

    def read_sea(self,filename):
        assert os.path.exists(filename),f'SEA file {filename} not found.'
        self.Coords['gro'].read_sea(filename)

    def calculate_sea(self):
        ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
        n=self.name
        boxsize=np.array(self.Coords['gro'].maxspan())+2*np.ones(3)
        pfs.checkout('nvt-sea.mdp')
        for ex in ['top','itp','gro']:
            pfs.checkout(f'{n}-p.{ex}')
        logging.info(f'hot md running...output to {n}-p-sea')
        grompp_and_mdrun(gro=f'{n}-p',top=f'{n}-p',
                        mdp='nvt-sea',out=f'{n}-p-sea',boxSize=boxsize)
        self.Coords['gro'].D['atoms']['sea-idx']=analyze_sea(f'{n}-p-sea')
        self.Coords['gro'].write_sea(f'{n}-p.sea')

    def parameterize(self,outname='',**kwargs):
        pfs.checkout(f'{self.name}.mol2')
        assert os.path.exists(f'{self.name}.mol2'),f'Cannot parameterize monomer {self.name} without {self.name}.mol2'
        if outname=='':
            outname=f'{self.name}-p'
        GAFFParameterize(self.name,outname,**kwargs)
        self.read_topology(f'{outname}.top')
        self.read_coords(f'{outname}.gro')
        self.read_coords(f'{outname}.mol2')
        assert self.cstale=='',f'Error: {self.cstale} coords are stale'

    def add_reactive_atom(self,ra):
        # add this reactive atom if its z is greater than any z
        # for an ra with the same name
        assert type(ra)==Atom,f'Argument to add_reactive_atom must be type Atom not {type(ra)}'
        found=False
        #logging.info(f'add_reactive_atom at monomer {self.name}({id(self)}), currently with {len(self.reactive_atoms)} reactive atoms.')
        for i in range(len(self.reactive_atoms)):
            if self.reactive_atoms[i].name==ra.name:
                #logging.info(f'   found {self.reactive_atoms[i].name} at {i}; z {self.reactive_atoms[i].z}')
                if ra.z>self.reactive_atoms[i].z:
                    logging.info(f'Overwriting zmax of {ra.name}')
                    self.reactive_atoms[i].z=ra.z
                found=True
        if not found:
            logging.info(f'Adding reactive atom {ra.name}(zmax={ra.z}) to monomer {self.name}')
            self.reactive_atoms.append(ra)

    def update_atom_specs(self,oldmol2):
        ''' Atom specifications from the configuration file refer to atom names in the
            user-provided mol2 files.  After processing via ambertools, the output mol2
            files have renamed atoms (potentially) in the same order as the atoms in the
            user-provided mol2.  This method updates the user-provided atom specifications
            so that they reflect the atom names generated by ambertools. This should not
            be necessary if antechamber obeys the atom ordering in the input mol2 file. '''
        oa=oldmol2.D['atoms']
        na=self.Coords['mol2'].D['atoms']
        for ra in self.reactive_atoms:
            assert type(ra)==Atom
            n=ra.name
            idx=oa[oa['atomName']==n]['globalIdx'].values[0]
            ra.name=na[na['globalIdx']==idx]['atomName'].values[0]
            if n!=ra.name:
                logging.info('Good thing you did this: the atom names changed!')

        # for c in self.capping_bonds:
        #     p=c.pairnames
        #     newpairnames=[]
        #     for n in p:
        #         idx=oa[oa['atomName']==n]['globalIdx'].values[0]
        #         nn=na[na['globalIdx']==idx]['atomName'].values[0]
        #         newpairnames.append(nn)
        #     c.pairnames=newpairnames

    def __str__(self):
        s=self.name+'\n'
        for r,a in self.reactive_atoms.items():
            s+=f'   reactive atom: {r}:'+str(a)+'\n'
        # for c in self.capping_bonds:
        #     s+='   cap: '+str(c)+'\n'
        return s

class Molecule:
    def __init__(self,name=''):
        self.name=name
        self.Topology=None
        self.Coords={}

    def __str__(self):
        restr=f'{self.name} '
        if not self.Topology:
            restr+=f'(empty topology) '
        if len(self.Coords)==0:
            restr+=f'(empty coords) '
        return restr+'\n'

    def populate_from_monomer(self,monomer=None):
        if monomer:
            self.name=monomer.name
            self.Topology=monomer.Topology
            self.Coords=monomer.Coords

    @classmethod
    def generate_from_monomer(cls,monomer=None):
        inst=cls('empty')
        if monomer:
            inst=cls()
            inst.populate_from_monomer(monomer)
        return inst

    @classmethod
    def generate_from_reaction(cls,R,moldict):
        assert type(moldict)==dict
        inst=cls(R.template['product'])
        inst.merge(deepcopy(moldict[R.template['reactants'][0]]))
        for name in R.template['reactants'][1:]:
            inst.merge(deepcopy(moldict[name]))
        # TODO -- get atoms and make new bond
        return inst

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

    def parameterize(self):
        pass
        '''
        1. write mol2
        2. Gaff-parameterize -> gro/top
        3. minimize->new gro
        4. read top->Topology, read gro->Coord
        '''

    def read_topology(self,filename):
        assert os.path.exists(filename),f'Topology file {filename} not found.'
        if self.Topology:
            logging.warning(f'overwriting topology of molecule {self.name} from file {filename}')
        self.Topology=Topology.read_gro(filename)

    def read_coords(self,filename):
        assert os.path.exists(filename),f'Coordinate file {filename} not found.'
        if self.Coords:
            logging.warning(f'overwriting coordinates of molecule {self.name} from file {filename}')
        basename,ext=os.path.splitext(filename)
        if ext=='mol2':
            self.Coords=Coordinates.read_mol2(filename)
        elif ext=='gro':
            self.Coords=Coordinates.read_gro(filename)
        else:
            raise Exception(f'Coordinate filename extension {ext} is not recognized.')

    def update_topology(self,t):
        self.Topology.merge(t)

    def update_coords(self,c):
        self.Coords.copy_coords(c)

class Reaction:
    ''' reactions may only be initialized after monomers '''
    def __init__(self,jsondict):
        self.jsondict=jsondict
        ''' attributes expected in cfg file '''
        self.name=jsondict.get('name','unnamed reaction')
        self.atoms=jsondict.get('atoms','empty reaction')
        self.bonds=jsondict.get('bonds','empty reaction')
        self.template=jsondict.get('template','empty reaction')
        self.restrictions=jsondict.get('restrictions',None)
        self.stage=jsondict.get('stage','cure')

        ''' attributes set up by other methods'''
        self.molecules={}

    def pass_reactive_atoms(self,monomer_dict):
        ''' a reaction identifies reactive atoms in a monomer, so we must
            parse the reaction to add this reactive atom to the monomer '''
        for i,at in self.atoms.items():
            # 'monomer' entry in atoms dict is the name of the monomer to which this atom belongs
            assert at['monomer'] in monomer_dict,f'Reaction {self.name} references non-existent monomer {at["monomer"]}'
            logging.info(f'Reaction "{self.name}" seeks to add reactive atom {at["atom"]} to monomer {at["monomer"]}')
            monomer_dict[at['monomer']].add_reactive_atom(Atom(at))

    def __str__(self):
        return f'Reaction "{self.name}"'

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