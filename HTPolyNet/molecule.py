from itertools import combinations_with_replacement, product
import os
from copy import deepcopy
from re import A
import pandas as pd
import numpy as np
import logging

from HTPolyNet.topocoord import TopoCoord

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
    def __init__(self,jsondict):
        self.jsondict=jsondict
        self.name=jsondict.get('name','')
        self.atoms=jsondict.get('atoms',{})
        self.bonds=jsondict.get('bonds',{})
        self.reactants=jsondict.get('reactants',{})
        self.product=jsondict.get('product','')
        self.restrictions=jsondict.get('restrictions',{})
        self.stage=jsondict.get('stage','')
        self.probability=1.0

    def __str__(self):
        retstr=f'Reaction "{self.name}"\n'
        for i,r in self.reactants.items():
            retstr+=f'reactant {i}: {r}\n'
        retstr+=f'product {self.product}\n'
        for i,a in self.atoms.items():
            retstr+=f'atom {i}: {a}\n'
        for b in self.bonds:
            retstr+=f'bond {b}\n'
        return retstr

    def get_bond_atom_globalIdx(self,bonddict,mol,moldict):
        A,B=[self.atoms[i] for i in bonddict['atoms']]
        # logging.debug(f'get_bond_atom_globalIdx: mol {mol.name}')
        # logging.debug(f'  -> bond\n    {self.reactants}\n    {bonddict}:\n    A {A}\n    B {B}')
        Aidx,Aresid=self.get_atom_globalIndx(A,mol,moldict)
        Bidx,Bresid=self.get_atom_globalIndx(B,mol,moldict)
        # logging.debug(f'     A: {Aidx}  B: {Bidx}')
        return (Aidx,Bidx),(Aresid,Bresid),(A['atom'],B['atom'])

    def get_atom_globalIndx(self,A,mol,moldict):
        aname=A['atom']
        aN=A['reactant']
        areactantname=self.reactants[aN]
        aresid_inreactant=A['resid']
        # logging.debug(f'determining globalIdx of aname {aname} in reactant {aN}({areactantname}):{aresid_inreactant}')
        aresid=0
        for k,v in self.reactants.items():
            if v==areactantname:
                aresid+=aresid_inreactant
                break
            else:
                # logging.debug(f'preseq {v} {moldict[v].sequence}')
                aresid+=len(moldict[v].sequence)
        # logging.debug(f'reaction {self.name}: asking for globalIdx of resNum {aresid} atomName {aname}')
        idx=mol.TopoCoord.get_gro_attribute_by_attributes('globalIdx',{'resNum':aresid,'atomName':aname})
        return idx,aresid

    def get_raz(self,moldict={}):  # [resname][atomname]=[list of z-values detected in this reaction]
        reactantdict=self.reactants
        razdict={}
        for ltr,atomrec in self.atoms.items():
            reactant_name=reactantdict[atomrec['reactant']]
            resid=atomrec['resid']
            an=atomrec['atom']
            z=atomrec['z']
            if len(moldict)==0:
                resname=reactant_name
            else:
                # logging.debug(f'grabbing resname of {resid} from position {resid-1} of {reactant_name} {moldict[reactant_name].sequence}')
                resname=moldict[reactant_name].sequence[resid-1]
            if not resname in razdict:
                razdict[resname]={}
            if not an in razdict[resname]:
                razdict[resname][an]=[]
            razdict[resname][an].append(z)
        return razdict

class Molecule:
    def __init__(self,name='',generator=None):
        self.name=name
        self.TopoCoord=TopoCoord()
        self.generator=generator
        self.sequence=[]
        self.origin=None
        self.reaction_bonds=[]
        self.symmetry_relateds=[]

    # def __str__(self):
    #     restr=f'{self.name} '
    #     if not self.Topology:
    #         restr+=f'(empty topology) '
    #     if len(self.Coords)==0:
    #         restr+=f'(empty coords) '
    #     return restr+'\n'

    def set_origin(self,value):
        self.origin=value

    def get_origin(self):
        return self.origin

    # def num_atoms(self):
    #     if hasattr(self.Coords,'N'):
    #         return self.Coords.N
    #     else:
    #         return self.Coords.shape[0]

    def previously_parameterized(self):
        rval=True
        for ext in ['mol2','top','itp','gro']:
            rval=rval and pfs.exists(os.path.join('molecules/parameterized',f'{self.name}.{ext}'))
        return rval

    def parameterize(self,outname='',**kwargs):
        assert os.path.exists(f'{self.name}.mol2'),f'Cannot parameterize molecule {self.name} without {self.name}.mol2 as input'
        if outname=='':
            outname=f'{self.name}'
        GAFFParameterize(self.name,outname,**kwargs)
        self.load_top_gro(f'{outname}.top',f'{outname}.gro',mol2filename=f'{outname}.mol2')
        #assert self.cstale=='',f'Error: {self.cstale} coords are stale'

    def calculate_sea(self):
        ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
        n=self.name
        boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
        pfs.checkout('mdp/nvt-sea.mdp')
        for ex in ['top','itp','gro']:
            pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        logging.info(f'Hot md running...output to {n}-sea')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='nvt-sea',out=f'{n}-sea',boxSize=boxsize)
        sea_srs=analyze_sea(f'{n}-sea')
        self.set_gro_attribute('sea-idx',sea_srs)

    def minimize(self,outname='',**kwargs):
        if outname=='':
            outname=f'{self.name}'
        n=self.name
        boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
        pfs.checkout('mdp/em-single-molecule.mdp')
        if 'checkout_required' in kwargs:
            for ex in ['top','itp','gro']:
                pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp='em-single-molecule',out=f'{outname}',boxSize=boxsize)
        self.TopoCoord.read_gro(f'{n}.gro')

#     def propagate_z(self,reactions,mdict):
#         ''' Assign all z values in Coords.A dataframes for each molecule '''
#         logging.debug(f'propagate_z called for {self.name}')
#         adf=self.TopoCoord.gro_DataFrame('atoms')
#         self.set_gro_attribute('z',0)
#         if self.generator:
#             R=self.generator
#             logging.debug(f'{self.name} was generated by {R.name}')
#             resids=list(set(adf['resNum']))
#             resnames=[]
#             logging.debug(f'{self.name} has resids {resids}')
#             for resid in resids:
#                 resnames.append(list(set(adf[adf['residue']==resid]['residue']))[0])
#             logging.debug(f'{self.name} has resnames {resnames}')
#             for resid,resname in zip(resids,resnames):
#                 logging.debug(f'Product {self.name} taking z\'s from residue {resname}')
#                 mol=mdict[resname]
#                 radf=mol.TopoCoord.gro_DataFrame('atoms')
#                 zs=radf[radf['z']>0]
#                 for i,r in zs.iterrows():
#                     atomName=r['atomName']
#                     resNum=r['resNum']
#                     resName=r['resName']
#                     z=r['z']
#                     self.TopoCoord.set_gro_attribute_by_attributes('z',z,{'atomName':atomName,'resName':resName})
#             ''' correct from self's bondlist '''
# #            bondlist=self.TopoCoord.Topology.bondlist
#             reactive_atoms_idx=list(adf[adf['z']>0]['globalIdx'])
#             # logging.debug(f'reactive_atoms_idx {reactive_atoms_idx}')
#             for i in range(len(reactive_atoms_idx)):
#                 ix=reactive_atoms_idx[i]
#                 for j in range(i,len(reactive_atoms_idx)):
#                     jx=reactive_atoms_idx[j]
#                     iz=self.TopoCoord.get_gro_attribute_by_attributes('z',{'globalIdx':ix})
#                     jz=self.TopoCoord.get_gro_attribute_by_attributes('z',{'globalIdx':jx})
#                     if self.TopoCoord.are_bonded(ix,jx):
#                         self.TopoCoord.set_gro_attribute_by_attributes('z',iz-1,{'globalIdx':ix})
#                         self.TopoCoord.set_gro_attribute_by_attributes('z',jz-1,{'globalIdx':jx})
#             # logging.debug(f'after propagate_z product {self.name}:\n{self.Coords.A.to_string()}')
#         else:  # this molecule is not a product -- was generated using an input mol2
#             ''' find any reaction in which this is a reactant as long as it is not a product '''
#             atoms=[]
#             for R in reactions:
#                 if self.name in R.reactants.values():
#                     ridx=[k for k,v in R.reactants.items() if v==self.name][0]
#                     for a in R.atoms.values(): # find all reactive atoms in this reactant
#                         aa=deepcopy(a)
#                         if aa['reactant']==ridx:  # this is an atom in this molecule
#                             aa['reactant']=R.reactants[ridx]
#                             if not aa in atoms:
#                                 atoms.append(aa)
#             logging.debug(f'atoms in reactions for which {self.name} is a reactant:')
#             logging.debug(f'{atoms}')
#             for a in atoms:
#                 z=a['z']  # this is the z value listed in the cfg file
#                 atomName=a['atom']
#                 resNum=a['resid']
#                 resName=a['reactant']
#                 if mdict[resName].generator:
#                     logging.debug(f'{self.name} is not a precursor -- it is generated by {mdict[resName].generator.name}')
#                 assert resName==self.name
#                 Aidx=self.TopoCoord.get_gro_attribute_by_attributes('globalIdx',{'atomName':atomName,'resName':resName,'resNum':resNum})
#                 self.TopoCoord.set_gro_attribute_by_attributes('z',z,{'globalIdx':Aidx})
#         # logging.debug(f'after propagate_z on {self.name}:\n{self.Coords.A.to_string()}')

    def generate(self,outname='',available_molecules={},**kwargs):
        logging.info(f'Generating {self.name}.mol2 for parameterization')
        if outname=='':
            outname=f'{self.name}'
        if self.generator:
            R=self.generator
            assert type(R)==Reaction,'HTPolyNet only recognizes Reaction-type generators at the moment'
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2.')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            can_react=all([a in available_molecules for a in R.reactants.values()])
            if not can_react:
                raise Exception(f'Cannot generate {self.name} because required reactants have not been generated')

            composite_mol=Molecule()
            shifts=[(0,0,0)]  # atom, bond, resid
            for n,ri in R.reactants.items():
                # logging.debug(f'adding {available_molecules[ri].name} to composite:\n{available_molecules[ri].TopoCoord.Coordinates.A.to_string()}')
                shifts.append(composite_mol.merge(deepcopy(available_molecules[ri])))
            self.TopoCoord=deepcopy(composite_mol.TopoCoord)
            self.set_sequence()
            self.set_reaction_bonds(available_molecules=available_molecules)
            reactantName=R.product
            # logging.debug(f'Generation of {self.name}: composite molecule has {len(self.sequence)} resids')
            # logging.debug(f'generation of {self.name}: composite molecule:\n{composite_mol.TopoCoord.Coordinates.A.to_string()}')
            idx_mapper=self.make_bonds()
            # nrb=[]
            # for b in self.reaction_bonds:
            #     at,rr,nn=b
            #     i,j=at
            #     i=idx_mapper[i]
            #     j=idx_mapper[j]
            #     nrb.append(((i,j),rr,nn))
            self.TopoCoord.write_mol2(filename=f'{self.name}.mol2',molname=self.name)
        else:
            logging.info(f'Using input molecules/inputs/{self.name}.mol2 as a generator.')
            pfs.checkout(f'molecules/inputs/{self.name}.mol2')
            # self.set_reaction_bonds(available_molecules=available_molecules)
            reactantName=self.name
            # self.sequence.append(self.name)
        self.parameterize(outname,**kwargs)
        self.minimize(outname,**kwargs)
        self.set_sequence()
        self.TopoCoord.set_gro_attribute('reactantName',reactantName)

    def set_reaction_bonds(self,available_molecules={}):
        R=self.generator
        self.reaction_bonds=[]
        if R:
            for bond in R.bonds:
                (Aidx,Bidx),(aresid,bresid),(Aname,Bname)=R.get_bond_atom_globalIdx(bond,self,available_molecules)
                self.reaction_bonds.append(((Aidx,Bidx),(aresid,bresid),(Aname,Bname)))
            logging.debug(f'{R.name} reaction_bonds\n{self.reaction_bonds}')

    def set_sequence(self):
        adf=self.TopoCoord.gro_DataFrame('atoms')
        self.sequence=[]
        current_resid=0
        for i,r in adf.iterrows():
            ri=r['resNum']
            rn=r['resName']
            if ri!=current_resid:
                current_resid=ri
                self.sequence.append(rn)
        # logging.debug(f'{self.name} sequence: {self.sequence}')

    def idx_mappers(self,otherTC,other_bond):
        seq_res_is_bystander=[False for _ in self.sequence]
        # logging.debug(f'idx_mappers begins: template name {self.name}')
        # get system resName, resNum, and atomNames for each of the two bonded atoms
        i_idx,j_idx=other_bond
        i_resName,i_resNum,i_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i_idx})
        j_resName,j_resNum,j_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j_idx})
        # Either of atom i or j *could* bond to a third or even fourth residue that
        # is not the residue of the other.  Any such "bystander" residue is asserted
        # to be represented in the template molecule, and must therefore be included
        # along with the residues of i and j among the residues from which interactions
        # can be mapped.
        neighbors_of_i=otherTC.partners_of(i_idx)
        resid_bystanders_of_ij=[]  # bystanders bonded to atom i that are not in resid of atom j
        # logging.debug(f'neighbors of i {i_idx} {i_resName} {i_resNum} {i_atomName}:')
        for xx in neighbors_of_i:
            x_resName,x_resNum,x_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xx})
            if x_resNum!=i_resNum and x_resNum!=j_resNum:
                resid_bystanders_of_ij.append((x_resNum,x_resName))
                # logging.debug(f'{xx} {x_resName} {x_resNum} {x_atomName}')
        neighbors_of_j=otherTC.partners_of(j_idx)
        # logging.debug(f'neighbors of j {j_idx} {j_resName} {j_resNum} {j_atomName}:')
        resid_bystanders_of_ji=[] # bystanders bonded to atom j that are not in resid of atom i
        for xx in neighbors_of_j:
            x_resName,x_resNum,x_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xx})
            if x_resNum!=i_resNum and x_resNum!=j_resNum:
                resid_bystanders_of_ji.append((x_resNum,x_resName))
                # logging.debug(f'{xx} {x_resName} {x_resNum} {x_atomName}')

        resid_bystanders=[*resid_bystanders_of_ij,*resid_bystanders_of_ji]
        # logging.debug(f'idx_mappers: other_bond {other_bond} resid_bystanders {resid_bystanders}')
        temp_bystanders=[]
        inst_bystanders=[]
        for xx in resid_bystanders:
            rn,rname=xx
            if not rname in self.sequence:
                logging.error(f'secondary neighbor {rn} {rname}: no res in pattern sequence {self.sequence}')
                raise Exception('this is a bug')
            for i,rnm in enumerate(self.sequence):
                ib=seq_res_is_bystander[i]
                if rnm==rname and not ib:
                    seq_res_is_bystander[i]=True
                    temp_bystanders.append(i+1) # resids start at 1 not 0!!!
                    inst_bystanders.append(rn)
                    break
            else:
                logging.error(f'secondary neighbor {rn} {rname}: no available res in pattern sequence {self.sequence} {seq_res_is_bystander}')
                    
        temp2inst={}
        inst2temp={}
        temp_iresid=-1
        temp_jresid=-1
        # identify the template bond represented by the other_bond parameter
        for b in self.reaction_bonds:
            (Aidx,Bidx),(aresid,bresid),(Aname,Bname)=b
            Aresname=self.sequence[aresid-1]
            Bresname=self.sequence[bresid-1]
            # logging.debug(f'idx_mappers: {Aresname} {aresid} {Bresname} {bresid}')
            if (i_atomName,i_resName)==(Aname,Aresname):
                temp_iresid=aresid
                temp_jresid=bresid
                break # found it -- stop looking
            elif (i_atomName,i_resName)==(Bname,Bresname):
                temp_iresid=bresid
                temp_jresid=aresid
                break
        if temp_iresid==-1:
            logging.error(f'Mappers using template {self.name} unable to map from instance bond {i_resName}-{i_resNum}-{i_atomName}---{j_resName}-{j_resNum}-{j_atomName}')
            raise Exception

        # use dataframe merges to create globalIdx maps
        instdf=otherTC.Coordinates.A
        tempdf=self.TopoCoord.Coordinates.A
        inst2temp={}
        temp2inst={}
        for inst,temp in zip([i_resNum,j_resNum,*inst_bystanders],[temp_iresid,temp_jresid,*temp_bystanders]):
            idf=instdf[instdf['resNum']==inst][['globalIdx','atomName']].copy()
            # logging.debug(f'idf res {inst}:\n{idf.to_string()}')
            tdf=tempdf[tempdf['resNum']==temp][['globalIdx','atomName']].copy()
            # logging.debug(f'tdf res {temp}:\n{tdf.to_string()}')
            tdf=tdf.merge(idf,on='atomName',how='inner',suffixes=('_template','_instance'))
            # logging.debug(f'merged\n{tdf.to_string()}')
            for i,r in tdf.iterrows():
                temp_idx=r['globalIdx_template']
                inst_idx=r['globalIdx_instance']
                # logging.debug(f't {temp_idx} <-> i {inst_idx}')
                inst2temp[inst_idx]=temp_idx
                assert not temp_idx in temp2inst,f'Error: temp_idx {temp_idx} already claimed in temp2inst; bug'
                temp2inst[temp_idx]=inst_idx
        assert len(inst2temp)==len(temp2inst),f'Error: could not establish two-way dict of atom globalIdx'
        return (inst2temp,temp2inst)

    def get_angles_dihedrals(self,bond):
        ai,aj=bond
        d=self.TopoCoord.Topology.D['angles']
        ad=d[((d.ai==ai)&(d.aj==aj))|
             ((d.ai==aj)&(d.aj==ai))|
             ((d.aj==ai)&(d.ak==aj))|
             ((d.aj==aj)&(d.ak==ai))].copy()
        d=self.TopoCoord.Topology.D['dihedrals']
        td=d[((d.ai==ai)&(d.aj==aj))|
             ((d.ai==aj)&(d.aj==ai))|
             ((d.aj==ai)&(d.ak==aj))|
             ((d.aj==aj)&(d.ak==ai))|
             ((d.ak==ai)&(d.al==aj))|
             ((d.ak==aj)&(d.al==ai))].copy()
        check=True
        for a in ['ai','aj','ak','al']:
            check=check and td[a].isnull().values.any()
        if check:
            logging.error('NAN in molecule/dihedrals')
            raise Exception

        d=self.TopoCoord.Topology.D['pairs']
        paird=pd.DataFrame()
        for ai,al in zip(td.ai,td.al):
            tpair=d[((d.ai==ai)&(d.aj==al))|
                    ((d.ai==al)&(d.aj==ai))].copy()
            paird=pd.concat((paird,tpair),ignore_index=True)
        check=True
        for a in ['ai','aj']:
            check=check and paird[a].isnull().values.any()
        if check:
            logging.error('NAN in molecule/pairs')
            raise Exception
        return ad,td,paird

    def label_ring_atoms(self):
        cycles=self.TopoCoord.ring_detector()
        adf=self.TopoCoord.gro_DataFrame('atoms')
        self.TopoCoord.set_gro_attribute('cycle-idx',np.zeros(adf.shape[0]).astype(int))
        cidx=1
        for l,cl in cycles.items():
            for c in cl:
                for idx in c:
                    self.TopoCoord.set_gro_attribute_by_attributes('cycle-idx',cidx,{'globalIdx':idx})
                cidx+=1
        # logging.debug(f'label_ring_atoms for {self.name}:\n{adf.to_string()}')

    def get_resname(self,internal_resid):
        # logging.debug(f'{self.name} sequence: {self.sequence}')
        return self.sequence[internal_resid-1]

    def inherit_attribute_from_reactants(self,attribute,available_molecules,increment=True):
        adf=self.TopoCoord.Coordinates.A
        ordered_attribute_idx=[]
        curr_max=0
        logging.debug(f'{self.name}({adf.shape[0]}) inheriting {attribute} from {self.sequence}')
        for i,r in enumerate(self.sequence):
            namesinres=list(adf[adf['resNum']==(i+1)]['atomName'])
            rdf=available_molecules[r].TopoCoord.Coordinates.A
            x=list(rdf[rdf['atomName'].isin(namesinres)]['sea-idx'])
            logging.debug(f'{r}->{len(x)}')
            if increment:
                x=[y+curr_max for y in x]
                curr_max=max(x)
            ordered_attribute_idx.extend(x)
        assert len(ordered_attribute_idx)==adf.shape[0]
        adf[attribute]=ordered_attribute_idx

    def merge(self,other):
        self.TopoCoord.merge(other.TopoCoord)

    def load_top_gro(self,topfilename,grofilename,mol2filename=''):
        self.TopoCoord=TopoCoord(topfilename=topfilename,grofilename=grofilename,mol2filename=mol2filename)

    def analyze_sea_topology(self):
        self.TopoCoord.analyze_sea_topology()

    def set_gro_attribute(self,attribute,srs):
        self.TopoCoord.set_gro_attribute(attribute,srs)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
        self.TopoCoord.read_gro_attributes(grxfilename,attribute_list=attribute_list)

    def write_gro_attributes(self,attribute_list,grxfilename):
        self.TopoCoord.write_gro_attributes(attribute_list,grxfilename)

    # def update_topology(self,t):
    #     self.Topology.merge(t)

    # def update_coords(self,c):
    #     self.Coords.copy_coords(c)

    def make_bonds(self):
        bonds=[]
        hs_from_tr=[]
        skip_H=[]
        for i,B in enumerate(self.reaction_bonds):
            (aidx,bidx),(aresid,bresid),(aname,bname)=B
            logging.debug(f'generating {self.name} bond {i} {aresid}:{aname}-{bresid}:{bname}')
            bonds.append((aidx,bidx))
            if aresid!=bresid:
                # transrot identifies the two sacrificial H's
                hxi,hxj=self.transrot(aidx,aresid,bidx,bresid)
                skip_H.append(i)
                hs_from_tr.append(hxi)
                hs_from_tr.append(hxj)
        # make_bonds returns list of sacrificial H idx's derived
        # from bonds whose indices are not in skip_H
        idx_scratch=self.TopoCoord.make_bonds(bonds,skip_H=skip_H)
        # the H's identified in transrot plus those identified by make_bonds
        # must be deleted
        idx_scratch.extend(hs_from_tr)
        return self.TopoCoord.delete_atoms(idx_scratch)

    # def new_bond(self,other,at_idx=-1,from_idx=-1):
    #     if self!=other:
    #         hxi,hxj=self.transrot(other,at_idx,from_idx)
    #         idx_shift=self.TopoCoord.num_atoms()
    #         from_idx+=idx_shift
    #         hxj+=idx_shift
    #         idx_scratch0=[hxi,hxj]
    #         self.merge(other)
    #     bond=(at_idx,from_idx)
    #     idx_scratch=self.TopoCoord.make_bonds([bond])
    #     if self!=other:
    #         logging.debug(f'Making {self.name}: transrot wants to delete {idx_scratch0} and make_bonds wants {idx_scratch}')
    #         idx_mapper=self.TopoCoord.delete_atoms(idx_scratch0)
    #     else:
    #         idx_mapper=self.TopoCoord.delete_atoms(idx_scratch)
    #     return idx_mapper[at_idx],idx_mapper[from_idx]

    # def add_bonds(self,pairs=[]):
    #     self.Topology.add_bonds(pairs,enumerate_others=False)

    # def delete_atoms(self,idx=[],return_idx_of=[]):
    #     new_idx=self.Topology.delete_atoms(idx,return_idx_of=return_idx_of)
    #     self.Coords.delete_atoms(idx)
    #     return new_idx

    def transrot(self,at_idx,at_resid,from_idx,from_resid):
        # Rotate and translate 
        # logging.debug('transrot for building {self.name} from {at_resid}:{at_idx} -- {from_resid}:{from_idx}')
        if at_resid==from_resid:
            return
        TC=self.TopoCoord
        ATC=TopoCoord()
        BTC=TopoCoord()
        C=TC.gro_DataFrame('atoms')
        ATC.Coordinates.A=C[C['resNum']==at_resid].copy()
        BTC.Coordinates.A=C[C['resNum']==from_resid].copy()
        mypartners=TC.partners_of(at_idx)
        otpartners=TC.partners_of(from_idx)
        # logging.debug(f'Partners of {at_idx} {mypartners}')
        # logging.debug(f'Partners of {from_idx} {otpartners}')
        myHpartners={k:v for k,v in zip(mypartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners={k:v for k,v in zip(otpartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        myHighestH={k:v for k,v in myHpartners.items() if v==max([k for k in myHpartners.values()])}
        otHighestH={k:v for k,v in otHpartners.items() if v==max([k for k in otHpartners.values()])}
        assert len(myHighestH)==1
        assert len(otHighestH)==1
        # logging.debug(f'Highest-named H partner of {at_idx} is {myHighestH}')
        # logging.debug(f'Highest-named H partner of {from_idx} is {otHighestH}')
        assert len(myHpartners)>0,f'Error: atom {at_idx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {from_idx} does not have a deletable H atom!'

        Ri=TC.get_R(at_idx)
        Rj=TC.get_R(from_idx)
        overall_maximum=(-1.e9,-1,-1)
        coord_trials={}
        for myH,myHnm in myHpartners.items():  # keys are globalIdx's, values are names
            coord_trials[myH]={}
            Rh=TC.get_R(myH)
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            for otH,otHnm in otHpartners.items():
                # logging.debug(f'{self.name}: Considering {myH} {otH}')
                coord_trials[myH][otH]=deepcopy(BTC)
                # logging.debug(f'\n{coord_trials[myH][otH].Coordinates.A.to_string()}')
                Rk=coord_trials[myH][otH].get_R(otH)
                # logging.debug(f'{self.name}:    otH {otH} Rk {Rk}')
                Rkj=Rk-Rj
                Rkj*=1.0/np.linalg.norm(Rkj)
                #Rhk=Rh-Rk
                #rhk=np.linalg.norm(Rhk)
                cp=np.cross(Rkj,Rih)
                c=np.dot(Rkj,Rih)
                v=np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
                v2=np.dot(v,v)
                I=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                # R is the rotation matrix that will rotate donb to align with accb
                R=I+v+v2/(1.+c)
                # logging.debug(f'{self.name}: R:\n{R}')
                # rotate translate all donor atoms!
                coord_trials[myH][otH].rotate(R)
                Rk=coord_trials[myH][otH].get_R(otH)
                # overlap the two H atoms by translation
                Rik=Rh-Rk
                coord_trials[myH][otH].translate(Rik)
                # coord_trials[myH][otH].Coordinates.write_mol2(f'{self.name}-{myH}-{otH}.mol2')
                minD=TC.minimum_distance(coord_trials[myH][otH],self_excludes=[myH],other_excludes=[otH])
                # logging.debug(f'{self.name}: minD {minD}')
                if minD>overall_maximum[0]:
                    overall_maximum=(minD,myH,otH)
        # logging.debug(f'{self.name}: overall_maximum {overall_maximum}')
        minD,myH,otH=overall_maximum
        BTC=coord_trials[myH][otH]
        TC.overwrite_coords(BTC)
        TC.swap_atom_names(myH,list(myHighestH.keys())[0])
        TC.swap_atom_names(otH,list(otHighestH.keys())[0])
        return myH,otH

    def atoms_w_same_attribute_as(self,find_dict={},same_attribute='',return_attribute=''):
        att_val=self.TopoCoord.get_gro_attribute_by_attributes(same_attribute,find_dict)
        return self.TopoCoord.get_gro_attributelist_by_attributes(return_attribute,{same_attribute:att_val})