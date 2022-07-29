from itertools import chain, product
import os
from copy import deepcopy
import pandas as pd
import numpy as np
import logging

from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.bondtemplate import BondTemplate,BondTemplateList,ReactionBond,ReactionBondList
from HTPolyNet.coordinates import _dfrotate
from HTPolyNet.ambertools import GAFFParameterize
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.gromacs import grompp_and_mdrun, mdp_library

logger=logging.getLogger(__name__)

def _rotmat(axis,radians):
    R=np.identity(3)
    sr=np.sin(radians)
    cr=np.cos(radians)
    if axis==2:
        R[0][0]=cr
        R[0][1]=-sr
        R[1][0]=sr
        R[1][1]=cr
    elif axis==1:
        R[0][0]=cr
        R[0][2]=sr
        R[2][0]=-sr
        R[2][2]=cr
    elif axis==0:
        R[1][1]=cr
        R[1][2]=-sr
        R[2][1]=sr
        R[2][2]=cr
    return R

class Reaction:
    def __init__(self,jsondict={}):
        self.jsondict=jsondict
        self.name=jsondict.get('name','')
        self.atoms=jsondict.get('atoms',{})
        self.bonds=jsondict.get('bonds',{})
        self.reactants=jsondict.get('reactants',{})
        self.product=jsondict.get('product','')
        self.restrictions=jsondict.get('restrictions',{})
        self.stage=jsondict.get('stage','')
        self.probability=jsondict.get('probability',1.0)
        self.symmetry_versions=[]
        
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

ReactionList = list[Reaction]

def is_reactant(name:str,reaction_list:ReactionList,stage='cure'):
    reactants=[]
    for r in reaction_list:
        if r.stage==stage:
            for v in r.reactants.values():
                if not v in reactants:
                    reactants.append(v)
    return name in reactants

def yield_bonds(R:Reaction,TC:TopoCoord,resid_mapper):
    nreactants=len(R.reactants)
    for bondrec in R.bonds:
        atom_keys=bondrec['atoms']
        order=bondrec['order']
        assert len(atom_keys)==2
        atomrecs=[R.atoms[x] for x in atom_keys]
        atom_names=[x['atom'] for x in atomrecs]
        in_reactant_resids=[x['resid'] for x in atomrecs]
        if nreactants==1:
            in_product_resids=[resid_mapper[0][in_reactant_resids[x]] for x in [0,1]]
        else:
            in_product_resids=[resid_mapper[x][in_reactant_resids[x]] for x in [0,1]]
        atom_idx=[TC.get_gro_attribute_by_attributes('globalIdx',{'resNum':in_product_resids[x],'atomName':atom_names[x]}) for x in [0,1]]
        bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames=TC.get_bystanders(atom_idx)
        oneaway_resids,oneaway_resnames,oneaway_atomidx,onewaway_atomnames=TC.get_oneaways(atom_idx)
        yield ReactionBond(atom_idx,in_product_resids,order,bystander_resids,bystander_atomidx,oneaway_resids,oneaway_atomidx)

class Molecule:
    def __init__(self,name='',generator:Reaction=None):
        self.name=name
        self.TopoCoord=TopoCoord()
        self.generator:Reaction=generator
        self.sequence=[]
        self.origin=None
        self.reaction_bonds:ReactionBondList=[]
        self.bond_templates:BondTemplateList=[]
        self.symmetry_relateds=[]
        self.stereocenters=[]
        self.stereoisomers=[]
        self.zrecs=[]
        self.is_reactant=False

    def set_origin(self,value):
        self.origin=value

    def get_origin(self):
        return self.origin

    def update_zrecs(self,zrecs):
        for zr in zrecs:
            if not zr in self.zrecs:
                self.zrecs.append(zr)

    def initialize_molecule_cycles(self):
        TC=self.TopoCoord
        TC.idx_lists['cycle']=[]
        cycle_dict=TC.Topology.detect_cycles()
        for l,cs_of_l in cycle_dict.items():
            TC.idx_lists['cycle'].extend(cs_of_l)
        TC.reset_grx_attributes_from_idx_list('cycle')

    def initialize_monomer_grx_attributes(self):
        TC=self.TopoCoord
        TC.set_gro_attribute('z',0)
        TC.set_gro_attribute('nreactions',0)
        for att in ['sea_idx','chain','chain_idx','cycle','cycle_idx']:
            TC.set_gro_attribute(att,-1)
        # set symmetry class indices
        sea_idx=1
        logger.debug(f'{self.name}: symmetry_relateds {self.symmetry_relateds}')
        for s in self.symmetry_relateds:
            logger.debug(f'sea_idx {sea_idx} set for set {s}')
            for atomName in s:
                logger.debug(f'{atomName} {sea_idx}')
                TC.set_gro_attribute_by_attributes('sea_idx',sea_idx,{'atomName':atomName})
            sea_idx+=1
        # set z and nreactions
        idx=[]
        for zr in self.zrecs:
            an=zr['atom']
            rnum=zr['resid']
            z=zr['z']
            TC.set_gro_attribute_by_attributes('z',z,{'atomName':an,'resNum':rnum})
            idx.append(TC.get_gro_attribute_by_attributes('globalIdx',{'atomName':an,'resNum':rnum}))
            for sr in self.symmetry_relateds:
                a,b=sr
                if a==an:
                    idx.append(TC.get_gro_attribute_by_attributes('globalIdx',{'atomName':b,'resNum':rnum}))
                    TC.set_gro_attribute_by_attributes('z',z,{'atomName':b,'resNum':rnum})
                elif b==an:
                    idx.append(TC.get_gro_attribute_by_attributes('globalIdx',{'atomName':a,'resNum':rnum}))
                    TC.set_gro_attribute_by_attributes('z',z,{'atomName':a,'resNum':rnum})

        # set chain, chain_idx
        TC.idx_lists['chain']=[]
        pairs=product(idx,idx)
        for i,j in pairs:
            if i<j:
                if TC.are_bonded(i,j):
                    # this monomer has two atoms capable of reacting
                    # that are bound to each other -- this means that
                    # the two originated as a double-bond.
                    # *If* there is one with three hydrogens
                    # (remember this is an activated monomer)
                    # then it is the "tail"; the other is the "head".
                    i_nH=TC.count_H(i)
                    j_nH=TC.count_H(j)
                    if i_nH==3 and j_nH!=3:
                        # i is the tail
                        entry=[j,i]
                    elif i_nH!=3 and j_nH==3:
                        # j is the tail
                        entry=[i,j]
                    else:
                        logger.warning(f'In molecule {self.name}, cannot identify head and tail atoms in reactive double bond\nAssuming {j} is head and {i} is tail')
                        entry=[j,i]
                    # logger.debug(f'Adding {entry} to chainlist of {self.name}')
                    TC.idx_lists['chain'].append(entry)
        TC.reset_grx_attributes_from_idx_list('chain')
        # set cycle, cycle_idx
        self.initialize_molecule_cycles()

    def previously_parameterized(self):
        rval=True
        for ext in ['mol2','top','itp','gro']:
            rval=rval and pfs.exists(os.path.join('molecules/parameterized',f'{self.name}.{ext}'))
        return rval

    def parameterize(self,outname='',input_structure_format='mol2',**kwargs):
        assert os.path.exists(f'{self.name}.{input_structure_format}'),f'Cannot parameterize molecule {self.name} without {self.name}.{input_structure_format} as input'
        if outname=='':
            outname=f'{self.name}'
        GAFFParameterize(self.name,outname,input_structure_format=input_structure_format,**kwargs)
        self.load_top_gro(f'{outname}.top',f'{outname}.gro',mol2filename=f'{outname}.mol2',wrap_coords=False)

    def minimize(self,outname='',**kwargs):
        if outname=='':
            outname=f'{self.name}'
        n=self.name
        boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix=mdp_library['minimize-single-molecule']
        pfs.checkout(f'mdp/{mdp_prefix}.mdp')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp=mdp_prefix,out=f'{outname}',boxSize=boxsize)
        self.TopoCoord.read_gro(f'{n}.gro',wrap_coords=False)

    def center_coords(self,new_boxsize:np.ndarray=None):
        if type(new_boxsize)==np.ndarray:
            if new_boxsize.shape==(3,):
                box_vectors=new_boxsize*np.identity(3,dtype=float)
                logger.debug(f'{box_vectors}')
            elif new_boxsize.shape==(3,3):
                box_vectors=new_boxsize
            self.TopoCoord.Coordinates.box=box_vectors
        center=self.TopoCoord.Coordinates.box.diagonal()/2.0
        gc=self.TopoCoord.Coordinates.geometric_center()
        addme=center-gc
        self.TopoCoord.Coordinates.translate(addme)

    def generate(self,outname='',available_molecules={},**kwargs):
        # logger.info(f'Generating {self.name}.mol2 for parameterization')
        if outname=='':
            outname=f'{self.name}'
        if self.generator:
            R=self.generator
            self.TopoCoord=TopoCoord()
            logger.debug(f'Using reaction {R.name} to generate {self.name}.mol2.')
            isf='mol2'
            resid_mapper=[]
            for ri in R.reactants.values():
                new_reactant=deepcopy(available_molecules[ri])
                new_reactant.TopoCoord.write_mol2(filename=f'{self.name}-reactant{ri}-prebonding.mol2',molname=self.name)
                rnr=len(new_reactant.sequence)
                shifts=self.TopoCoord.merge(new_reactant.TopoCoord)
                resid_mapper.append({k:v for k,v in zip(range(1,rnr+1),range(1+shifts[2],1+rnr+shifts[2]))})
            # logger.debug(f'{self.name}: resid_mapper {resid_mapper}')
            # logger.debug(f'{self.TopoCoord.idx_lists}')
            # logger.debug(f'\n{self.TopoCoord.Coordinates.A.to_string()}')
            logger.debug(f'composite prebonded molecule in box {self.TopoCoord.Coordinates.box}')
            self.TopoCoord.write_mol2(filename=f'{self.name}-prebonding.mol2',molname=self.name)
            self.set_sequence()
            bonds_to_make=list(yield_bonds(R,self.TopoCoord,resid_mapper))
            # logger.debug(f'Generation of {self.name}: composite molecule has {len(self.sequence)} resids')
            # logger.debug(f'generation of {self.name}: composite molecule:\n{composite_mol.TopoCoord.Coordinates.A.to_string()}')
            idx_mapper=self.make_bonds(bonds_to_make)
            self.TopoCoord.set_gro_attribute('reactantName',R.product)
            self.TopoCoord.set_gro_attribute('sea_idx',-1) # turn off symmetry-equivalence for multimers
            self.write_gro_attributes(['z','nreactions','reactantName','sea_idx','cycle','cycle_idx','chain','chain_idx'],f'{R.product}.grx')
            self.TopoCoord.write_mol2(filename=f'{self.name}.mol2',molname=self.name)
        else:
            # this is a monomer; we need an input structure file to feed antechamber
            input_structure_formats=['mol2','pdb']
            isf=None
            for isf in input_structure_formats:
                if pfs.exists(f'molecules/inputs/{self.name}.{isf}'):
                    logger.debug(f'Using input molecules/inputs/{self.name}.{isf} as a generator.')
                    pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
                    break
            assert isf,'Error: no valid input structure file found'

        reactantName=self.name
        self.parameterize(outname,input_structure_format=isf,**kwargs)
        self.minimize(outname,**kwargs)
        self.set_sequence()
        self.TopoCoord.set_gro_attribute('reactantName',reactantName)
        if not self.generator:
            self.initialize_monomer_grx_attributes()
            self.write_gro_attributes(['z','nreactions','reactantName','sea_idx','cycle','cycle_idx','chain','chain_idx'],f'{reactantName}.grx')
        else:
            grx=f'{reactantName}.grx'
            if (os.path.exists(grx)):
                self.TopoCoord.read_gro_attributes(grx)
                #self.reset_chains_from_attributes()
        # logger.debug(f'{self.name} gro\n{self.TopoCoord.Coordinates.A.to_string()}')
        self.prepare_new_bonds(available_molecules=available_molecules)


    # def set_reaction_bonds(self,available_molecules={}):
    def prepare_new_bonds(self,available_molecules={}):
        # logger.debug(f'set_reaction_bonds: molecules {list(available_molecules.keys())}')
        R=self.generator
        if not R:
            return
        self.reaction_bonds=[]
        self.bond_templates=[]
        TC=self.TopoCoord
        # logger.debug(f'prepare_new_bonds {self.name}: chainlists {TC.idx_lists["chain"]}')
        for bondrec in R.bonds:
            atom_keys=bondrec['atoms']
            order=bondrec['order']
            assert len(atom_keys)==2
            atomrecs=[R.atoms[x] for x in atom_keys]
            atom_names=[x['atom'] for x in atomrecs]
            reactant_keys=[x['reactant'] for x in atomrecs]
            in_reactant_resids=[x['resid'] for x in atomrecs]
            if reactant_keys[0]==reactant_keys[1]:  # this is an intraresidue bond
                reactant_names=[R.reactants[reactant_keys[0]]]
            else:
                reactant_names=[R.reactants[x] for x in reactant_keys]
            reactant_sequences=[available_molecules[x].sequence for x in reactant_names]
            product_sequence=[]
            for seq in reactant_sequences:
                product_sequence.extend(seq)
            # logger.debug(f'product_sequence {product_sequence}')
            sequence_residue_idx_origins=[0,0]
            if len(reactant_sequences)==2:
                sequence_residue_idx_origins[1]=len(reactant_sequences[0])
            in_product_resids=[in_reactant_resids[x]+sequence_residue_idx_origins[x] for x in [0,1]]
            # logger.debug(f'in_product_resids {in_product_resids}')
            in_product_resnames=[product_sequence[in_product_resids[x]-1] for x in [0,1]]
            atom_idx=[TC.get_gro_attribute_by_attributes('globalIdx',{'resNum':in_product_resids[x],'atomName':atom_names[x]}) for x in [0,1]]
            logger.debug(f'{R.name} names {atom_names} in_product_resids {in_product_resids} idx {atom_idx}')
            bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames=TC.get_bystanders(atom_idx)
            oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames=TC.get_oneaways(atom_idx)
            # logger.debug(f'{self.name} bystanders {bystander_resids} {bystander_resnames} {bystander_atomidx} {bystander_atomnames}')
            # logger.debug(f'{self.name} oneaways {oneaway_resids} {oneaway_resnames} {oneaway_atomidx} {oneaway_atomnames}')
            self.reaction_bonds.append(ReactionBond(atom_idx,in_product_resids,order,bystander_resids,bystander_atomidx,oneaway_resids,oneaway_atomidx))
            intraresidue=in_product_resids[0]==in_product_resids[1]
            self.bond_templates.append(BondTemplate(atom_names,in_product_resnames,intraresidue,order,bystander_resnames,bystander_atomnames,oneaway_resnames,oneaway_atomnames))

    def set_sequence(self):
        """set_sequence Establish the sequence-list (residue names in order) based on resNum attributes in atom list
        """
        adf=self.TopoCoord.gro_DataFrame('atoms')
        self.sequence=[]
        current_resid=0
        for i,r in adf.iterrows():
            ri=r['resNum']
            rn=r['resName']
            if ri!=current_resid:
                current_resid=ri
                self.sequence.append(rn)
        # logger.debug(f'{self.name} sequence: {self.sequence}')

    def idx_mappers(self,otherTC:TopoCoord,other_bond,bystanders,oneaways,uniq_atom_idx:set):
        assert len(other_bond)==2
        assert len(bystanders)==2
        assert len(oneaways)==2
        ut=uniq_atom_idx.copy()
        logger.debug(f'Template name {self.name}')
        i_idx,j_idx=other_bond
        i_resName,i_resNum,i_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i_idx})
        j_resName,j_resNum,j_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j_idx})
        logger.debug(f'i_idx {i_idx} i_resName {i_resName} i_resNum {i_resNum} i_atomName {i_atomName}')
        logger.debug(f'j_idx {j_idx} j_resName {j_resName} j_resNum {j_resNum} j_atomName {j_atomName}')
        # identify the template bond represented by the other_bond parameter
        ij=[]
        for RB,BT in zip(self.reaction_bonds,self.bond_templates):
            temp_resids=RB.resids
            temp_iname,temp_jname=BT.names
            temp_iresname,temp_jresname=BT.resnames
            temp_bystander_resids=RB.bystander_resids
            temp_oneaway_resids=RB.oneaway_resids
            logger.debug(f'temp_iresname {temp_iresname} temp_iname {temp_iname}')
            logger.debug(f'temp_jresname {temp_jresname} temp_jname {temp_jname}')
            if (i_atomName,i_resName)==(temp_iname,temp_iresname):
                ij=[0,1]
                break # found it -- stop looking
            elif (i_atomName,i_resName)==(temp_jname,temp_jresname):
                ij=[1,0]
                break
        assert len(ij)==2,f'Mappers using template {self.name} unable to map from instance bond {i_resName}-{i_resNum}-{i_atomName}---{j_resName}-{j_resNum}-{j_atomName}'
        inst_resids=[i_resNum,j_resNum]
        inst_resids=[inst_resids[ij[x]] for x in [0,1]]
        inst_bystander_resids=[bystanders[ij[x]] for x in [0,1]]
        inst_oneaway_resids=[oneaways[ij[x]] for x in [0,1]]
        assert all([len(inst_bystander_resids[x])==len(temp_bystander_resids[x]) for x in [0,1]]),f'Error: bystander count mismatch'
           # use dataframe merges to create globalIdx maps
        instdf=otherTC.Coordinates.A
        tempdf=self.TopoCoord.Coordinates.A
        inst2temp={}
        temp2inst={}
        # logger.debug(f'inst resids from {[i_resNum,j_resNum,*inst_bystanders]}') 
        # logger.debug(f'temp resids from {[temp_iresid,temp_jresid,*temp_bystanders]}')
        for inst,temp in zip([*inst_resids,*inst_bystander_resids[0],*inst_bystander_resids[1],*inst_oneaway_resids],
                             [*temp_resids,*temp_bystander_resids[0],*temp_bystander_resids[1],*temp_oneaway_resids]):
            if inst and temp:  # None's in the bystander lists and oneaways lists should be ignored
                logger.debug(f'map inst resid {inst} to template resid {temp}')
                idf=instdf[instdf['resNum']==inst][['globalIdx','atomName']].copy()
                # logger.debug(f'idf res {inst}:\n{idf.to_string()}')
                tdf=tempdf[tempdf['resNum']==temp][['globalIdx','atomName']].copy()
                # logger.debug(f'tdf res {temp}:\n{tdf.to_string()}')
                tdf=tdf.merge(idf,on='atomName',how='inner',suffixes=('_template','_instance'))
                # logger.debug(f'merged\n{tdf.to_string()}')
                for i,r in tdf.iterrows():
                    temp_idx=r['globalIdx_template']
                    inst_idx=r['globalIdx_instance']
                    # logger.debug(f't {temp_idx} <-> i {inst_idx}')
                    # only map template atoms that are identified in the passed in set
                    if temp_idx in ut:
                        ut.remove(temp_idx)
                        temp2inst[temp_idx]=inst_idx
                        inst2temp[inst_idx]=temp_idx
                    if temp_idx in temp2inst and temp2inst[temp_idx]!=inst_idx:
                        raise Exception(f'Error: temp_idx {temp_idx} already claimed in temp2inst; bug')
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
            logger.error('NAN in molecule/dihedrals')
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
            logger.error('NAN in molecule/pairs')
            raise Exception
        return ad,td,paird

    def label_cycle_atoms(self):
        self.TopoCoord.label_cycle_atoms()

    def get_resname(self,internal_resid):
        # logger.debug(f'{self.name} sequence: {self.sequence}')
        return self.sequence[internal_resid-1]

    def inherit_attribute_from_reactants(self,attribute,available_molecules,increment=True,no_increment_if_negative=True):
        adf=self.TopoCoord.Coordinates.A
        ordered_attribute_idx=[]
        curr_max=0
        # logger.debug(f'{self.name}({adf.shape[0]}) inheriting {attribute} from {self.sequence}')
        # logger.debug(f'available molecules {list(available_molecules.keys())}')
        for i,r in enumerate(self.sequence):
            '''
            for this residue number, read the list of unique atom names
            '''
            namesinres=list(adf[adf['resNum']==(i+1)]['atomName'])
            '''
            access coordinates of standalone residue template with this name 'r' on the list of available molecules
            '''
            rdf=available_molecules[r].TopoCoord.Coordinates.A
            '''
            get the attribute values from residue template
            '''
            x=list(rdf[rdf['atomName'].isin(namesinres)][attribute])
            # logger.debug(f'{r}->{len(x)}')
            '''
            increment these attribute value based on residue number in this molecule
            '''
            if increment:
                i_x=[]
                for y in x:
                    if y>0 or (y<0 and not no_increment_if_negative):
                        i_x.append(y+curr_max)
                    else:
                        i_x.append(y)
                curr_max=max(i_x)
            ordered_attribute_idx.extend(i_x)
        assert len(ordered_attribute_idx)==adf.shape[0]
        adf[attribute]=ordered_attribute_idx

    def merge(self,other):
        shifts=self.TopoCoord.merge(other.TopoCoord)
        return shifts

    def load_top_gro(self,topfilename,grofilename,mol2filename='',**kwargs):
        wrap_coords=kwargs.get('wrap_coords',True)
        self.TopoCoord=TopoCoord(topfilename=topfilename,grofilename=grofilename,mol2filename=mol2filename,wrap_coords=wrap_coords)
        # logger.debug(f'box: {self.TopoCoord.Coordinates.box}')

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

    def make_bonds(self,bondrecs:ReactionBondList):
        bonds=[]
        hs_from_tr=[]
        skip_H=[]
        TC=self.TopoCoord
        # logger.debug(f'{self.name} make_bonds: chainlists: {TC.idx_lists["chain"]}')
        # for cl in TC.idx_lists['chain']:
        #     nl=[self.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl]
        #     logger.debug(f'{nl}')
        for i,B in enumerate(bondrecs):
            # logger.debug(f'{self.name} reaction bond {i}: {str(B)} ')
            aidx,bidx=B.idx
            aname,bname=[TC.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in [aidx,bidx]]
            aresid,bresid=B.resids
            bystanders=B.bystander_resids
            oneaways=B.oneaway_resids
            order=B.order
            logger.debug(f'generating {self.name} bond {aresid}:{aname}:{aidx}-{bresid}:{bname}:{bidx} order {order}')
            # logger.debug(f'bystander resids {bystanders}')
            # logger.debug(f'oneaway resids {oneaways}')
            bonds.append((aidx,bidx,order))
            if aresid!=bresid:
                # transrot identifies the two sacrificial H's
                cresids=[]
                if len(oneaways)==2 and oneaways[1]!=None:
                    cresids=[oneaways[1]]
                if len(bystanders)==2 and any(bystanders[1]):
                    cresids.extend(bystanders[1])
                # cresids=[x for x in oneaways if x]
                # cresids.extend([x for x in bystanders if x])
                # logger.debug(f'cresids {cresids}')
                hxi,hxj=self.transrot(aidx,aresid,bidx,bresid,connected_resids=cresids)
                skip_H.append(i)
                hs_from_tr.append(hxi)
                hs_from_tr.append(hxj)
        # make_bonds returns list of sacrificial H idx's derived
        # from bonds whose indices are not in skip_H
        idx_scratch=TC.make_bonds(bonds,skip_H=skip_H)
        # the H's identified in transrot plus those identified by make_bonds
        # must be deleted
        idx_scratch.extend(hs_from_tr)
        idx_mapper=TC.delete_atoms(idx_scratch)
        #TC.remap_idx_list('chain',idx_mapper)
        for B in bondrecs:
#        for i,B in enumerate(self.reaction_bonds):
            aidx,bidx=B.idx
            # (aidx,bidx),(aresid,bresid),(aoresids,boresids),(aname,bname),order=B
            aidx=idx_mapper[aidx]
            bidx=idx_mapper[bidx]
            TC.decrement_gro_attribute_by_attributes('z',{'globalIdx':aidx})
            TC.decrement_gro_attribute_by_attributes('z',{'globalIdx':bidx})
            TC.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':aidx})
            TC.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':bidx})
        self.initialize_molecule_cycles()
        # cb=self.TopoCoord.checkbox()
        # logger.debug(f'checkbox: {cb}')
        return idx_mapper

    def transrot(self,at_idx,at_resid,from_idx,from_resid,connected_resids=[]):
        # Rotate and translate
        if at_resid==from_resid:
            return
        TC=self.TopoCoord
        ATC=TopoCoord()
        BTC=TopoCoord()
        C=TC.gro_DataFrame('atoms')
        ATC.Coordinates.A=C[C['resNum']==at_resid].copy()
        bresids=connected_resids.copy()
        bresids.append(from_resid)
        BTC.Coordinates.A=C[C['resNum'].isin(bresids)].copy()
        NONROT=C[~C['resNum'].isin(bresids)].shape[0]
        logger.debug(f'{self.TopoCoord.Coordinates.A.shape[0]} atoms')
        logger.debug(f'holding {at_resid} ({NONROT})')
        logger.debug(f'rotating/translating {bresids} ({BTC.Coordinates.A.shape[0]})')
        assert self.TopoCoord.Coordinates.A.shape[0]==(NONROT+BTC.Coordinates.A.shape[0])
        mypartners=TC.partners_of(at_idx)
        otpartners=TC.partners_of(from_idx)
        # logger.debug(f'Partners of {at_idx} {mypartners}')
        # logger.debug(f'Partners of {from_idx} {otpartners}')
        myHpartners={k:v for k,v in zip(mypartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners={k:v for k,v in zip(otpartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        myHighestH={k:v for k,v in myHpartners.items() if v==max([k for k in myHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
        otHighestH={k:v for k,v in otHpartners.items() if v==max([k for k in otHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
        assert len(myHighestH)==1
        assert len(otHighestH)==1
        # logger.debug(f'Highest-named H partner of {at_idx} is {myHighestH}')
        # logger.debug(f'Highest-named H partner of {from_idx} is {otHighestH}')
        assert len(myHpartners)>0,f'Error: atom {at_idx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {from_idx} does not have a deletable H atom!'

        Ri=TC.get_R(at_idx)
        Rj=TC.get_R(from_idx)
        overall_maximum=(-1.e9,-1,-1)
        coord_trials={}
        for myH,myHnm in myHpartners.items():  # keys are globalIdx's, values are names
            coord_trials[myH]:dict[TopoCoord]={}
            Rh=TC.get_R(myH)
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            for otH,otHnm in otHpartners.items():
                # logger.debug(f'{self.name}: Considering {myH} {otH}')
                coord_trials[myH][otH]=deepcopy(BTC)
                # logger.debug(f'\n{coord_trials[myH][otH].Coordinates.A.to_string()}')
                Rk=coord_trials[myH][otH].get_R(otH)
                # logger.debug(f'{self.name}:    otH {otH} Rk {Rk}')
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
                # logger.debug(f'{self.name}: R:\n{R}')
                # rotate translate all donor atoms!
                coord_trials[myH][otH].rotate(R)
                Rk=coord_trials[myH][otH].get_R(otH)
                # overlap the two H atoms by translation
                Rik=Rh-Rk
                coord_trials[myH][otH].translate(Rik)
                minD=TC.minimum_distance(coord_trials[myH][otH],self_excludes=[myH],other_excludes=[otH])
                # logger.debug(f'{self.name}: minD {minD}')
                if minD>overall_maximum[0]:
                    overall_maximum=(minD,myH,otH)
        # logger.debug(f'{self.name}: overall_maximum {overall_maximum}')
        minD,myH,otH=overall_maximum
        BTC=coord_trials[myH][otH]
        TC.overwrite_coords(BTC)
        TC.swap_atom_names(myH,list(myHighestH.keys())[0])
        TC.swap_atom_names(otH,list(otHighestH.keys())[0])
        return myH,otH

    def atoms_w_same_attribute_as(self,find_dict={},same_attribute='',return_attribute=''):
        att_val=self.TopoCoord.get_gro_attribute_by_attributes(same_attribute,find_dict)
        return self.TopoCoord.get_gro_attributelist_by_attributes(return_attribute,{same_attribute:att_val})

    def flip_stereocenter(self,idx):
        """flip_stereocenter Flips stereochemistry of atom at idx

        :param idx: global index of atom
        :type idx: int
        """
        TC=self.TopoCoord
        A=TC.Coordinates.A
        logger.debug(f'{self.name} flipping on {idx}')
        ligand_idx=self.TopoCoord.Topology.bondlist.partners_of(idx)
        if len(ligand_idx)!=4:
            logger.debug(f'Atom {idx} cannot be a stereocenter; it only has {len(ligand_idx)} ligands.')
            return
        # determine the two lightest ligands
        branches={}
        for n in ligand_idx:
            bl=deepcopy(self.TopoCoord.Topology.bondlist)
            branches[n]=bl.half_as_list([idx,n],3)
        ligand_idx.sort(key=lambda x: len(branches[x]))
        a=ligand_idx[0]
        aset=list(set(list(chain.from_iterable(branches[a]))))
        b=ligand_idx[1]
        bset=list(set(list(chain.from_iterable(branches[b]))))
        # translate origin to location of stereocenter
        O=TC.get_R(idx)
        TC.translate(-1*O)
        rO=TC.get_R(idx)
        ra=TC.get_R(a)
        rb=TC.get_R(b)
        adf=A[A['globalIdx'].isin(aset)]
        bdf=A[A['globalIdx'].isin(bset)]
        # rotate the two branches to swap them
        rOa=rO-ra
        rOa*=1.0/np.linalg.norm(rOa)
        rOb=rO-rb
        rOb*=1.0/np.linalg.norm(rOb)
        cp=np.cross(rOa,rOb)
        c=np.dot(rOa,rOb)
        v=np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
        v2=np.dot(v,v)
        I=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        R=I+v+v2/(1.+c)
        _dfrotate(adf,R)
        A.loc[A['globalIdx'].isin(aset),['posX','posY','posZ']]=adf[['posX','posY','posZ']]
        cp=np.cross(rOb,rOa)
        v=np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
        v2=np.dot(v,v)
        R=I+v+v2/(1.+c)
        _dfrotate(bdf,R)
        A.loc[A['globalIdx'].isin(bset),['posX','posY','posZ']]=bdf[['posX','posY','posZ']]
        # translate back to original coordinate frame
        TC.translate(O)

    def rotate_bond(self,a,b,deg):
        """Rotates all atoms in molecule on b-side of a-b bond by deg degrees

        :param a: index of a
        :type a: int
        :param b: index of b
        :type b: int
        :param deg: angle of rotation (degrees)
        :type deg: float
        """
        TC=self.TopoCoord
        A=TC.Coordinates.A
        branch=TopoCoord()
        bl=deepcopy(self.TopoCoord.Topology.bondlist)
        branchidx=bl.half_as_list((a,b),99)
        ra=TC.get_R(a)
        rb=TC.get_R(b)
        Rab=ra-rb
        rab=Rab/np.linalg.norm(Rab)
        O=rb
        TC.translate(-1*O)
        branch.Coordinates.A=A[A['globalIdx'].isin(branchidx)].copy()
        # do stuff
        rx,ry,rz=rab
        caz=rx/(rx**2+ry**2)**(1/2)
        saz=ry/(rx**2+ry**2)**(1/2)
        az=np.acos(caz)
        if saz<0:
            az=2*np.pi-az
        R1=_rotmat(2,az)
        ay=np.acos(rab/rz) # (must live on 0<ay<pi bc polar angle)
        R2=_rotmat(1,ay)
        R3=_rotmat(2,deg/180.*np.pi)
        R4=_rotmat(1,-ay)
        R5=_rotmat(2,-az)
        R=np.matmult(R2,R1)
        R=np.matmult(R3,R)
        R=np.matmult(R4,R)
        R=np.matmult(R5,R)
        branch.rotate(R)
        bdf=branch.Coordinates.A
        A.loc[A['globalIdx'].isin(branchidx),['posX','posY','posZ']]=bdf[['posX','posY','posZ']]
        TC.translate(O)

    def sea_of(self,idx):
        clu=self.atoms_w_same_attribute_as(find_dict={'globalIdx':idx},
                                                same_attribute='sea_idx',
                                                return_attribute='globalIdx')
        return list(clu)

MoleculeDict = dict[str,Molecule]
MoleculeList = list[Molecule]