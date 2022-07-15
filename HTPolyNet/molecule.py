from itertools import chain, product
import os
from copy import deepcopy
import pandas as pd
import numpy as np
import logging

from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.coordinates import _dfrotate
from HTPolyNet.ambertools import GAFFParameterize
import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.gromacs import grompp_and_mdrun, mdp_library

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

def get_base_reactants(mname,mdict):
    mol=mdict[mname]
    if not mol.generator:
        return [mol]
    mine=[]
    for r in mol.generator.reactants.values():
        mine.extend(get_base_reactants(r,mdict))
    return mine

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
        Aidx,Aresid,Aoresids=self.get_atom_globalIndx(A,mol,moldict)
        Bidx,Bresid,Boresids=self.get_atom_globalIndx(B,mol,moldict)
        # logging.debug(f'     A: {Aidx},{Aresid}  B: {Bidx},{Bresid}')
        return (Aidx,Bidx),(Aresid,Bresid),(Aoresids,Boresids),(A['atom'],B['atom'])

    def get_atom_globalIndx(self,A,mol,moldict):
        # logging.debug(f'get_atom_globalIndx moldict {list(moldict.keys())}')
        aname=A['atom']
        aN=A['reactant']
        areactantname=self.reactants[aN]
        aresid_inreactant=A['resid']
        # logging.debug(f'determining globalIdx of aname {aname} in reactant {aN}({areactantname}):{aresid_inreactant}')
        # generate the apparent resid in the composite molecule by visiting reactants in order
        resids={}
        cresid=0
        for k,v in self.reactants.items():
            resids[k]=[]
            for i in range(1,len(moldict[v].sequence)+1):
                resids[k].append(i+cresid)
            cresid+=len(moldict[v].sequence)
        aresid=0
        for k,v in self.reactants.items():
            if k==aN:
                aresid+=aresid_inreactant
                break
            else:
                # logging.debug(f'preseq {v} {moldict[v].sequence}')
                aresid+=len(moldict[v].sequence)
        resids[aN].remove(aresid)
        # logging.debug(f'reaction {self.name}: asking for globalIdx of resNum {aresid} atomName {aname}')
        idx=mol.TopoCoord.get_gro_attribute_by_attributes('globalIdx',{'resNum':aresid,'atomName':aname})
        return idx,aresid,resids[aN]

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

class ReactionBond:
    def __init__(self,idx,names,resids,order,bystanders,oneaways):
        self.idx=idx
        self.names=names
        self.resids=resids
        self.bystanders=bystanders
        self.oneaways=oneaways
        self.order=order
    def __str__(self):
        return f'bond {self.idx} order {self.order} bystander-resids {self.bystanders} oneaway-resids {self.oneaways}'

ReactionBondList = list[ReactionBond]

class Molecule:
    def __init__(self,name='',generator:Reaction=None):
        self.name=name
        self.TopoCoord=TopoCoord()
        self.generator:Reaction=generator
        self.sequence=[]
        self.origin=None
        self.reaction_bonds:ReactionBondList=[]
        self.symmetry_relateds=[]
        self.stereocenters=[]
        self.stereoisomers=[]
        #self.chains=[]
        self.zrecs=[]

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
        # TC.set_gro_attribute('sea-idx',-1)
        for att in ['sea-idx','chain','chain-idx','cycle','cycle-idx']:
            TC.set_gro_attribute(att,-1)
        # set symmetry class indices
        sea_idx=1
        logging.debug(f'Initialize grx {self.name}: symmetry_relateds {self.symmetry_relateds}')
        for s in self.symmetry_relateds:
            for a,b in s:
                TC.set_gro_attribute_by_attributes('sea-idx',sea_idx,{'globalIdx':a})
                TC.set_gro_attribute_by_attributes('sea-idx',sea_idx,{'globalIdx':b})
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

        # set chain, chain-idx
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
                        logging.warning(f'In molecule {self.name}, cannot identify head and tail atoms in reactive double bond\nAssuming {j} is head and {i} is tail')
                        entry=[j,i]
                    # logging.debug(f'Adding {entry} to chainlist of {self.name}')
                    TC.idx_lists['chain'].append(entry)
        TC.reset_grx_attributes_from_idx_list('chain')
        # set cycle, cycle-idx
        self.initialize_molecule_cycles()

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

    # def calculate_sea(self,sea_thresh=0.1,sea_temperature=1000,sea_nsteps=50000):
    #     ''' use a hot gromacs run to establish symmetry-equivalent atoms '''
    #     n=self.name
    #     boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
    #     mdp_prefix=mdp_library['sea']
    #     pfs.checkout(f'mdp/{mdp_prefix}.mdp')
    #     mdp_modify(f'{mdp_prefix}.mdp',{'ref_t':sea_temperature,'gen-temp':sea_temperature})
    #     for ex in ['top','itp','gro']:
    #         pfs.checkout(f'molecules/parameterized/{n}.{ex}')
    #     TC=TopoCoord(topfilename=f'{n}.top',grofilename=f'{n}.gro')
    #     dihdf=TC.Topology.D['dihedrals']
    #     dihdf.loc[:,'c0']=0.0
    #     dihdf.loc[:,'c1']=0.0
    #     dihdf.loc[:,'c2']=dihdf.loc[:,'c2'].fillna(1)
    #     TC.write_top(f'{n}-noodly.top')
    #     logging.info(f'Hot md running...output to {n}-sea')
    #     grompp_and_mdrun(gro=f'{n}',top=f'{n}-noodly',
    #                     mdp=mdp_prefix,out=f'{n}-sea',nsteps=sea_nsteps,boxSize=boxsize)
    #     sea_srs=analyze_sea(f'{n}-sea',thresh=sea_thresh)
    #     self.set_gro_attribute('sea-idx',sea_srs)

    def minimize(self,outname='',**kwargs):
        if outname=='':
            outname=f'{self.name}'
        n=self.name
        boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
        mdp_prefix=mdp_library['minimize-single-molecule']
        pfs.checkout(f'mdp/{mdp_prefix}.mdp')
        if 'checkout_required' in kwargs:
            for ex in ['top','itp','gro']:
                pfs.checkout(f'molecules/parameterized/{n}.{ex}')
        grompp_and_mdrun(gro=f'{n}',top=f'{n}',
                        mdp=mdp_prefix,out=f'{outname}',boxSize=boxsize)
        self.TopoCoord.read_gro(f'{n}.gro')

    def generate(self,outname='',available_molecules={},**kwargs):
        # logging.info(f'Generating {self.name}.mol2 for parameterization')
        if outname=='':
            outname=f'{self.name}'
        if self.generator:
            R=self.generator
            self.TopoCoord=TopoCoord()
            # assert type(R)==Reaction,'HTPolyNet only recognizes Reaction-type generators at the moment'
            logging.info(f'Using reaction {R.name} to generate {self.name}.mol2.')
            # this molecule is to be generated using a reaction
            # check to make sure this reactions reactants are among the available molecules
            # can_react=all([a in available_molecules for a in R.reactants.values()])
            # logging.debug(f'can_react: {can_react}')
            # assert can_react,f'Cannot generate {self.name} because required reactants have not been generated'
            # composite_mol=Molecule()
            shifts=[(0,0,0)]  # atom, bond, resid
            for ri in R.reactants.values():
                # logging.debug(f'adding {available_molecules[ri].name} to composite:\n{available_molecules[ri].TopoCoord.Coordinates.A.to_string()}')
                new_reactant=deepcopy(available_molecules[ri])
                shifts.append(self.TopoCoord.merge(new_reactant.TopoCoord))
            # self.TopoCoord=deepcopy(composite_mol.TopoCoord)
            #self.chains=deepcopy(composite_mol.chains)
            # logging.debug(f'composite mol\n{self.TopoCoord.Coordinates.A.to_string()}')
            self.set_sequence()
            self.prepare_new_bonds(available_molecules=available_molecules)
            # self.set_reaction_bonds(available_molecules=available_molecules)
            # logging.debug(f'Generation of {self.name}: composite molecule has {len(self.sequence)} resids')
            # logging.debug(f'generation of {self.name}: composite molecule:\n{composite_mol.TopoCoord.Coordinates.A.to_string()}')
            idx_mapper=self.make_bonds()
            # self.update_chains(idx_mapper)
            self.TopoCoord.set_gro_attribute('reactantName',R.product)
            self.write_gro_attributes(['z','nreactions','reactantName','sea-idx','cycle','cycle-idx','chain','chain-idx'],f'{R.product}.grx')
            # nrb=[]
            # for b in self.reaction_bonds:
            #     at,rr,nn=b
            #     i,j=atmake_bonds
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
        self.prepare_new_bonds(available_molecules=available_molecules)
        self.TopoCoord.set_gro_attribute('reactantName',reactantName)
        if not self.generator:
            self.initialize_monomer_grx_attributes()
            self.write_gro_attributes(['z','nreactions','reactantName','sea-idx','cycle','cycle-idx','chain','chain-idx'],f'{reactantName}.grx')
        else:
            grx=f'{reactantName}.grx'
            if (os.path.exists(grx)):
                self.TopoCoord.read_gro_attributes(grx)
                #self.reset_chains_from_attributes()
        # logging.debug(f'{self.name} gro\n{self.TopoCoord.Coordinates.A.to_string()}')

    # def set_reaction_bonds(self,available_molecules={}):
    def prepare_new_bonds(self,available_molecules={}):
        # logging.debug(f'set_reaction_bonds: molecules {list(available_molecules.keys())}')
        if len(self.reaction_bonds)>0:
            # assume prepare_new_bonds was called during construction
            return
        R=self.generator
        self.reaction_bonds=[]
        TC=self.TopoCoord
        # logging.debug(f'prepare_new_bonds {self.name}: chainlists {TC.idx_lists["chain"]}')
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
            sequence_residue_idx_origins=[0,0]
            if len(reactant_sequences)==2:
                sequence_residue_idx_origins[1]=len(reactant_sequences[0])
            in_product_resids=[in_reactant_resids[x]+sequence_residue_idx_origins[x] for x in [0,1]]
            atom_idx=[TC.get_gro_attribute_by_attributes('globalIdx',{'resNum':in_product_resids[x],'atomName':atom_names[x]}) for x in [0,1]]
            # logging.debug(f'{R.name} names {atom_names} in_product_resids {in_product_resids} idx {atom_idx}')
            bystanders=[[],[]]
            for x in [0,1]:
                atom_bystanders=self.TopoCoord.interresidue_partners_of(atom_idx[x])
                if atom_idx[1-x] in atom_bystanders: # if new bond has not yet formed
                    atom_bystanders.remove(atom_idx[1-x])
                bystanders[x]=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':y}) for y in atom_bystanders]
            chains=[TC.get_gro_attribute_by_attributes('chain',{'globalIdx':x}) for x in atom_idx]
            chain_idx=[TC.get_gro_attribute_by_attributes('chain-idx',{'globalIdx':x}) for x in atom_idx]
            oneaways=[None,None]
            # assert chain_idx[0]==0 or chain_idx[1]==0 # one must be a head
            if chains != [-1,-1]:
                # logging.debug(f'finding oneaways: chains: {chains} chain_idx: {chain_idx}')
                # intraresidue_chainmates=[TC.get_gro_attribute_by_attributes('globalIdx',{'chain':chains[x],'resNum':in_product_resids[x]}) for x in [0,1]]
                chain_idx_lists=[TC.idx_lists['chain'][chains[x]] for x in [0,1]]
                # logging.debug(f'chain_idx_lists: {chain_idx_lists}')
                if chain_idx[0]==0 or chain_idx[0]>chain_idx[1]:
                    # find next member of chainlist that is *not* in *this* residue
                    # logging.debug(f'querying resids from {chain_idx_lists[0][1:]}')
                    for idx in chain_idx_lists[0][chain_idx[0]+1:]:
                        nresid=TC.get_gro_attribute_by_attributes('resNum',{'globalIdx':idx})
                        if nresid!=in_product_resids[0]:
                            oneaways[0]=nresid
                    # logging.debug(f'querying resids from {chain_idx_lists[1][chain_idx[1]-1::-1]}')
                    for idx in chain_idx_lists[1][chain_idx[1]-1::-1]:
                        nresid=TC.get_gro_attribute_by_attributes('resNum',{'globalIdx':idx})
                        if nresid!=in_product_resids[1]:
                            oneaways[1]=nresid
                elif chain_idx[1]==0 or chain_idx[1]>chain_idx[0]:
                    for idx in chain_idx_lists[1][chain_idx[1]+1:]:
                        nresid=TC.get_gro_attribute_by_attributes('resNum',{'globalIdx':idx})
                        if nresid!=in_product_resids[1]:
                            oneaways[1]=nresid
                    for idx in chain_idx_lists[0][chain_idx[1]::-1]:
                        nresid=TC.get_gro_attribute_by_attributes('resNum',{'globalIdx':idx})
                        if nresid!=in_product_resids[0]:
                            oneaways[0]=nresid
            # logging.debug(f'prepare_new_bonds {self.name} oneaways {oneaways}')
            self.reaction_bonds.append(ReactionBond(atom_idx,atom_names,in_product_resids,order,bystanders,oneaways))

            # (Aidx,Bidx),(aresid,bresid),(Aoresids,Boresids),(Aname,Bname)=R.get_bond_atom_globalIdx(bond,self,available_molecules)
            # rb.bond=(Aidx,Bidx)
            # rb=ReactionBond([Aidx,Bidx],)

            # self.reaction_bonds.append(((Aidx,Bidx),(aresid,bresid),(Aoresids,Boresids),(Aname,Bname),bond['order']))
        # outstr='\n'+"\n".join([str(x) for x in self.reaction_bonds])
        # logging.debug(f'{R.name} reaction_bonds {outstr}')

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
        # logging.debug(f'{self.name} sequence: {self.sequence}')

    def idx_mappers(self,otherTC:TopoCoord,other_bond,bystanders,oneaways):
        assert len(other_bond)==2
        assert len(bystanders)==2
        assert len(oneaways)==2
        logging.debug(f'idx_mappers begins: template name {self.name}')
        i_idx,j_idx=other_bond
        i_resName,i_resNum,i_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i_idx})
        j_resName,j_resNum,j_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j_idx})
        logging.debug(f'idx_mappers: i_idx {i_idx} i_resName {i_resName} i_resNum {i_resNum} i_atomName {i_atomName}')
        logging.debug(f'idx_mappers: j_idx {j_idx} j_resName {j_resName} j_resNum {j_resNum} j_atomName {j_atomName}')
        # identify the template bond represented by the other_bond parameter
        ij=[]
        for B in self.reaction_bonds:
            temp_resids=B.resids
            temp_iname,temp_jname=B.names
            temp_iresname=self.sequence[temp_resids[0]-1]
            temp_jresname=self.sequence[temp_resids[1]-1]
            temp_bystanders=B.bystanders
            temp_oneaways=B.oneaways
            logging.debug(f'idx_mappers: temp_iresname {temp_iresname} temp_iname {temp_iname}')
            logging.debug(f'idx_mappers: temp_jresname {temp_jresname} temp_jname {temp_jname}')
            if (i_atomName,i_resName)==(temp_iname,temp_iresname):
                ij=[0,1]
                break # found it -- stop looking
            elif (i_atomName,i_resName)==(temp_jname,temp_jresname):
                ij=[1,0]
                break
        assert len(ij)==2,f'Mappers using template {self.name} unable to map from instance bond {i_resName}-{i_resNum}-{i_atomName}---{j_resName}-{j_resNum}-{j_atomName}'
        inst_resids=[i_resNum,j_resNum]
        inst_resids=[inst_resids[ij[x]] for x in [0,1]]
        inst_bystanders=[bystanders[ij[x]] for x in [0,1]]
        inst_oneaways=[oneaways[ij[x]] for x in [0,1]]
        assert all([len(inst_bystanders[x])==len(temp_bystanders[0]) for x in [0,1]]),f'Error: bystander count mismatch'
           # use dataframe merges to create globalIdx maps
        instdf=otherTC.Coordinates.A
        tempdf=self.TopoCoord.Coordinates.A
        inst2temp={}
        temp2inst={}
        # logging.debug(f'inst resids from {[i_resNum,j_resNum,*inst_bystanders]}') 
        # logging.debug(f'temp resids from {[temp_iresid,temp_jresid,*temp_bystanders]}')
        for inst,temp in zip([*inst_resids,*inst_bystanders[0],*inst_bystanders[1],*inst_oneaways],
                             [*temp_resids,*temp_bystanders[0],*temp_bystanders[1],*temp_oneaways]):
            if inst and temp:  # None's in the bystander lists and oneaways lists should be ignored
                logging.debug(f'map inst resid {inst} to template resid {temp}')
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
                    if temp_idx in temp2inst and temp2inst[temp_idx]!=inst_idx:
                        raise Exception(f'Error: temp_idx {temp_idx} already claimed in temp2inst; bug')
                    temp2inst[temp_idx]=inst_idx
        assert len(inst2temp)==len(temp2inst),f'Error: could not establish two-way dict of atom globalIdx'
        return (inst2temp,temp2inst)

    def x_idx_mappers(self,otherTC,other_bond):
        seq_res_is_bystander=[False for _ in self.sequence]
        logging.debug(f'idx_mappers begins: template name {self.name}')
        # get system resName, resNum, and atomNames for each of the two bonded atoms
        i_idx,j_idx=other_bond
        i_resName,i_resNum,i_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i_idx})
        j_resName,j_resNum,j_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j_idx})
        logging.debug(f'idx_mappers: i_idx {i_idx} i_resName {i_resName} i_resNum {i_resNum} i_atomName {i_atomName}')
        logging.debug(f'idx_mappers: j_idx {j_idx} j_resName {j_resName} j_resNum {j_resNum} j_atomName {j_atomName}')
        # Either of atom i or j *could* bond to a third or even fourth residue that
        # is not the residue of the other.  Any such "bystander" residue is asserted
        # to be represented in the template molecule, and must therefore be included
        # along with the residues of i and j among the residues from which interactions
        # can be mapped.
        neighbors_of_i=otherTC.partners_of(i_idx)
        neighbors_of_i.remove(j_idx)
        resid_bystanders_of_ij=[]  # bystanders bonded to atom i that are not in resid of atom j
        resid_bystanders_of_iij=[] # bystanders bonded to _a neighbor of_ atom i that are not in resid of atom _i or_ j
        logging.debug(f'neighbors of i {i_idx} {i_resName} {i_resNum} {i_atomName}:')
        for xx in neighbors_of_i:
            neighbors_of_in=otherTC.partners_of(xx)
            # neighbors_of_in.remove(i_idx)
            x_resName,x_resNum,x_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xx})
            # logging.debug(f'-> {xx} {x_resName} {x_resNum} {x_atomName} n {neighbors_of_in}')
            if x_resNum!=i_resNum and x_resNum!=j_resNum and (x_resNum,x_resName) not in resid_bystanders_of_ij:
                resid_bystanders_of_ij.append((x_resNum,x_resName))
            for xxx in neighbors_of_in:
                xx_resName,xx_resNum,xx_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xxx})
                # logging.debug(f'   -> {xxx} {xx_resName} {xx_resNum} {xx_atomName}')
                if xx_resNum!=i_resNum and xx_resNum!=j_resNum and (xx_resNum,xx_resName) not in resid_bystanders_of_iij:
                    resid_bystanders_of_iij.append((xx_resNum,xx_resName))
        neighbors_of_j=otherTC.partners_of(j_idx)
        neighbors_of_j.remove(i_idx)
        # logging.debug(f'neighbors of j {j_idx} {j_resName} {j_resNum} {j_atomName}:')
        resid_bystanders_of_ji=[]  # bystanders bonded to atom j that are not in resid of atom i
        resid_bystanders_of_jji=[] # bystanders bonded to _a neighbor of_ atom j that are not in resid of atom i _or j_
        for xx in neighbors_of_j:
            neighbors_of_jn=otherTC.partners_of(xx)
            # neighbors_of_jn.remove(j_idx)
            x_resName,x_resNum,x_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xx})
            # logging.debug(f'-> {xx} {x_resName} {x_resNum} {x_atomName} n {neighbors_of_jn}')
            if x_resNum!=i_resNum and x_resNum!=j_resNum and (x_resNum,x_resName) not in resid_bystanders_of_ji:
                resid_bystanders_of_ji.append((x_resNum,x_resName))
                # logging.debug(f'{xx} {x_resName} {x_resNum} {x_atomName}')
            for xxx in neighbors_of_jn:
                xx_resName,xx_resNum,xx_atomName=otherTC.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':xxx})
                # logging.debug(f'   -> {xxx} {xx_resName} {xx_resNum} {xx_atomName}')
                if xx_resNum!=i_resNum and xx_resNum!=j_resNum and (xx_resNum,xx_resName) not in resid_bystanders_of_jji:
                    resid_bystanders_of_jji.append((xx_resNum,xx_resName))

        # logging.debug(f'idx_mappers: resid_bystanders_of_ij {resid_bystanders_of_ij}')
        # logging.debug(f'idx_mappers: resid_bystanders_of_iij {resid_bystanders_of_iij}')
        # logging.debug(f'idx_mappers: resid_bystanders_of_ji {resid_bystanders_of_ji}')
        # logging.debug(f'idx_mappers: resid_bystanders_of_jji {resid_bystanders_of_jji}')

        resid_bystanders=[*resid_bystanders_of_ij,*resid_bystanders_of_iij,*resid_bystanders_of_ji,*resid_bystanders_of_jji]
        # logging.debug(f'idx_mappers: other_bond {other_bond} resid_bystanders {resid_bystanders}')
        temp_bystanders=[]
        inst_bystanders=[]
        # for xx in resid_bystanders:
            # rn,rname=xx
            # if not rname in self.sequence:
            #     logging.error(f'secondary neighbor {rn} {rname}: no res in pattern sequence {self.sequence}')
            #     raise Exception('this is a bug')
            # for i,rnm in enumerate(self.sequence):
            #     ib=seq_res_is_bystander[i]
            #     if rnm==rname and not ib: # won't work if all resnames are same!!!
            #         seq_res_is_bystander[i]=True
            #         temp_bystanders.append(i+1) # resids start at 1 not 0!!!
            #         inst_bystanders.append(rn)
            #         break
            # else:
            #     logging.error(f'secondary neighbor {rn} {rname}: no available res in pattern sequence {self.sequence} {seq_res_is_bystander}')

        temp2inst={}
        inst2temp={}
        temp_iresid=-1
        temp_jresid=-1
        # identify the template bond represented by the other_bond parameter
        for b in self.reaction_bonds:
            (Aidx,Bidx),(aresid,bresid),(Aoresids,Boresids),(Aname,Bname),order=b
            Aresname=self.sequence[aresid-1]
            Bresname=self.sequence[bresid-1]
            # logging.debug(f'idx_mappers: reaction_bond {b}')
            # logging.debug(f'idx_mappers: {Aresname} {aresid} {Bresname} {bresid}')
            if (i_atomName,i_resName)==(Aname,Aresname):
                temp_iresid=aresid
                temp_jresid=bresid
                # logging.debug(f'idx_mappers: temp_iresid {temp_iresid} temp_jresid {temp_jresid}')
                break # found it -- stop looking
            elif (i_atomName,i_resName)==(Bname,Bresname):
                temp_iresid=bresid
                temp_jresid=aresid
                # logging.debug(f'idx_mappers: temp_iresid {temp_iresid} temp_jresid {temp_jresid}')
                break
        if temp_iresid==-1:
            logging.error(f'Mappers using template {self.name} unable to map from instance bond {i_resName}-{i_resNum}-{i_atomName}---{j_resName}-{j_resNum}-{j_atomName}')
            raise Exception

        iapp_inst_bystanders=list(set([*resid_bystanders_of_ij,*resid_bystanders_of_iij]))
        assert len(iapp_inst_bystanders)==len(Aoresids)
        # logging.debug(f'idx_mappers: iapp_inst_bystanders {iapp_inst_bystanders} Aoresids {Aoresids}')
        for ib,tb in zip(iapp_inst_bystanders,Aoresids):
            inst_bystanders.append(ib[0])
            temp_bystanders.append(tb)
        japp_inst_bystanders=list(set([*resid_bystanders_of_ji,*resid_bystanders_of_jji]))
        assert len(japp_inst_bystanders)==len(Boresids),f'idx_mappers: japp_inst_bystanders {japp_inst_bystanders} Boresids {Boresids}'
        # logging.debug(f'idx_mappers: japp_inst_bystanders {japp_inst_bystanders} Boresids {Boresids}')
        for ib,tb in zip(japp_inst_bystanders,Boresids):
            inst_bystanders.append(ib[0])
            temp_bystanders.append(tb)
        # logging.debug(f'idx_mappers: inst_bystanders {inst_bystanders}')
        # logging.debug(f'idx_mappers: temp_bystanders {temp_bystanders}')

        # use dataframe merges to create globalIdx maps
        instdf=otherTC.Coordinates.A
        tempdf=self.TopoCoord.Coordinates.A
        inst2temp={}
        temp2inst={}
        # logging.debug(f'inst resids from {[i_resNum,j_resNum,*inst_bystanders]}')
        # logging.debug(f'temp resids from {[temp_iresid,temp_jresid,*temp_bystanders]}')
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
                if temp_idx in temp2inst and temp2inst[temp_idx]!=inst_idx:
                    raise Exception(f'Error: temp_idx {temp_idx} already claimed in temp2inst; bug')
                temp2inst[temp_idx]=inst_idx
        # if len(inst2temp)!=len(temp2inst):
            # logging.debug(f'len(inst2temp) {len(inst2temp)} != len(temp2inst) {len(temp2inst)}')
            # logging.debug(f'inst2temp:')
            # for k,v in inst2temp.items():
                # logging.debug(f'   {k} <-> {v}')
            # logging.debug(f'temp2inst:')
            # for k,v in temp2inst.items():
                # logging.debug(f'   {k} <-> {v}')
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

    def label_cycle_atoms(self):
        self.TopoCoord.label_cycle_atoms()

    def get_resname(self,internal_resid):
        # logging.debug(f'{self.name} sequence: {self.sequence}')
        return self.sequence[internal_resid-1]

    def inherit_attribute_from_reactants(self,attribute,available_molecules,increment=True,no_increment_if_negative=True):
        adf=self.TopoCoord.Coordinates.A
        ordered_attribute_idx=[]
        curr_max=0
        # logging.debug(f'{self.name}({adf.shape[0]}) inheriting {attribute} from {self.sequence}')
        # logging.debug(f'available molecules {list(available_molecules.keys())}')
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
            # logging.debug(f'{r}->{len(x)}')
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

    def load_top_gro(self,topfilename,grofilename,mol2filename=''):
        self.TopoCoord=TopoCoord(topfilename=topfilename,grofilename=grofilename,mol2filename=mol2filename)

    # def analyze_sea_topology(self):
    #     self.TopoCoord.analyze_sea_topology()

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
        TC=self.TopoCoord
        # logging.debug(f'{self.name} make_bonds: chainlists: {TC.idx_lists["chain"]}')
        # for cl in TC.idx_lists['chain']:
        #     nl=[self.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl]
        #     logging.debug(f'{nl}')
        for i,B in enumerate(self.reaction_bonds):
            # logging.debug(f'{self.name} reaction bond {i}: {str(B)} ')
            aidx,bidx=B.idx
            aname,bname=B.names
            aresid,bresid=B.resids
            bystanders=B.bystanders
            oneaways=B.oneaways
            order=B.order
            # (aidx,bidx),(aresid,bresid),(aoresids,boresids),(aname,bname),order=B
            logging.debug(f'generating {self.name} bond {i} {aresid}:{aname}:{aidx}-{bresid}:{bname}:{bidx} order {order}')
            bonds.append((aidx,bidx,order))
            if aresid!=bresid:
                # transrot identifies the two sacrificial H's
                hxi,hxj=self.transrot(aidx,aresid,bidx,bresid,connected_resids=[oneaways[1],*bystanders[1]])
                skip_H.append(i)
                hs_from_tr.append(hxi)
                hs_from_tr.append(hxj)
        # make_bonds returns list of sacrificial H idx's derived
        # from bonds whose indices are not in skip_H
        idx_scratch=self.TopoCoord.make_bonds(bonds,skip_H=skip_H)
        # the H's identified in transrot plus those identified by make_bonds
        # must be deleted
        idx_scratch.extend(hs_from_tr)
        idx_mapper=self.TopoCoord.delete_atoms(idx_scratch)
        #TC.remap_idx_list('chain',idx_mapper)
        for i,B in enumerate(self.reaction_bonds):
            aidx,bidx=B.idx
            # (aidx,bidx),(aresid,bresid),(aoresids,boresids),(aname,bname),order=B
            aidx=idx_mapper[aidx]
            bidx=idx_mapper[bidx]
            self.TopoCoord.decrement_gro_attribute_by_attributes('z',{'globalIdx':aidx})
            self.TopoCoord.decrement_gro_attribute_by_attributes('z',{'globalIdx':bidx})
            self.TopoCoord.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':aidx})
            self.TopoCoord.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':bidx})
        self.initialize_molecule_cycles()
        return idx_mapper

    def transrot(self,at_idx,at_resid,from_idx,from_resid,connected_resids=[]):
        # Rotate and translate
        # logging.debug('transrot for building {self.name} from {at_resid}:{at_idx} -- {from_resid}:{from_idx}')
        if at_resid==from_resid:
            return
        TC=self.TopoCoord
        ATC=TopoCoord()
        BTC=TopoCoord()
        C=TC.gro_DataFrame('atoms')
        ATC.Coordinates.A=C[C['resNum']==at_resid].copy()
        bresids=connected_resids[:]
        bresids.append(from_resid)
        BTC.Coordinates.A=C[C['resNum'].isin(bresids)].copy()
        mypartners=TC.partners_of(at_idx)
        otpartners=TC.partners_of(from_idx)
        # logging.debug(f'Partners of {at_idx} {mypartners}')
        # logging.debug(f'Partners of {from_idx} {otpartners}')
        myHpartners={k:v for k,v in zip(mypartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners={k:v for k,v in zip(otpartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        myHighestH={k:v for k,v in myHpartners.items() if v==max([k for k in myHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
        otHighestH={k:v for k,v in otHpartners.items() if v==max([k for k in otHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
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

    def flip_stereocenter(self,idx):
        """flip_stereocenter Flips stereochemistry of atom at idx

        :param idx: global index of atom
        :type idx: int
        """
        TC=self.TopoCoord
        A=TC.Coordinates.A
        ligand_idx=self.TopoCoord.Topology.bondlist.partners_of(idx)
        if len(ligand_idx)!=4:
            logging.debug(f'Atom {idx} cannot be a stereocenter; it only has {len(ligand_idx)} ligands.')
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
                                                same_attribute='sea-idx',
                                                return_attribute='globalIdx')
        return list(clu)