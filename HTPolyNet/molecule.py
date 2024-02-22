"""

.. module:: molecule
   :synopsis: manages generation of molecular templates
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import pandas as pd
import numpy as np
import logging
import shutil
from itertools import product
from copy import deepcopy

import HTPolyNet.projectfilesystem as pfs
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.bondtemplate import BondTemplate,BondTemplateList,ReactionBond,ReactionBondList
from HTPolyNet.coordinates import GRX_ATTRIBUTES
from HTPolyNet.ambertools import GAFFParameterize
from HTPolyNet.gromacs import mdp_modify,gro_from_trr
from HTPolyNet.command import Command
from HTPolyNet.reaction import Reaction, ReactionList, reaction_stage, generate_product_name, reactant_resid_to_presid

logger=logging.getLogger(__name__)

def _rotmat(axis,radians):
    """_rotmat generates a rotation matrix that rotates coordinates around axis by radians

    :param axis: 0=x,1=y,2=z
    :type axis: int
    :param radians: angle measure
    :type radians: float
    :return: rotation matrix
    :rtype: np.ndarray(3,float)
    """
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

def yield_bonds(R:Reaction,TC:TopoCoord,resid_mapper):
    """yield_bonds for each bond pattern in the reaction R, yield the specific bond represented by this pattern in the molecule's topology

    :param R: a Reaction
    :type R: Reaction
    :param TC: molecule topology and coordinates
    :type TC: TopoCoord
    :param resid_mapper: dictionary that maps in-reactant resids to in-product resids
    :type resid_mapper: dict
    :yield: a reaction bond
    :rtype: ReactionBond
    """
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

def yield_bonds_as_df(R:Reaction,TC:TopoCoord,resid_mapper):
    """yield_bonds_as_df returns a pandas DataFrame identifying all reaction bonds obtained by matching the reaction template bonds to instances in the molecule topology

    :param R: a Reaction
    :type R: Reaction
    :param TC: molecule topology and coordinates
    :type TC: TopoCoord
    :param resid_mapper: dictionary that maps in-reactant resids to in-product resids
    :type resid_mapper: dict
    :return: dataframe with all reaction bonds
    :rtype: pd.DataFrame
    """
    bdf=pd.DataFrame()
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
        ai,aj=[TC.get_gro_attribute_by_attributes('globalIdx',{'resNum':in_product_resids[x],'atomName':atom_names[x]}) for x in [0,1]]
        ri,rj=in_product_resids
        row_dict={'ai':[ai],'aj':[aj],'ri':[ri],'rj':[rj],'order':[order],'reactantName':[R.product]}
        bdf=pd.concat((bdf,pd.DataFrame(row_dict)),ignore_index=True)
    return bdf

class Molecule:
    def __init__(self,name='',generator:Reaction=None,origin:str=None):
        self.name=name
        self.parentname=name # stereoisomer parent
        self.TopoCoord=TopoCoord()
        self.generator:Reaction=generator
        self.sequence=[]
        self.origin=origin
        self.reaction_bonds:ReactionBondList=[]
        self.bond_templates:BondTemplateList=[]
        self.symmetry_relateds=[]
        self.stereocenters=[] # list of atomnames
        self.stereoisomers={}
        self.nconformers=0
        self.conformers_dict={}
        self.conformers=[] # just a list of gro file basenames
        self.zrecs=[]
        self.is_reactant=False

    @classmethod
    def New(cls,mol_name,generator:Reaction,molrec={}):
        """New generates a new, partially populated Molecule based on directives in the configuration input

        :param mol_name: name of molecule
        :type mol_name: str
        :param generator: reaction that generates this molecule, if applicable
        :type generator: Reaction
        :param molrec: dictionary of directive for this molecule, defaults to {}
        :type molrec: dict, optional
        :return: a new Molecule object
        :rtype: Molecule
        """
        M=cls(name=mol_name)
        M.generator=generator
        if not molrec: return M
        M.symmetry_relateds=molrec.get('symmetry_equivalent_atoms',[])
        M.stereocenters=molrec.get('stereocenters',[])
        # expand list of stereocenters if any are in symmetry sets
        extra_stereocenters=[]
        for stc in M.stereocenters:
            for sc in M.symmetry_relateds:
                if stc in sc:
                    sc_copy=sc.copy()
                    sc_copy.remove(stc)
                    if not sc_copy in extra_stereocenters:
                        extra_stereocenters.extend(sc_copy)
        M.stereocenters.extend(extra_stereocenters)
        logger.debug(f'{M.name} stereocenters: {M.stereocenters}')
        # generate shells for new stereoisomers
        M.create_new_stereoisomers()
        logger.debug(f'{M.name} stereoisomers: {[s.name for s in M.stereoisomers.values()]}')
        M.conformers_dict=molrec.get('conformers',{})
        logger.debug(f'{M.name} conformers_dict {M.conformers_dict}')
        return M

    def set_origin(self,value):
        """set_origin sets the value of the origin member

        :param value: value
        :type value: anything
        """
        self.origin=value

    def get_origin(self):
        """get_origin returns the value of the origin member

        :return: value of origin member
        :rtype: anything
        """
        return self.origin

    def update_zrecs(self,zrecs,moldict):
        """update_zrecs updates the "z-records" based on z's declared in the input configuration file

        :param zrecs: zrecs extracted from configuration file for this molecule
        :type zrecs: dict
        :param moldict: dictionary of available molecules
        :type moldict: dicts
        """
        def replace_if_greater(rec,D,matchattr=[],maxattr=[]):
            if not matchattr or not maxattr: return
            for r in D:
                matched=all([r[a]==rec[a] for a in matchattr])
                if matched:
                    replace=all([r[a]<rec[a] for a in maxattr])
                    if replace:
                        D.remove(r)
                        D.append(rec)
                    return True
            return False
        logger.debug(f'Update zrecs in {self.name} from {zrecs}')
        seq=self.sequence
        for zr in zrecs:
            resid=zr['resid']-1
            rname=seq[resid]
            target=moldict[rname]
            logger.debug(f'{target.name} {target.zrecs} ->')
            found=replace_if_greater(zr,target.zrecs,matchattr=['resid','atom'],maxattr=['z'])
            if not found: target.zrecs.append(zr)
            logger.debug(f'-> {target.name} {target.zrecs}')

    def determine_sequence(self,moldict):
        """determine_sequence recursively determine the sequence of a molecule using the network of reactions that must be executed to generate it from primitives

        :param moldict: dictionary of available molecules
        :type moldict: dict
        :return: list of resnames in order of sequence
        :rtype: list
        """
        if not self.generator: return [self.parentname]
        R:Reaction=self.generator
        thisseq=[]
        for rid,rname in R.reactants.items():
            parentname=moldict[rname].parentname
            thisseq.extend(moldict[parentname].determine_sequence(moldict))
            # logger.debug(thisseq)
        return thisseq
    
    def set_sequence_from_moldict(self,moldict):
        """set_sequence_from_moldict set the sequence of this molecule using the recursive determine_sequence method

        :param moldict: dictionary of available molecules
        :type moldict: dict
        :return: self
        :rtype: Molecule
        """
        self.sequence=self.determine_sequence(moldict)
        return self

    def set_sequence_from_coordinates(self):
        """set_sequence Establish the sequence-list (residue names in order) based on resNum attributes in atom list
        """
        adf=self.TopoCoord.gro_DataFrame('atoms')
        trial_sequence=[]
        current_resid=0
        for i,r in adf.iterrows():
            ri=r['resNum']
            rn=r['resName']
            if ri!=current_resid:
                current_resid=ri
                trial_sequence.append(rn)
        assert trial_sequence==self.sequence,f'trial {trial_sequence} seq {self.sequence}'
        return self

    def create_new_stereoisomers(self):
        """create_new_stereoisomers generate new molecules to hold stereoisomers of self

        :return: None if no action taken
        :rtype: none
        """
        if self.generator: return  # we only consider stereoisomers on monomers
        if not self.stereocenters: return
        basename=self.name+'-S'
        b=[[0,1] for _ in range(len(self.stereocenters))]
        sseq=product(*b)
        next(sseq) # skip the unmodified; it's the base molecule
        for x in sseq:
            s=''.join([str(_) for _ in x])
            mname=f'{basename}{s}'
            self.stereoisomers[mname]=Molecule.New(mname,None)
            self.stereoisomers[mname].parentname=self.name
            
    def initialize_molecule_rings(self):
        """initialize_molecule_rings generates the dictionary of rings

        the dictionary of rings is keyed on ring size
        """
        TC=self.TopoCoord
        TC.Topology.detect_rings()
        logger.debug(f'Detected {len(TC.Topology.rings)} unique rings.')
        logger.debug('Done')

    def initialize_monomer_grx_attributes(self):
        """initialize_monomer_grx_attributes initializes all GRX attributes of atoms in molecule
        """
        logger.debug(f'{self.name}')
        TC=self.TopoCoord
        TC.set_gro_attribute('z',0)
        TC.set_gro_attribute('nreactions',0)
        TC.set_gro_attribute('molecule',1)
        TC.set_gro_attribute('molecule_name',self.name)
        for att in ['sea_idx','bondchain','bondchain_idx']:
            TC.set_gro_attribute(att,-1)
        # set symmetry class indices
        sea_idx=1
        logger.debug(f'{self.name}: symmetry_relateds {self.symmetry_relateds}')
        for s in self.symmetry_relateds:
            # logger.debug(f'sea_idx {sea_idx} set for set {s}')
            for atomName in s:
                # logger.debug(f'{atomName} {sea_idx}')
                TC.set_gro_attribute_by_attributes('sea_idx',sea_idx,{'atomName':atomName})
            sea_idx+=1
        # set z and nreactions
        idx=[]
        for zr in self.zrecs:
            an=zr['atom']
            rnum=zr['resid']
            if rnum!=1: continue
            z=zr['z']
            logger.debug(f'{self.name} setting z for {an} {rnum} {z}')
            TC.set_gro_attribute_by_attributes('z',z,{'atomName':an,'resNum':rnum})
            idx.append(TC.get_gro_attribute_by_attributes('globalIdx',{'atomName':an,'resNum':rnum}))
            for sr in self.symmetry_relateds:
                # logger.debug(f'{self.name}: setting z for {an}, considering sr {sr}')
                if an in sr:
                    for bn in sr:
                        if bn==an: continue
                        # logger.debug(f'{self.name}: setting z for {bn}')
                        idx.append(TC.get_gro_attribute_by_attributes('globalIdx',{'atomName':bn,'resNum':rnum}))
                        TC.set_gro_attribute_by_attributes('z',z,{'atomName':bn,'resNum':rnum})

        TC.idx_lists['bondchain']=[]
        pairs=product(idx,idx)
        for i,j in pairs:
            if i<j:
                iname=TC.get_gro_attribute_by_attributes('atomName',{'globalIdx':i})
                jname=TC.get_gro_attribute_by_attributes('atomName',{'globalIdx':j})
                if TC.are_bonded(i,j) and iname.startswith('C') and jname.startswith('C'):
                    # this monomer has two carbon atoms capable of reacting
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
                        logger.warning(f'In molecule {self.name}, cannot identify bonded reactive head and tail atoms\nAssuming {j} is head and {i} is tail')
                        entry=[j,i]
                    # logger.debug(f'Adding {entry} to chainlist of {self.name}')
                    TC.idx_lists['bondchain'].append(entry)
        TC.reset_grx_attributes_from_idx_list('bondchain')
        self.initialize_molecule_rings()

    def previously_parameterized(self):
        """previously_parameterized if a gro file exists in the project molecule/parameterized directory for this molecule, return True

        :return: True if gro file found
        :rtype: bool
        """
        rval=True
        for ext in ['gro']:
            rval=rval and pfs.exists(os.path.join('molecules/parameterized',f'{self.name}.{ext}'))
        return rval

    def parameterize(self,outname='',input_structure_format='mol2',**kwargs):
        """parameterize manages GAFF parameterization of this molecule

        :param outname: output file basename, defaults to ''
        :type outname: str, optional
        :param input_structure_format: input structure format, defaults to 'mol2' ('pdb' is other possibility)
        :type input_structure_format: str, optional
        """
        assert os.path.exists(f'{self.name}.{input_structure_format}'),f'Cannot parameterize molecule {self.name} without {self.name}.{input_structure_format} as input'
        if outname=='':
            outname=f'{self.name}'
        GAFFParameterize(self.name,outname,input_structure_format=input_structure_format,**kwargs)
        self.load_top_gro(f'{outname}.top',f'{outname}.gro',mol2filename=f'{outname}.mol2',wrap_coords=False)
        self.initialize_molecule_rings()
        self.TopoCoord.write_tpx(f'{outname}.tpx')

    def minimize(self,outname='',**kwargs):
        """minimize manages invocation of vacuum minimization

        :param outname: output file basename, defaults to ''
        :type outname: str, optional
        """
        if outname=='':
            outname=f'{self.name}'
        self.TopoCoord.vacuum_minimize(outname,**kwargs)

    def relax(self,relax_dict):
        """relax manages invocation of MD relaxations

        :param relax_dict: dictionary of simulation directives
        :type relax_dict: dict
        """
        deffnm=relax_dict.get('deffnm',f'{self.name}-relax')
        nsteps=relax_dict.get('nsteps',10000)
        temperature=relax_dict.get('temperature',10000)
        n=self.name
        boxsize=np.array(self.TopoCoord.maxspan())+2*np.ones(3)
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix='single-molecule-nvt'
        pfs.checkout(f'mdp/{mdp_prefix}.mdp')
        mdp_modify(f'{mdp_prefix}.mdp',{'nsteps':nsteps,'gen-vel':'yes','ref_t':temperature,'gen-temp':temperature})
        logger.info(f'In vacuo equilibration of {self.name}.gro for {nsteps} steps at {temperature} K')
        self.TopoCoord.grompp_and_mdrun(out=deffnm,mdp=mdp_prefix,boxSize=boxsize)

    def center_coords(self,new_boxsize:np.ndarray=None):
        """center_coords wrapper for the TopoCoord.center_coords method

        :param new_boxsize: new box size, defaults to None
        :type new_boxsize: np.ndarray, optional
        """
        self.TopoCoord.center_coords(new_boxsize)

    def generate(self,outname='',available_molecules={},**kwargs):
        """generate manages generating topology and coordinates for self

        :param outname: output file basename, defaults to ''
        :type outname: str, optional
        :param available_molecules: dictionary of available molecules, defaults to {}
        :type available_molecules: dict, optional
        """
        logger.debug(f'{self.name}.generate() begins')
        if outname=='': outname=f'{self.name}'
        do_minimization=True
        GAFF_dict=kwargs.get('GAFF',{})
        if GAFF_dict: do_minimization=GAFF_dict.get('minimize_molecules',True)
        do_parameterization=False
        if self.generator:
            R=self.generator
            if R.stage in [reaction_stage.cure,reaction_stage.param,reaction_stage.cap]: do_parameterization=True
            self.TopoCoord=TopoCoord()
            logger.debug(f'Using reaction {R.name} ({str(R.stage)}) to generate {self.name} parent {self.parentname}')
            isf='mol2'
            resid_mapper=[]
            for ri in R.reactants.values():
                logger.debug(f'Adding {ri}')
                new_reactant=deepcopy(available_molecules[ri])
                new_reactant.TopoCoord.write_mol2(filename=f'{self.name}-reactant{ri}-prebonding.mol2',molname=self.name)
                rnr=len(new_reactant.sequence)
                shifts=self.TopoCoord.merge(new_reactant.TopoCoord)
                # for ln in self.TopoCoord.Coordinates.A.head().to_string().split('\n'): logger.debug(ln)
                resid_mapper.append({k:v for k,v in zip(range(1,rnr+1),range(1+shifts[2],1+rnr+shifts[2]))})
            # logger.debug(f'{self.name}: resid_mapper {resid_mapper}')
            # logger.debug(f'{self.TopoCoord.idx_lists}')
            # logger.debug(f'\n{self.TopoCoord.Coordinates.A.to_string()}')
            # logger.debug(f'composite prebonded molecule in box {self.TopoCoord.Coordinates.box}')
            self.TopoCoord.write_mol2(filename=f'{self.name}-prebonding.mol2',molname=self.name)
            self.set_sequence_from_coordinates()
            # bonds_to_make=list(yield_bonds(R,self.TopoCoord,resid_mapper))
            bdf=yield_bonds_as_df(R,self.TopoCoord,resid_mapper)
            # logger.debug(f'Generation of {self.name}: composite molecule has {len(self.sequence)} resids')
            # logger.debug(f'generation of {self.name}: composite molecule:\n{composite_mol.TopoCoord.Coordinates.A.to_string()}')
            self.make_bonds(bdf,available_molecules,R.stage)
            # self.TopoCoord.set_gro_attribute('reactantName',R.product)
            self.TopoCoord.set_gro_attribute('sea_idx',-1) # turn off symmetry-equivalence for multimers
            self.TopoCoord.set_gro_attribute('molecule',1)
            self.TopoCoord.set_gro_attribute('molecule_name',self.name)
            self.write_gro_attributes(GRX_ATTRIBUTES,f'{R.product}.grx')
            # if pfs.exists(f'molecules/inputs/{self.name}.mol2'): # an override structure is present
            #     logger.debug(f'Using override input molecules/inputs/{self.name}.{isf} as a generator')
            #     pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
            # else:
            self.TopoCoord.write_mol2(filename=f'{self.name}.mol2',molname=self.name)
            if not do_parameterization:
                self.TopoCoord.write_gro(f'{self.name}.gro',grotitle=self.name)
                self.TopoCoord.write_top(f'{self.name}.top')
            # if pfs.exists(f'molecules/inputs/{self.name}.pdb'): # an override structure is present
            #     isf='pdb'
            #     logger.debug(f'Using override input molecules/inputs/{self.name}.{isf} as a generator')
            #     pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
        else:
            # this molecule was not assigned a generator: implies it is a monomer and must have a parameterization
            input_structure_formats=['mol2','pdb']
            isf=None
            for isf in input_structure_formats:
                if pfs.exists(f'molecules/inputs/{self.name}.{isf}'):
                    logger.debug(f'Using input molecules/inputs/{self.name}.{isf} as a generator')
                    pfs.checkout(f'molecules/inputs/{self.name}.{isf}')
                    break
            assert isf,'Error: no valid input structure file found'
            do_parameterization=True

        reactantName=self.name
        if do_parameterization:
            self.parameterize(outname,input_structure_format=isf,**kwargs)
        else:
            if self.name!=self.parentname:
                logger.info(f'Built {self.name} using topology of {self.parentname}; copying {self.parentname}.top to {self.name}.top')
                self.load_top_gro(f'{self.parentname}.top',f'{self.name}.gro',tpxfilename=f'{self.parentname}.tpx',wrap_coords=False)
                shutil.copy(f'{self.parentname}.top',f'{self.name}.top')
                shutil.copy(f'{self.parentname}.grx',f'{self.name}.grx')
                shutil.copy(f'{self.parentname}.tpx',f'{self.name}.tpx')

        if do_minimization:
            self.minimize(outname,**kwargs)
        self.set_sequence_from_coordinates()
        if not self.generator:
            self.TopoCoord.set_gro_attribute('reactantName',reactantName)
            self.initialize_monomer_grx_attributes()
            self.write_gro_attributes(GRX_ATTRIBUTES,f'{reactantName}.grx')
        else:
            if do_parameterization or self.name!=self.parentname:
                grx=f'{reactantName}.grx'
                if (os.path.exists(grx)):
                    self.TopoCoord.read_gro_attributes(grx)
                #self.reset_chains_from_attributes()
        # logger.debug(f'{self.name} gro\n{self.TopoCoord.Coordinates.A.to_string()}')
        self.prepare_new_bonds(available_molecules=available_molecules)

        # for ln in self.TopoCoord.Coordinates.A.head().to_string().split('\n'): logger.debug(ln)
        logger.info(f'{self.name}: {self.get_molecular_weight():.2f} g/mol')
        logger.debug('Done.')

    def get_molecular_weight(self):
        """get_molecular_weight returns the molecular weight of self

        :return: _description_
        :rtype: _type_
        """
        mass=self.TopoCoord.total_mass(units='gromacs') # g
        return mass

    def prepare_new_bonds(self,available_molecules={}):
        """prepare_new_bonds populates the bond templates and reaction bonds for self

        :param available_molecules: dictionary of available molecules, defaults to {}
        :type available_molecules: dict, optional
        """
        # logger.debug(f'set_reaction_bonds: molecules {list(available_molecules.keys())}')
        R=self.generator
        if not R:
            return
        self.reaction_bonds=[]
        self.bond_templates=[]
        TC=self.TopoCoord
        # logger.debug(f'prepare_new_bonds {self.name}: chainlists {TC.idx_lists["bondchain"]}')
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
            # logger.debug(f'{R.name} names {atom_names} in_product_resids {in_product_resids} idx {atom_idx}')
            bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames=TC.get_bystanders(atom_idx)
            oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames=TC.get_oneaways(atom_idx)
            # logger.debug(f'{self.name} bystanders {bystander_resids} {bystander_resnames} {bystander_atomidx} {bystander_atomnames}')
            # logger.debug(f'{self.name} oneaways {oneaway_resids} {oneaway_resnames} {oneaway_atomidx} {oneaway_atomnames}')
            self.reaction_bonds.append(ReactionBond(atom_idx,in_product_resids,order,bystander_resids,bystander_atomidx,oneaway_resids,oneaway_atomidx))
            intraresidue=in_product_resids[0]==in_product_resids[1]
            self.bond_templates.append(BondTemplate(atom_names,in_product_resnames,intraresidue,order,bystander_resnames,bystander_atomnames,oneaway_resnames,oneaway_atomnames))

    def idx_mappers(self,otherTC:TopoCoord,other_bond,bystanders,oneaways,uniq_atom_idx:set):
        """idx_mappers computes the mapping dictionary from molecule template index to instance index in the other TopoCoord

        :param otherTC: the other TopoCoord
        :type otherTC: TopoCoord
        :param other_bond: 2 atom indices of the bond in the other TopoCoord
        :type other_bond: list-like container of length 2
        :param bystanders: bystander lists, one for each reacting atom
        :type bystanders: lists
        :param oneaways: oneaways, one for each atom in the bond
        :type oneaways: list (2)
        :param uniq_atom_idx: set of unique atoms in template that must be mapped to instance
        :type uniq_atom_idx: set
        :raises Exception: if there is a buggy double-counting of one or more indexes
        :return: two-way dictionaries of index mappers instance<->template
        :rtype: tuple of two dictionaries
        """
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
            # logger.debug(f'temp_iresname {temp_iresname} temp_iname {temp_iname}')
            # logger.debug(f'temp_jresname {temp_jresname} temp_jname {temp_jname}')
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
        """get_angles_dihedrals returns copies of selections from the Topology interaction-type dataframes that contain the two atoms indicated in the bond

        :param bond: 2-element list-like container of ints
        :type bond: list-like container
        :raises Exception: dies if a NaN is found in any selection
        :return: tuple of the three dataframe selection copies for angles, dihedrals, and 1-4 pairs 
        :rtype: tuple
        """
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

    def get_resname(self,internal_resid):
        """get_resname returns the residue name at position internal_resid in the molecule's sequence

        :param internal_resid: molecule-internal residue index
        :type internal_resid: int
        :return: residue name
        :rtype: str
        """
        return self.sequence[internal_resid-1]

    # def inherit_attribute_from_reactants(self,attribute,available_molecules,increment=True,no_increment_if_negative=True):
    #     """inherit_attribute_from_reactants populate certain atom attributes in molecule from its constituent reactants

    #     :param attribute: _description_
    #     :type attribute: _type_
    #     :param available_molecules: _description_
    #     :type available_molecules: _type_
    #     :param increment: _description_, defaults to True
    #     :type increment: bool, optional
    #     :param no_increment_if_negative: _description_, defaults to True
    #     :type no_increment_if_negative: bool, optional
    #     """
    #     adf=self.TopoCoord.Coordinates.A
    #     ordered_attribute_idx=[]
    #     curr_max=0
    #     # logger.debug(f'{self.name}({adf.shape[0]}) inheriting {attribute} from {self.sequence}')
    #     # logger.debug(f'available molecules {list(available_molecules.keys())}')
    #     for i,r in enumerate(self.sequence):
    #         '''
    #         for this residue number, read the list of unique atom names
    #         '''
    #         namesinres=list(adf[adf['resNum']==(i+1)]['atomName'])
    #         '''
    #         access coordinates of standalone residue template with this name 'r' on the list of available molecules
    #         '''
    #         rdf=available_molecules[r].TopoCoord.Coordinates.A
    #         '''
    #         get the attribute values from residue template
    #         '''
    #         x=list(rdf[rdf['atomName'].isin(namesinres)][attribute])
    #         # logger.debug(f'{r}->{len(x)}')
    #         '''
    #         increment these attribute value based on residue number in this molecule
    #         '''
    #         if increment:
    #             i_x=[]
    #             for y in x:
    #                 if y>0 or (y<0 and not no_increment_if_negative):
    #                     i_x.append(y+curr_max)
    #                 else:
    #                     i_x.append(y)
    #             curr_max=max(i_x)
    #         ordered_attribute_idx.extend(i_x)
    #     assert len(ordered_attribute_idx)==adf.shape[0]
    #     adf[attribute]=ordered_attribute_idx

    def merge(self,other):
        """merge merges TopoCoord from other into self's TopoCoord

        :param other: another Molecule
        :type other: Molecule
        :return: a shift tuple (returned by Coordinates.merge())
        :rtype: tuple
        """
        shifts=self.TopoCoord.merge(other.TopoCoord)
        return shifts

    def load_top_gro(self,topfilename,grofilename,tpxfilename='',mol2filename='',**kwargs):
        """load_top_gro generate a new TopoCoord member object for this molecule by reading in a Gromacs topology file and a Gromacs gro file

        :param topfilename: Gromacs topology file
        :type topfilename: str
        :param grofilename: Gromacs gro file
        :type grofilename: str
        :param tpxfilename: extended topology file, defaults to ''
        :type tpxfilename: str, optional
        :param mol2filename: alternative coordinate mol2 file, defaults to ''
        :type mol2filename: str, optional
        """
        self.TopoCoord=TopoCoord(topfilename=topfilename,grofilename=grofilename,tpxfilename=tpxfilename,mol2filename=mol2filename,**kwargs)

    def set_gro_attribute(self,attribute,srs):
        """set_gro_attribute sets attribute of atoms to srs (drillst through to Coordinates.set_atomset_attributes())

        :param attribute: name of attribute
        :type attribute: str
        :param srs: scalar or list-like attribute values in same ordering as self.A
        :type srs: scalar or list-like
        """
        self.TopoCoord.set_gro_attribute(attribute,srs)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
        """Read attributes from file into self.TopoCoord.Coordinates.A

        :param grxfilename: name of input file
        :type grxfilename: str
        :param attribute_list: list of attributes to take, defaults to [] (take all)
        :type attribute_list: list, optional
        """
        self.TopoCoord.read_gro_attributes(grxfilename,attribute_list=attribute_list)

    def write_gro_attributes(self,attribute_list,grxfilename):
        """Writes atomic attributes to a file

        :param attributes_list: list of attributes to write
        :type attributes_list: list
        :param grxfilename: name of output file
        :type grxfilename: str
        """
        self.TopoCoord.write_gro_attributes(attribute_list,grxfilename)

    def make_bonds(self,bdf:pd.DataFrame,moldict,stage):
        """make_bonds adds new bonds to the molecule's topology and deletes any sacrificial hydrogens

        :param bdf: pandas dataframe identifying new bonds
        :type bdf: pd.DataFrame
        :param moldict: dictionary of available molecular templates
        :type moldict: dict
        :param stage: enumerated parameter indicating reaction_stage
        :type stage: reaction_stage(Enum)
        """
        TC=self.TopoCoord
        explicit_sacrificial_Hs={}
        for i,r in bdf.iterrows():
            aname,bname=[TC.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in [r.ai,r.rj]]
            logger.debug(f'generating {self.name} bond {r.ri}:{aname}:{r.ai}-{r.rj}:{bname}:{r.aj} order {r.order}')
            if r.ri!=r.rj:
                resid_sets=TC.get_resid_sets([r.ai,r.aj])
                hxi,hxj=self.transrot(r.ai,r.ri,r.aj,r.rj,connected_resids=resid_sets[1])
                explicit_sacrificial_Hs[i]=[hxi,hxj]
        if stage in [reaction_stage.cure, reaction_stage.param, reaction_stage.cap]:
            template_source='ambertools'
        else:
            template_source='internal'  # signals that a template molecule should be identified to parameterize this bond
        TC.update_topology_and_coordinates(bdf,moldict,explicit_sacH=explicit_sacrificial_Hs,template_source=template_source)
        self.initialize_molecule_rings()

    def transrot(self,at_idx,at_resid,from_idx,from_resid,connected_resids=[]):
        """transrot given a composite molecule, translate and rotate the piece downstream of the yet-to-be created bond specified by (at_idx,at_resid) and (from_idx,from_resid) to minimize steric overlaps and identify the best two sacrificial hydrogens

        :param at_idx: global index of left-hand atom in new bond
        :type at_idx: int
        :param at_resid: global index of left-hand residue
        :type at_resid: int
        :param from_idx: global index of right-hand atom in new
        :type from_idx: int
        :param from_resid: global index of right-hand residue
        :type from_resid: int
        :param connected_resids: list of all other resids attached to right-hand residue that should move with it, defaults to []
        :type connected_resids: list, optional
        :return: 2-tuple containing global indices 
        :rtype: tuple
        """
        # Rotate and translate
        if at_resid==from_resid:
            return # should never happen but JIC
        logger.debug(f'{self.name} connected resids {connected_resids}')
        TC=self.TopoCoord
        ATC=TopoCoord()
        BTC=TopoCoord()
        C=TC.gro_DataFrame('atoms')
        ATC.Coordinates.A=C[C['resNum']==at_resid].copy()
        bresids=connected_resids.copy()
        bresids.append(from_resid)
        BTC.Coordinates.A=C[C['resNum'].isin(bresids)].copy()
        for ln in BTC.Coordinates.A.to_string().split('\n'): logger.debug(ln)
        NONROT=C[~C['resNum'].isin(bresids)].shape[0]
        logger.debug(f'{self.TopoCoord.Coordinates.A.shape[0]} atoms')
        logger.debug(f'holding {at_resid} ({NONROT})')
        logger.debug(f'rotating/translating {bresids} ({BTC.Coordinates.A.shape[0]})')
        assert self.TopoCoord.Coordinates.A.shape[0]==(NONROT+BTC.Coordinates.A.shape[0])
        mypartners=TC.partners_of(at_idx)
        otpartners=TC.partners_of(from_idx)
        logger.debug(f'Partners of {at_idx} {mypartners}')
        logger.debug(f'Partners of {from_idx} {otpartners}')
        myHpartners={k:v for k,v in zip(mypartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in mypartners]) if v.startswith('H')}
        otHpartners={k:v for k,v in zip(otpartners,[C[C['globalIdx']==i]['atomName'].values[0] for i in otpartners]) if v.startswith('H')}
        myHighestH={k:v for k,v in myHpartners.items() if v==max([k for k in myHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
        otHighestH={k:v for k,v in otHpartners.items() if v==max([k for k in otHpartners.values()],key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0'))}
        assert len(myHighestH)==1
        assert len(otHighestH)==1
        logger.debug(f'Highest-named H partner of {at_idx} is {myHighestH}')
        logger.debug(f'Highest-named H partner of {from_idx} is {otHighestH}')
        assert len(myHpartners)>0,f'Error: atom {at_idx} does not have a deletable H atom!'
        assert len(otHpartners)>0,f'Error: atom {from_idx} does not have a deletable H atom!'

        Ri=TC.get_R(at_idx)
        Rj=TC.get_R(from_idx)
        logger.debug(f'Ri {at_idx} {Ri} {type(Ri)} {Ri.dtype}')
        logger.debug(f'Rj {from_idx} {Rj} {type(Rj)} {Rj.dtype}')
        overall_maximum=(-1.e9,-1,-1)
        coord_trials={}
        for myH,myHnm in myHpartners.items():  # keys are globalIdx's, values are names
            coord_trials[myH]={}
            Rh=TC.get_R(myH)
            logger.debug(f'  Rh {myH} {Rh} {Rh.dtype}')
            Rih=Ri-Rh
            Rih*=1.0/np.linalg.norm(Rih)
            for otH,otHnm in otHpartners.items():
                logger.debug(f'{self.name}: Considering {myH} {otH}')
                coord_trials[myH][otH]=deepcopy(BTC)
                # logger.debug(f'\n{coord_trials[myH][otH].Coordinates.A.to_string()}')
                Rk=coord_trials[myH][otH].get_R(otH)
                logger.debug(f'{self.name}:    otH {otH} Rk {Rk} {Rk.dtype}')
                Rkj=Rk-Rj
                Rkj*=1.0/np.linalg.norm(Rkj)
                logger.debug(f'Rkj {Rkj} {Rkj.dtype} Rih {Rih} {Rih.dtype}')
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
                logger.debug(f'{self.name}: minD {minD}')
                if minD>overall_maximum[0]:
                    overall_maximum=(minD,myH,otH)
        logger.debug(f'{self.name}: overall_maximum {overall_maximum}')
        minD,myH,otH=overall_maximum
        BTC=coord_trials[myH][otH]
        TC.overwrite_coords(BTC)
        TC.swap_atom_names(myH,list(myHighestH.keys())[0])
        TC.swap_atom_names(otH,list(otHighestH.keys())[0])
        return myH,otH

    def atoms_w_same_attribute_as(self,find_dict={},same_attribute='',return_attribute=''):
        """atoms_w_same_attribute_as returns a list of atom attributes named in the return_attribute parameter from atoms that share an attribute named in the same_attribute parameter with the atom identified by the find_dict parameter

        :param find_dict: dictionary of attribute:value pairs that should uniquely identify an atom, defaults to {}
        :type find_dict: dict, optional
        :param same_attribute: name of attribute used to screen atoms, defaults to ''
        :type same_attribute: str, optional
        :param return_attribute: attribute value to return a list of from the screened atoms, defaults to ''
        :type return_attribute: str, optional
        :return: list of attribute values
        :rtype: list
        """
        att_val=self.TopoCoord.get_gro_attribute_by_attributes(same_attribute,find_dict)
        return self.TopoCoord.get_gro_attributelist_by_attributes(return_attribute,{same_attribute:att_val})

    def flip_stereocenters(self,idxlist):
        """flip_stereocenters flips stereochemistry of atoms in idxlist

        :param idxlist: global indices of chiral atoms
        :type idxlist: list
        """
        self.TopoCoord.flip_stereocenters(idxlist)

    def rotate_bond(self,a,b,deg):
        """rotate_bond rotates all atoms in molecule on b-side of a-b bond by deg degrees

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

    # def sea_of(self,idx):
    #     clu=self.atoms_w_same_attribute_as(find_dict={'globalIdx':idx},
    #                                             same_attribute='sea_idx',
    #                                             return_attribute='globalIdx')
    #     return list(clu)

    def generate_stereoisomers(self):
        """generate_stereoisomers generates list of Molecule shells, one for each stereoisomer

        :return: only returns if no stereoisomers need to be generated
        :rtype: None
        """
        if self.TopoCoord.Topology.D['atoms'].shape[0]==0: return  # self has not yet acquired topology/coordinates
        if len(self.stereoisomers)==0: return
        flip=[[0,1] for _ in range(len(self.stereocenters))]
        st_idx=[self.TopoCoord.get_gro_attribute_by_attributes('globalIdx',{'atomName':n}) for n in self.stereocenters]
        P=product(*flip)
        next(P) # one with no flips is the original molecule, so skip it
        for p in P:
            si_name=self.name+'-S'+''.join([str(_) for _ in p])
            if not si_name in self.stereoisomers:
                logger.debug(f'{si_name} not found in dict of stereoisomers of {self.name}')
            logger.debug(f'Stereocenter sequence {p} generates stereoisomer {si_name}')
            M=self.stereoisomers[si_name]
            M.origin=self.origin
            M.TopoCoord=deepcopy(self.TopoCoord)
            fsc=[st_idx[i] for i in range(len(self.stereocenters)) if p[i]]
            M.flip_stereocenters(fsc)
            M.TopoCoord.write_gro(f'{si_name}.gro')

    def generate_conformers(self):
        """generate_conformers generates this molecule's conformer instances using either gromacs or obabel
        """
        # only generates gro files
        default_gromacs_params={'ensemble': 'nvt', 'temperature': 600, 'ps': 100, 'begin_at': 50}
        # if self.nconformers==0: return
        cd=self.conformers_dict
        if not cd: return
        logger.debug(f'{self.name} conformer_dict {cd}')
        self.nconformers=cd['count']
        minimize=cd.get('minimize',False)
        generator=cd.get('generator',{})
        if generator=='obabel' and not minimize:
            logger.debug(f'Confomers generated by obabel should be energy-minimized.  Indicate "mimimize: True" in the confomers directive for {self.name}')
            minimize=True
        if not generator: return
        logger.info(f'Generating {self.nconformers*(1+len(self.stereoisomers))} conformers for {self.name}')
        gronames=[f'{self.name}']
        nd=cd.get('nzeros',2)
        for si in self.stereoisomers:
            gronames.append(f'{si}')
        for gro in gronames:
            pfx=f'{gro}-C'
            if generator['name']=='obabel':
                compfile=f'{gro}-obabel-confs.gro'
                c=Command(f'obabel -igro {gro}.gro -O {compfile} --conformer --nconf {self.nconformers} --writeconformers')
                out,err=c.run()
                c=Command(f'wc -l {compfile}')
                out,err=c.run()
                tok=out.split()
                lpf=int(tok[0])//self.nconformers
                c=Command(f'split -d -n {nd} -l {lpf} {compfile} {pfx} --additional-suffix=".gro"')
                out,err=c.run()
            elif generator['name']=='gromacs':
                params=generator.get('params',default_gromacs_params)
                begin_at=params.get('begin_at',0.0)
                compfile=f'{gro}-gromacs-confs'
                TC=self.TopoCoord if gro==self.name else self.stereoisomers[gro].TopoCoord
                TC.vacuum_simulate(outname=f'{compfile}',nsamples=cd['count'],params=params)
                gro_from_trr(compfile,nzero=nd,outpfx=pfx,b=begin_at)
            # os.remove(f'{gro}-confs.gro')
            fmt=r'{A}{B:0'+str(nd)+r'd}'  # the trjconv command in gro_from_trr must generate these files
            cfnl=[fmt.format(A=pfx,B=x) for x in range(self.nconformers)]
            for mname in cfnl:
                assert os.path.exists(f'{mname}.gro'),f'Error: Conformer coordinates file {mname}.gro not found'
            logger.debug(f'Conformer coordinate filenames {cfnl}')
            self.conformers.extend(cfnl)
        if minimize:
            saveTC=deepcopy(self.TopoCoord)
            for mname in self.conformers:
                self.TopoCoord.read_gro(f'{mname}.gro')
                logger.info(f'Minimizing conformer {mname}')
                self.TopoCoord.vacuum_minimize(outname=f'{mname}')
            self.TopoCoord=saveTC

MoleculeDict = dict[str,Molecule]
MoleculeList = list[Molecule]

def generate_stereo_reactions(RL:ReactionList,MD:MoleculeDict):
    """generate_stereo_reactions scans the list of reactions and generates any additional reactions in which all possible stereoisomers are reactants

    :param RL: list of Reactions
    :type RL: ReactionList
    :param MD: dictionary of available Molecules
    :type MD: MoleculeDict
    :return: number of new reactions/molecular products created
    :rtype: int
    """
    # any reaction with one or more reactant with one or more stereoisomers 
    # generates new "build" reactions using the stereoisomer as a reactant
    # in place
    adds=0
    terminal_reactions=[]
    for R in RL:
        if R.stage not in [reaction_stage.param,reaction_stage.build]: continue
        Prod=MD[R.product]
        logger.debug(f'Stereos for {R.name} ({str(R.stage)})')
        reactant_stereoisomers={k:[r]+list(MD[r].stereoisomers.keys()) for k,r in R.reactants.items()}
        for k,v in reactant_stereoisomers.items():
            assert all([m in MD for m in v])  # all stereoisomers must be in the dict of molecules
        logger.debug(reactant_stereoisomers)
        reactant_keys=list(R.reactants.keys())
        isomer_lists=list(reactant_stereoisomers.values())
        combos=product(*isomer_lists)
        next(combos)
        # new_reactions=[]
        sidx=1
        for c in combos:
            nR=deepcopy(R)
            nR.name=R.name+f'S-{sidx}'
            nR.product=R.product+f'S-{sidx}'
            nR.stage=reaction_stage.build
            nR.reactants={k:v for k,v in zip(reactant_keys,c)}
            # MD[R.product].stereoisomers[nR.product]=Molecule.NewCopy(MD[R.product],nR.product)
            # add resulting product to global molecule dict so that it will be generated
            MD[nR.product]=Molecule.New(nR.product,nR)
            adds+=1
            MD[nR.product].sequence=MD[R.product].sequence
            MD[nR.product].parentname=R.product
            Prod.stereoisomers[nR.product]=MD[nR.product]
            logger.debug(c)                
            # new_reactions.append(nR)
            sidx+=1
            terminal_reactions.append(nR)
    RL.extend(terminal_reactions)
    return adds

def generate_symmetry_reactions(RL:ReactionList,MD:MoleculeDict):
    """generate_symmetry_reactions scans reaction list to generate any new reactions implied by symmetry-equivalent atoms

    :param RL: list of Reactions
    :type RL: ReactionList
    :param MD: dict of available molecules
    :type MD: MoleculeDict
    :return: number of new reactions/molecular products created
    :rtype: int
    """
    jdx=1
    terminal_reactions=[]
    tail_adds=0
    for R in RL:
        if R.stage not in [reaction_stage.param,reaction_stage.cure,reaction_stage.cap]: continue
        Prod=MD[R.product]
        logger.debug(f'Symmetry versions for {R.name} ({str(R.stage)})\n{str(R)}')
        # thisR_extra_reactions=[]
        # thisR_extra_molecules={}
        # logger.debug(f'  Product {R.product} resname sequence {prod_seq_resn}')
        sra_by_reactant={k:MD[rname].symmetry_relateds for k,rname in R.reactants.items()}
        logger.debug(f'  sra_by_reactant: {sra_by_reactant}')
        atom_options=[]
        for atom_key,atom_rec in R.atoms.items():
            this_atom_options=[]
            art=atom_rec['reactant']
            target_atom_name=atom_rec['atom']
            if art in sra_by_reactant:
                logger.debug(f'art {art} {sra_by_reactant[art]}')
                if len(sra_by_reactant[art])==0: sra_by_reactant[art]=[[target_atom_name]]
                for symm_set in sra_by_reactant[art]:
                    if target_atom_name in symm_set:
                        for atom_name in symm_set:
                            this_atom_options.append([atom_key,atom_name])
            atom_options.append(this_atom_options)
        logger.debug(f'  atom options: {atom_options}')
        if len(R.reactants)>1:
            olist=list(product(*atom_options))
        else:
            olist=list(zip(*atom_options))
        idx=1
        R.symmetry_versions=olist
        logger.debug(f'olist {olist}')
        if len(olist)==1: continue
        for P in olist[1:]:
            newR=deepcopy(R)
            newR.name=R.name+f'-S{idx}'
            logger.debug(f'Permutation {P}:')
            for pp in P:
                atomKey,atomName=pp
                newR.atoms[atomKey]['atom']=atomName
            pname=generate_product_name(newR)
            if len(pname)==0:
                pname=R.product+f'-{idx}'
            newR.product=pname 
            newR.stage=R.stage
            logger.debug(f'Primary:')
            for ln in str(newR).split('\n'): logger.debug(ln)
            terminal_reactions.append(newR)
            MD[newR.product]=Molecule(name=newR.product,generator=newR)
            MD[newR.product].set_origin('symmetry_product')
            MD[newR.product].set_sequence_from_moldict(MD)
            for rR in [x for x in RL if R.product in x.reactants.values()]:
                reactantKey=list(rR.reactants.keys())[list(rR.reactants.values()).index(R.product)]
                logger.debug(f'  New product {newR.product} must replace reactant {reactantKey}:{R.product} in {rR.name}')
                nooR=deepcopy(rR)
                nooR.stage=rR.stage
                nooR.name=rR.name+f'-{reactantKey}:S{jdx}'
                nooR.reactants[reactantKey]=newR.product
                # update any atom names to reflect origin of this reactant
                for naK,naRec in {k:v for k,v in nooR.atoms.items() if v['reactant']==reactantKey}.items():
                    na_resid=naRec['resid'] # resid of reactant atom in target reactant
                    na_name=naRec['atom']
                    for p in P:
                        oaK,oa_name=p 
                        oaRec=R.atoms[oaK]
                        oa_reactatnName=R.reactants[oaRec['reactant']]
                        oa_resid=oaRec['resid']
                        oa_resid_in_o_product=reactant_resid_to_presid(R,oa_reactatnName,oa_resid,RL)
                        # this atom is an atom in the permutation the resid in product matches
                        if na_resid == oa_resid_in_o_product:
                            nooR.atoms[naK]['resid']=oa_resid_in_o_product
                            nooR.atoms[naK]['atom']=oa_name
                noor_pname=generate_product_name(nooR)
                if noor_pname in MD: continue
                if len(noor_pname)==0:
                    noor_pname=rR.product+f'-{jdx}'
                nooR.product=noor_pname
                logger.debug(f'Secondary:')
                for ln in str(nooR).split('\n'): logger.debug(ln)
                jdx+=1
                RL.append(nooR)
                tail_adds+=1
                MD[nooR.product]=Molecule.New(nooR.product,nooR)
                MD[nooR.product].set_origin('symmetry_product')
                MD[nooR.product].set_sequence_from_moldict(MD)
            idx+=1
        logger.debug(f'Symmetry expansion of reaction {R.name} ends')

    RL.extend(terminal_reactions)
    return len(terminal_reactions)+tail_adds