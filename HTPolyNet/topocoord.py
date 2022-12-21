"""

.. module:: topocoords
   :synopsis: Class for jointly handling Topology and Coordinate objects

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from itertools import product
import pandas as pd
import logging
import numpy as np
from enum import Enum
import os
import shutil
from copy import deepcopy
from HTPolyNet.coordinates import Coordinates, GRX_ATTRIBUTES, GRX_GLOBALLY_UNIQUE, GRX_UNSET_DEFAULTS
from HTPolyNet.topology import Topology
from HTPolyNet.bondtemplate import BondTemplate,ReactionBond
from HTPolyNet.gromacs import grompp_and_mdrun,mdp_get, mdp_modify, gmx_energy_trace
import HTPolyNet.projectfilesystem as pfs

logger=logging.getLogger(__name__)

class BTRC(Enum):
    """Bond test return codes: bond tests are applied to those bond-candidates that are within search radius of each other

    :param Enum: inherits from Enum class
    :type Enum: class
    """
    passed = 0
    failed_pierce_ring = 1        # does this bond-candidate pierce a ring?
    failed_short_circuit = 2      # does this bond-candidate result in short-circuit?
    failed_bond_cycle = 3         # does this bond-candidate create a bondcycle on its own?
    unset = 99

class TopoCoord:
    """Container for Topology and Coordinates, along with methods that
        use either or both of them
    """
    def __init__(self,topfilename='',grofilename='',grxfilename='',mol2filename='',system_name='htpolynet',**kwargs):
        """Constructor method for TopoCoord.

        :param topfilename: name of Gromacs-format topology file (top), defaults to ''
        :type topfilename: str, optional
        :param grofilename: name of Gromacs-format coordinate file (gro), defaults to ''
        :type grofilename: str, optional
        :param mol2filename: name of SYBYL MOL2-format coordinate/bonds file, defaults to ''
        :type mol2filename: str, optional
        """
        wrap_coords=kwargs.get('wrap_coords',False)
        # self.basefilenames={}
        self.files={}
        self.files['gro']=os.path.abspath(grofilename)
        self.files['top']=os.path.abspath(topfilename)
        self.files['grx']=os.path.abspath(grxfilename)
        self.files['mol2']=os.path.abspath(mol2filename)
        # self.basefilenames['gro']=grofilename
        # self.basefilenames['top']=topfilename
        # self.basefilenames['grx']=grxfilename
        # self.basefilenames['mol2']=mol2filename
        self.grxattr=[]
        self.idx_lists={}
        self.idx_lists['chain']=[]
        self.idx_lists['cycle']=[]
        if grofilename!='':
            self.read_gro(grofilename,wrap_coords=wrap_coords)
        else:
            self.Coordinates=Coordinates()  # empty
        if topfilename!='':
            self.read_top(topfilename)
        else:
            self.Topology=Topology(system_name=system_name) # empty
        if mol2filename!='':
            self.read_mol2(mol2filename,**kwargs)
        if grxfilename!='':
            self.read_gro_attributes(grxfilename)

    @classmethod
    def from_top_gro(cls,top,gro):
        X=cls(topfilename=top,grofilename=gro)
        return X

    def make_bonds(self,pairs,explicit_sacH={}):
        """Adds new bonds to the global topology

        :param pairs: list of pairs of atom global indices indicating each new bond
        :type pairs: list
        :param skip_H: list of pairs of atom global indices to skip when identifying
            the sacrificial H atoms; likely these are identified during molecule-building to
            optimize mutual orientation and placement of the two reactant molecules, defaults
            to []
        :type skip_H: list, optional
        :return: list of indexes of atoms that must now be deleted (sacrifical H's)
        :rtype: list
        """
        idx_to_ignore=self.Coordinates.find_sacrificial_H(pairs,self.Topology,explicit_sacH=explicit_sacH)
        logger.debug(f'idx_to_ignore {idx_to_ignore}')
        self.Topology.add_bonds(pairs)
        self.chainlist_update(pairs,msg='TopoCoord.make_bonds')
        self.Topology.null_check(msg='add_bonds')
        rename=False if len(explicit_sacH)>0 else False
        idx_to_delete=self.Coordinates.find_sacrificial_H(pairs,self.Topology,explicit_sacH=explicit_sacH,rename=rename)
        assert type(idx_to_delete)==list
        return idx_to_delete

    def add_restraints(self,pairdf,typ=6):
        """Adds bonds of type typ to the topology from the dataframe pairdf

        :param pairdf: ai, aj, initital-distance
        :type pairdf: pandas.DataFrame
        :param typ: bond type to add, defaults to 6 (non-topological restraint)
        :type typ: int, optional
        """
        self.Topology.add_restraints(pairdf,typ=typ)

    def delete_atoms(self,atomlist):
        """Deletes atoms from both the Topology and Coordinates instances

        :param atomlist: list of global indexes of atoms to delete
        :type atomlist: list
        :return: old-to-new atom global index mapper dictionary resulting from reindexing
            remaining atoms to make sure global indexes are sequential
        :rtype: dict
        """
        # logger.debug(f'delete_atoms: {atomlist}')
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper=self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper)==dict
        # logger.debug(f'idx_mapper: {idx_mapper}')
        for list_name in ['chain','cycle']:
            # logger.debug(f'remapping idxs in stale {list_name} lists: {self.idx_lists[list_name]}')
            self.reset_idx_list_from_grx_attributes(list_name)
            # self.remap_idx_list(list_name,idx_mapper)
        # logger.debug(f'finished')
        return idx_mapper

    def count_H(self,idx):
        """count_H Counts the number of hydrogens bound to atom with index idx; 
        any atom whose name begins with the character 'H' or 'h' is assumed to be a hydrogen!

        :param idx: atom global index
        :type idx: int
        :return: number of H atoms bound to this atom
        :rtype: int
        """
        aneigh=self.partners_of(idx)
        aneighnames=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':y}) for y in aneigh]
        anH=sum([int(x.upper().startswith('H')) for x in aneighnames])
        return anH

    def map_from_templates(self,bdf,moldict,overcharge_threshhold=0.1):
        """Updates angles, pairs, dihedrals, atom types, and charges, based on product
            templates associated with each bond in 'bdf'

        :param bdf: dataframe, columns 'ai', 'aj', 'reactantName'
        :type bdf: pandas.DataFrame
        :param moldict: dictionary of template Molecules keyed by name
        :type moldict: dict
        :raises Exception: nan found in any attribute of any new system angle
        :raises Exception: nan found in any attribute of any new system dihedral
        :raises Exception: nan found in any attribute of any new system pair
            that came along with a dihedral
        :raises Exception: nan found in ai attribute of any template pair
        :raises Exception: nan found in aj attribute of any template pair
        :raises Exception: nan found in any system pair
        """
        atdf=self.Topology.D['atoms']
        # ij=self.Topology.D['bondtypes'].set_index(['i','j'])
        #mb=self.D['mol2_bonds']
        # bmi=self.Topology.D['bonds'].set_index(['ai','aj']).sort_index().index
        grodf=self.Coordinates.A
        grodf['old_reactantName']=grodf['reactantName'].copy()
        logger.debug(f'Mapping {bdf.shape[0]} bonds.')
        premapping_total_charge=self.Topology.total_charge()
        logger.debug(f'Must compensate for an overcharge of {premapping_total_charge:.4f}')
        mapped_inst_atoms=[]
        for b in bdf.itertuples():
            logger.debug(f'Mapping bond {b}')
            # for ln in b.to_string().split('\n'):
            #     logger.debug(f'  -> {ln}')
            bb=[b.ai,b.aj]
            order=b.order
            names=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in bb]
            resnames=[self.get_gro_attribute_by_attributes('resName',{'globalIdx':x}) for x in bb]
            logger.debug(f'{resnames}')
            resids=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in bb]
            # this is the product name of the reaction used to identify this bond
            product_name=b.reactantName
            if product_name in moldict:
                P=moldict[product_name]
                if P.is_reactant:
                    logger.debug(f'{P.name} is a reactant; updating \"reactantName\" attributes of {bb}')
                    # if the product of this reaction is also a reactant in a cure reaction, we should
                    # change the reactantName attribute of the two atoms to match this so either of
                    # these atoms can be found in a later bond search
                    for j in bb:
                        self.set_gro_attribute_by_attributes('reactantName',product_name,{'globalIdx':j})
                else:
                    logger.debug(f'{P.name} is not a reactant; no update of \"reactantName\" attributes')
            bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames=self.get_bystanders(bb)
            oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames=self.get_oneaways(bb)
            intraresidue=resids[0]==resids[1]
            BT=BondTemplate(names,resnames,intraresidue,order,bystander_resnames,bystander_atomnames,oneaway_resnames,oneaway_atomnames)
            RB=ReactionBond(bb,resids,order,bystander_resids,bystander_atomidx,oneaway_resids,oneaway_atomidx)
            logger.debug(f'apparent bond template {str(BT)}')
            logger.debug(f'apparent bond instance {str(RB)}')
            T,rb,reverse_bond=find_template(BT,moldict)
            if reverse_bond: RB.reverse()
            temp_i_idx,temp_j_idx=rb.idx
            d=T.TopoCoord.Topology.D['bonds']
            # copy all bond records matching these two bonds; should be only one!!
            d=d[((d.ai==temp_i_idx)&(d.aj==temp_j_idx))|
                ((d.ai==temp_j_idx)&(d.aj==temp_i_idx))].copy()
            assert d.shape[0]==1,f'Using {T.name} is sent inst-bond {i_idx}-{j_idx} which is claimed to map to {temp_i_idx}-{temp_j_idx}, but no such unique bond is found:\n{T.TopoCoord.Topology.D["bonds"].to_string()}'
            ''' check passed '''
            temp_angles,temp_dihedrals,temp_pairs=T.get_angles_dihedrals((temp_i_idx,temp_j_idx))
            logger.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            # determine the set of unique atoms in template that must be mapped to instance -- it is all
            # atoms involved in these bonded interactions
            uniq_temp_idx=set()
            for df in [temp_angles,temp_dihedrals,temp_pairs]:
                for label in ['ai','aj','ak','al']:
                    if label in df:
                        uniq_temp_idx=uniq_temp_idx.union(set(df[label].to_list()))
            logger.debug(f'Template atom indexes that must be mapped: {uniq_temp_idx}')
            # get the bidirectional instance<->template mapping dictionaries
            inst2temp,temp2inst=T.idx_mappers(self,RB.idx,RB.bystander_resids,RB.oneaway_resids,uniq_temp_idx)
            # some hard checks on compatibility of the dicts
            assert len(inst2temp)==len(temp2inst)
            check=True
            for k,v in inst2temp.items():
                check = check and (k == temp2inst[v])
            for k,v in temp2inst.items():
                check = check and (k == inst2temp[v])
            assert check,f'Error: bidirectional dicts are incompatible; bug\n{inst2temp}\b{temp2inst}'
            # logger.debug(f'inst2temp {inst2temp}')
            # logger.debug(f'temp2inst {temp2inst}')
            i_idx,j_idx=bb
            _temp_i_idx,_temp_j_idx=inst2temp[i_idx],inst2temp[j_idx]
            assert temp_i_idx==_temp_i_idx,f'mapping mismatch -- bug'
            assert temp_j_idx==_temp_j_idx,f'mapping mismatch -- bug'

            need_new_bond_parameters=False
            total_dcharge=0.0
            temp_atdf=T.TopoCoord.Topology.D['atoms']
            for temp_atom,inst_atom in temp2inst.items():
                assert inst_atom in atdf['nr'].values,f'Error: mapped atom {inst_atom} not found in [ atoms ]'
                inst_type,inst_charge,inst_name,inst_resn,inst_rnam=atdf[atdf['nr']==inst_atom][['type','charge','atom','resnr','residue']].values[0]
                temp_type,temp_charge,temp_name,temp_resn,temp_rnam=temp_atdf[temp_atdf['nr']==temp_atom][['type','charge','atom','resnr','residue']].values[0]
                # logger.debug(f'temp {temp_atom} {temp_name} {temp_rnam} {temp_resn} {temp_type} {temp_charge}')
                # logger.debug(f'inst {inst_atom} {inst_name} {inst_rnam} {inst_resn} {inst_type} {inst_charge}')
                if inst_type!=temp_type:
                    logger.debug(f'changing type of inst atom {inst_atom} ({inst_resn} {inst_rnam} {inst_name}) from {inst_type} to {temp_type}')
                    atdf.loc[atdf['nr']==inst_atom,'type']=temp_type
                    if inst_atom==b.ai or inst_atom==b.aj:
                        # changed type of one or both of the bond atoms
                        need_new_bond_parameters=True
                if inst_charge!=temp_charge:
                    logger.debug(f'charge {inst_atom} ({inst_resn} {inst_rnam} {inst_name}) from {inst_charge} to {temp_charge}')
                    atdf.loc[atdf['nr']==inst_atom,'charge']=temp_charge
                    dcharge=temp_charge-inst_charge
                    total_dcharge+=dcharge
            mapped_inst_atoms.extend(list(temp2inst.values()))
            if need_new_bond_parameters:
                self.Topology.reset_override_from_type('bonds','bondtypes',inst_idx=(b.ai,b.aj))

            # get all angle, dihedrals, and pairs from template that result from the existence of the specified bond
            # temp_angles,temp_dihedrals,temp_pairs=T.get_angles_dihedrals((temp_i_idx,temp_j_idx))
            # logger.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            # map from template atom indicies to system atom indicies in angles
            # logger.debug(f'Template angles:')
            # for ln in temp_angles.to_string().split('\n'):
            #     logger.debug(ln)
            inst_angles=temp_angles.copy()
            inst_angles.ai=temp_angles.ai.map(temp2inst)
            inst_angles.aj=temp_angles.aj.map(temp2inst)
            inst_angles.ak=temp_angles.ak.map(temp2inst)
            # logger.debug(f'Mapped instance angles:')
            # for ln in inst_angles.to_string().split('\n'):
            #     logger.debug(ln)
            # add new angles to the system topology
            d=self.Topology.D['angles']
            self.Topology.D['angles']=pd.concat((d,inst_angles),ignore_index=True)                            
            # hard check for any nan's in any atom index attribute in any angle
            d=self.Topology.D['angles']
            check=True
            for a in ['ai','aj','ak']:
                check=check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in angles')
                raise Exception

            # logger.debug(f'Template dihedrals:')
            # for ln in temp_dihedrals.to_string().split('\n'):
            #     logger.debug(ln)
            # map from template atom indicies to system atom indicies in dihedrals
            inst_dihedrals=temp_dihedrals.copy()
            inst_dihedrals.ai=temp_dihedrals.ai.map(temp2inst)
            inst_dihedrals.aj=temp_dihedrals.aj.map(temp2inst)
            inst_dihedrals.ak=temp_dihedrals.ak.map(temp2inst)
            inst_dihedrals.al=temp_dihedrals.al.map(temp2inst)
            # logger.debug(f'Mapped instance dihedrals:')
            # for ln in inst_dihedrals.to_string().split('\n'):
            #     logger.debug(ln)
            d=inst_dihedrals
            check=False
            for a in ['ai','aj','ak','al']:
                check=check or d[a].isnull().values.any()
            if check:
                logger.error(f'a {a} NAN in dihedrals\n{inst_dihedrals.to_string()}')
                raise Exception
            # add new dihedrals to global topology
            d=self.Topology.D['dihedrals']
            self.Topology.D['dihedrals']=pd.concat((d,inst_dihedrals),ignore_index=True)
            # hard check for no nans
            d=self.Topology.D['dihedrals']
            check=True
            for a in ['ai','aj','ak','al']:
                check=check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in dihedrals')
                raise Exception
            d=self.Topology.D['pairs']
            check=True
            for a in ['ai','aj']:
                check=check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in pairs premapping')
                raise Exception

            # double-hard check to make sure pairs can be mapped
            k=np.array(list(temp2inst.keys()))
            v=np.array(list(temp2inst.values()))
            if any(np.isnan(k)):
                logger.error('null in temp2inst keys')
            if any(np.isnan(v)):
                logger.error('null in temp2inst values')
            # logger.debug(f'temp_pairs:\n{temp_pairs.to_string()}')
            isin=[not x in temp2inst for x in temp_pairs.ai]
            if any(isin):
                for ii,jj in enumerate(isin):
                    if jj:
                        logger.error(f'atom ai {temp_pairs.ai.iloc[ii]} not in temp2inst')
            isin=[not x in temp2inst for x in temp_pairs.aj]
            if any(isin):
                for ii,jj in enumerate(isin):
                    if jj:
                        logger.error(f'atom aj {temp_pairs.aj.iloc[ii]} not in temp2inst')

            # map all ai attributes of all template pairs to global ai
            temp_pairs.ai=temp_pairs.ai.map(temp2inst)
            # check AGAIN for nans (I am afraid of nans)
            if temp_pairs.ai.isnull().values.any():
                logger.error('NAN in pairs ai')
                raise Exception
            # map all aj attributes of all template pairs to global ai
            temp_pairs.aj=temp_pairs.aj.map(temp2inst)
            # check YET AGAIN for nans (eek!)
            if temp_pairs.aj.isnull().values.any():
                logger.error('NAN in pairs aj')
                raise Exception
            # add these pairs to the topology
            # logger.debug(f'Concatenating this pairs to global pairs')
            # for ln in temp_pairs.to_string().split('\n'):
            #     logger.debug(ln)
            self.Topology.D['pairs']=pd.concat((d,temp_pairs),ignore_index=True)
            d=self.Topology.D['pairs']
            # check AGAIN for nans
            check=True
            for a in ['ai','aj']:
                check=check and d[a].isnull().values.any()
            if check:
                logger.error('NAN in pairs post mapping')
                raise Exception
            # return temp_pairs
        mapped_inst_atoms=list(set(mapped_inst_atoms))
        logger.debug(f'System overcharge after mapping: {self.Topology.total_charge():.4f}')
        self.adjust_charges(atoms=mapped_inst_atoms,overcharge_threshhold=overcharge_threshhold,msg=f'overcharge magnitude exceeds {overcharge_threshhold}')

    def enumerate_1_4_pairs(self,at_idx):
        """enumerate_1_4_pairs enumerate all 1-4 pair interactions resulting from new bonds in at_idx

        :param at_idx: list of 2-tuples, each containing global indices of atoms that bond to each other
        :type at_idx: list
        :return: dataframe of all pairs
        :rtype: pandas.DataFrame
        """
        pdf=self.Topology.D['pairs']
        # each of these bonds results in 1-4 pair interactions
        bl=self.Topology.bondlist
        pai=[]
        paj=[]
        pri=[]
        prj=[]
        prsource=[]
        for p in at_idx:
            # determine all nj--nk pairs from j-k bond, where nj are neighbors of j not including k, and nk are neighbors of k not including j
            # and all nnj--k pairs from all nj-j bonds, where nnj are next-nearest neighbors of j excluding j
            # and all j--nnk pairs from all k-nk bonds, where nnk are next-nearest neighbors of k excluding k
            j,k=p
            nj=bl.partners_of(j)
            nj.remove(k)
            nnj=[]
            for nn in nj:
                t=bl.partners_of(nn)
                t.remove(j)
                nnj.append(t)
            nk=bl.partners_of(k)
            nk.remove(j)
            nnk=[]
            for nn in nk:
                t=bl.partners_of(nn)
                t.remove(k)
                nnk.append(t)
            logger.debug(f'{j} {nj} {k} {nk}')
            logger.debug(f'{nnj} {nnk}')
            for pp in product(nj,nk):
                jj,kk=pp
                pai.append(jj)
                paj.append(kk)
                pri.append(j)
                prj.append(k)
                prsource.append('c')
                for nn in nnj:
                    for ii in nn:
                        pai.append(ii)
                        paj.append(k)
                        pri.append(j)
                        prj.append(k)
                        prsource.append('l')
                for nn in nnk:
                    for ll in nn:
                        pai.append(j)
                        paj.append(ll)
                        pri.append(j)
                        prj.append(k)
                        prsource.append('r')

        pi_df=pd.DataFrame({'ai':pai,'aj':paj,'source':prsource,'bi':pri,'bj':prj})
        pi_df.drop_duplicates(inplace=True,ignore_index=True)
        return pi_df


    def update_topology_and_coordinates(self,bdf,template_dict={},write_mapper_to=None,**kwargs):
        """update_topology_and_coordinates updates global topology and necessary atom attributes in the configuration to reflect formation of all bonds listed in `keepbonds`

        :param bdf: bonds dataframe, columns 'ai', 'aj', 'reactantName'
        :type bdf: pandas.DataFrame
        :param template_dict: dictionary of molecule templates keyed on molecule name
        :type template_dict: dict
        :return: 3-tuple: new topology file name, new coordinate file name, list of bonds with atom indices updated to reflect any atom deletions
        :rtype: 3-tuple
        
        """
        explicit_sacH=kwargs.get('explicit_sacH',{})
        template_source=kwargs.get('template_source','internal')
        overcharge_threshhold=kwargs.get('overcharge_threshhold',0.1)
        logger.debug(f'begins.')
        if bdf.shape[0]>0:
            assert bdf['ai'].dtype==int
            assert bdf['aj'].dtype==int
            # pull out just the atom index pairs (first element of each tuple)
            at_idx=[(int(x.ai),int(x.aj),x.order) for x in bdf.itertuples()]
            logger.debug(f'Making {len(at_idx)} bonds.')
            idx_to_delete=self.make_bonds(at_idx,explicit_sacH=explicit_sacH)
            logger.debug(f'Deleting {len(idx_to_delete)} atoms.')
            idx_mapper=self.delete_atoms(idx_to_delete) # will result in full reindexing
            # logger.debug(f'null check')
            self.Topology.null_check(msg='delete_atoms')
            # reindex all atoms in the list of bonds sent in, and write it out
            logger.debug(f'z-decrement, nreactions increment')
            ri_bdf=bdf.copy()
            ri_bdf.ai=ri_bdf.ai.map(idx_mapper)
            ri_bdf.aj=ri_bdf.aj.map(idx_mapper)
            at_idx=[(x.ai,x.aj) for x in ri_bdf.itertuples()]
            for idx_pair in at_idx:
                for idx in idx_pair:
                    self.decrement_gro_attribute_by_attributes('z',{'globalIdx':idx})
                    self.increment_gro_attribute_by_attributes('nreactions',{'globalIdx':idx})
            if template_source=='internal':
                logger.debug(f'calling map_from_templates')
                self.map_from_templates(ri_bdf,template_dict,overcharge_threshhold=overcharge_threshhold)
            logger.debug(f'1-4 pair update')
            pi_df=self.enumerate_1_4_pairs(at_idx)
            self.Topology.null_check(msg='update_topology_and_coordinates')
            if write_mapper_to:
                tdf=pd.DataFrame({'old':list(idx_mapper.keys()),'new':list(idx_mapper.values())})
                tdf.to_csv(write_mapper_to,sep=' ',index=False)
            logger.debug('finished')
            return ri_bdf,pi_df

    def read_top(self,topfilename):
        """read_top Creates a new Topology member by reading from a Gromacs-style top file.
            Just a wrapper for the read_gro method of Topology

        :param topfilename: name of topology file
        :type topfilename: str
        """
        self.files['top']=os.path.abspath(topfilename)
        self.Topology=Topology.read_gro(topfilename)

    def read_gro(self,grofilename,preserve_box=False,wrap_coords=False):
        """read_gro Creates a new Coordinates member by reading from a Gromacs-style coordinates
            file.  Just a wrapper for the read_gro method of Coordinates

        :param grofilename: name of gro file
        :type grofilename: str
        """
        self.files['gro']=os.path.abspath(grofilename)
        if preserve_box:
            savebox=self.Coordinates.box.copy()
        self.Coordinates=Coordinates.read_gro(grofilename,wrap_coords=wrap_coords)
        if preserve_box:
            self.Coordinates.box=savebox
        # logger.debug(f'box: {self.Coordinates.box}')

    def read_mol2(self,mol2filename,**kwargs):
        """read_mol2 Creates a new Coordinates member by reading from a SYBYL-style MOL2 file.
            A wrapper for read_mol2 from Coordinates, but also sets the 'mol2_bonds'
            dataframe in the Topology if the parameter ignore_bonds is False.  If
            the mol2_bonds dataframe is created, and the Topology already has a 'bonds' dataframe, a consistency check is peformed.

        :param mol2filename: name of mol2 file
        :type mol2filename: str
        """
        ignore_bonds=kwargs.get('ignore_bonds',False)
        overwrite_coordinates=kwargs.get('overwrite_coordinates',False)
        self.files['mol2']=os.path.abspath(mol2filename)
        temp_coords=Coordinates.read_mol2(mol2filename)
        if self.Coordinates.empty or overwrite_coordinates:
            save_box=self.Coordinates.box.copy()
            self.Coordinates=temp_coords
            self.Coordinates.box=save_box
            temp_coords=self.Coordinates
        if not ignore_bonds:
            self.Topology.D['mol2_bonds']=temp_coords.mol2_bonds.copy()
            if 'bonds' in self.Topology.D:
                self.Topology.bond_source_check()

    def swap_atom_names(self,ai,aj):
        """swap_atom_names Swaps the names of the two atoms with global indicies ai and aj.  This is used when
        automatically selected which of several possible sacrificial H's will actually be
        selected.  Surviving H's are renamed so it always appears that the H with "least
        important" name (lowest order if sorted) is the sacrificial H.  Why do we do this?  It
        gives us perfect control of the names of the atoms that survive a bond.

        :param ai: global index of first atom
        :type ai: int
        :param aj: global index of second atom
        :type aj: int
        """
        T=self.Topology.D['atoms']
        C=self.Coordinates.A
        l1=T.columns=='atom'
        l2=C.columns=='atomName'
        iname=T.iloc[ai-1,l1].values[0]
        jname=T.iloc[aj-1,l1].values[0]
        # logger.debug(f'Swapping names of atoms {ai}({iname}) and {aj}({jname})')
        tmpNm=T.iloc[ai-1,l1].values[0]
        T.iloc[ai-1,l1]=T.iloc[aj-1,l1]
        T.iloc[aj-1,l1]=tmpNm
        C.iloc[ai-1,l2]=C.iloc[aj-1,l2]
        C.iloc[aj-1,l2]=tmpNm

    def read_top_gro(self,topfilename,grofilename):
        """read_top_gro Wrapper for read_top and read_gro; generates new Topology and Coordinates members

        :param topfilename: name of topology file
        :type topfilename: str
        :param grofilename: name of coordinates file (Gromacs format)
        :type grofilename: str
        """
        self.read_top(topfilename)
        self.read_gro(grofilename)

    def write_top(self,topfilename):
        """write_top Write a Gromacs-format topology file; this will only write an in-line version,
            no itp; wrapper for Topology.to_file()

        :param topfilename: name of file to write
        :type topfilename: str
        """
        self.Topology.to_file(topfilename)
        self.files['top']=os.path.abspath(topfilename)

    def write_gro(self,grofilename,grotitle=''):
        """write_gro Write a Gromacs-format coordinate file; wrapper for Coordinates.write_gro()

        :param grofilename: name of file to write
        :type grofilename: str
        """
        self.Coordinates.write_gro(grofilename,grotitle=grotitle)
        self.files['gro']=os.path.abspath(grofilename)

    def write_top_gro(self,topfilename,grofilename):
        """write_top_gro Writes both a Gromacs top file and Gromacs coordinate file

        :param topfilename: name of topology file to write
        :type topfilename: str
        :param grofilename: name of coordinate file to write
        :type grofilename: str
        """
        self.write_top(topfilename)
        self.write_gro(grofilename)

    def return_bond_lengths(self,bdf):
        """return_bond_lengths Return the length of all bonds in list bonds

        :param bdf: bonds dataframe, 'ai','aj','reactantName'
        :type bonds: pandas.DataFrame
        :return: list of lengths parallel to bonds
        :rtype: list of floats
        """
        return self.Coordinates.return_bond_lengths(bdf)

    def add_length_attribute(self,bdf:pd.DataFrame,attr_name='length'):
        """add_length_attribute computes bond lengths based on bonds indicated by the parallel 'ai' and 'aj' columns of the parameter dataframe bdf and stores result in a new column called attr_name

        :param bdf: a pandas dataframe with 'ai' and 'aj' columns of atom indices indicating bonds
        :type bdf: pd.DataFrame
        :param attr_name: name of length attribute column, defaults to 'length'
        :type attr_name: str, optional
        """
        self.Coordinates.add_length_attribute(bdf,attr_name=attr_name)

    def copy_bond_parameters(self,bonds):
        """Generate and return a copy of a bonds dataframe that contains all bonds
           listed in bonds

        :param bonds: list of bonds, each a 2-tuple of global atom indices
        :type bonds: list
        :return: bonds dataframe
        :rtype: pandas.DataFrame
        """
        return self.Topology.copy_bond_parameters(bonds)

    def remove_restraints(self,pairsdf):
        """Remove all bonds represented in in pairdf.
        These are interpreted as non-topological
        restraints, so deleting these 'bonds' does 
        not influence angles or dihedrals

        :param pairdf: dataframe of pairs ['ai','aj']
        :type pairdf: pandas DataFrame
        """
        self.Topology.remove_restraints(pairsdf)

    def attenuate_bond_parameters(self,bonds,i,n,minimum_distance=0.0,init_colname='initial_distance'):
        """Alter the kb and b0 parameters for new crosslink bonds according to the values prior to
            relaxation (stored in lengths), their equilibrium values, and the ratio stage/max_stages.
            Let stage/max_stages be x, and 1/max_stages <= x <= 1.  The spring constant for each
            bond is multiplied by x and the distance is 1 xth of the way from its maximum value
            to its equilibrium value.

        :param bonds: bonds dataframe, 'ai', 'aj', 'initial_distance'
        :type bonds: pandas.DataFrame
        :param stage: index of stage in the series of post-bond-formation relaxation
        :type stage: int
        :param max_stages: total number of relaxation stages for this iteration
        :type max_stages: int
        :param minimum_distance: minimum bondlegth allowed, overriding type-specific b0
        :type lengths: float
        """
        self.Topology.attenuate_bond_parameters(bonds,i,n,minimum_distance=minimum_distance,init_colname=init_colname)

    def attenuate_pair_parameters(self,pairdf,i,n,draglimit_nm=0.3):
        """Alter the kb and b0 parameters for new pre-crosslink pairs according
            to the values prior to dragging (stored in pairdf['initial_distances']),
            the desired lower limit of interatomic distance 'draglimit_nm',
            and the ratio stage/max_stages.

        :param pairdf: pairs dataframe (['ai'],['aj'],['initial_distance'])
        :type pairdf: pandas.DataFrame
        :param stage: index of stage in the series of pre-bond-formation dragging
        :type stage: int
        :param max_stages: total number of drag stages for this iteration
        :type max_stages: int
        :param draglimit_nm: lower limit of interatomic distance requested from drag
        :type draglimit_nm: float
        """
        self.Topology.attenuate_pair_parameters(pairdf,i,n,minimum_distance=draglimit_nm)

    def copy_coords(self,other):
        """Copy coordinates and box size from other to self

        :param other: a TopoCoord instance
        :type other: TopoCoord
        """
        self.Coordinates.copy_coords(other.Coordinates)
        self.Coordinates.box=other.Coordinates.box.copy()
        self.files['gro']=os.path.abspath(other.files['gro'])

    def restore_bond_parameters(self,saved):
        """Retores saved bond parameters in df saved by overwriting

        :param saved: [ bonds ] dataframe
        :type saved: pandas.DataFrame
        """
        self.Topology.restore_bond_parameters(saved)

    def set_grx_attributes(self,attributes=[]):
        """set_grx_attributes override the global GRX_ATTRIBUTES

        :param attributes: new GRX attributes to use, defaults to []
        :type attributes: list, optional
        """
        if len(attributes)==0:
            self.grxattr=GRX_ATTRIBUTES
        else:
            self.grxattr=attributes
        logger.debug(f'grxattr set to {self.grxattr}')
    
    def write_gro_attributes(self,attributes_list,grxfilename):
        """write_gro_attributes Writes atomic attributes to a file

        :param attributes_list: list of attributes to write
        :type attributes_list: list
        :param grxfilename: name of output file
        :type grxfilename: str
        """
        self.Coordinates.write_atomset_attributes(attributes_list,grxfilename)
        self.files['grx']=os.path.abspath(grxfilename)

    def write_grx_attributes(self,grxfilename):
        """write_grx_attributes Writes GRX attributes to a file

        :param grxfilename: name of output file
        :type grxfilename: str
        """        
        self.write_gro_attributes(self.grxattr,grxfilename)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
        """Read attributes from file into self.Coordinates.A

        :param grxfilename: name of input file
        :type grxfilename: str
        :param attribute_list: list of attributes to take, defaults to [] (take all)
        :type attribute_list: list, optional
        """
        self.files['grx']=os.path.abspath(grxfilename)
        logger.debug(f'Reading {grxfilename}')
        attributes_read=self.Coordinates.read_atomset_attributes(grxfilename,attributes=attribute_list)
        if 'chain' in attributes_read and 'chain_idx' in attributes_read:
            self.reset_idx_list_from_grx_attributes('chain')
        if 'cycle' in attributes_read and 'cycle_idx' in attributes_read:
            self.reset_idx_list_from_grx_attributes('cycle')
        if attributes_read!=self.grxattr:
            self.grxattr=attributes_read

    def set_gro_attribute(self,attribute,srs):
        """set_gro_attribute sets attribute of atoms to srs (drillst through to Coordinates.set_atomset_attributes())

        :param attribute: name of attribute
        :type attribute: str
        :param srs: scalar or list-like attribute values in same ordering as self.A
        :type srs: scalar or list-like
        """
        self.Coordinates.set_atomset_attribute(attribute,srs)

    def set_gro_attribute_by_attributes(self,att_name,att_value,attribute_dict):
        """set_atom_attribute set the attributes named in name to values named in values (names||values) for the set of atoms specified in the attributes dict

        :param name: list of names of attributes
        :type name: list
        :param value: list of values of attributes to be set
        :type value: list
        :param attributes: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attributes: dict
        """
        self.Coordinates.set_atom_attribute(att_name,att_value,attribute_dict)

    def get_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """get_gro_attribute_by_attributes return values of attributes listed in name from atoms specified by attribute:value pairs in attribute_dict

        :param att_name: list of attributes whose values are to be returned
        :type att_name: list
        :param attribute_dict: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attribute_dict: dict
        :return: scalar or list of one or more return attribute values
        :rtype: list if name is a list; scalar otherwise
        """
        return self.Coordinates.get_atom_attribute(att_name,attribute_dict)

    def increment_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """increment_gro_attribute_by_attributes add one to attribute att_name of all atoms identified by attribute:value pairs in attribute_dict

        :param att_name: name of attribute to increment
        :type att_name: str
        :param attribute_dict: attribute:value pairs that specify atoms to which to apply this incrementation
        :type attribute_dict: dict
        """
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        val+=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def decrement_gro_attribute_by_attributes(self,att_name,attribute_dict):
        """decrement_gro_attribute_by_attributes subtract one from attribute att_name of all atoms identified by attribute:value pairs in attribute_dict

        :param att_name: name of attribute to increment
        :type att_name: str
        :param attribute_dict: attribute:value pairs that specify atoms to which to apply this incrementation
        :type attribute_dict: dict
        """
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        val-=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def get_gro_attributelist_by_attributes(self,attribute_list,attribute_dict):
        """get_atoms_w_attribute returns all rows of atoms dataframe and columns named in names of atoms identified by the attributes dict

        :param attribute_list: list of attributes to be in the rows that are returned
        :type sttribute_list: list
        :param attribute_dict: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attribute_dict: dict
        :return: a dataframe segment
        :rtype: pd.DataFrame
        """
        return self.Coordinates.get_atoms_w_attribute(attribute_list,attribute_dict)

    def get_R(self,idx):
        """get_R return the cartesian position of atom with global index idx

        :param idx: atom global index
        :type idx: int
        :return: position of atom
        :rtype: numpy.ndarray(3,float)
        """
        return self.Coordinates.get_R(idx)

    def rotate(self,R):
        """rotate applies rotation matrix R to all atom positions

        :param R: rotation matrix
        :type R: numpy.ndarray((3,3),float)
        """
        self.Coordinates.rotate(R)

    def translate(self,L):
        """translate applies translation vector L to all atom positions

        :param L: translation vector
        :type L: numpy.ndarray(3,float)
        """
        self.Coordinates.translate(L)

    def partners_of(self,i):
        """partners_of return list of atom indices of bonded partners of atom i

        :param i: global index of atom
        :type i: int
        :return: list of partner global indices
        :rtype: list
        """
        return self.Topology.bondlist.partners_of(i)

    def resid_partners_of(self,ri):
        """resid_partners_of return list of resid partners of resid ri

        :param ri: residue index
        :type ri: int
        :return: list of partner residue indices (two residues are partners if there is at least one interatomic bond joining them)
        :rtype: list
        """
        result=[]
        adf=self.Coordinates.A
        radf=adf[adf['resNum']==ri]
        for at in radf['globalIdx'].to_list():
            bl=self.Topology.bondlist.partners_of(at)
            for j in bl:
                theirresid=adf.iloc[j-1]['resNum']
                if theirresid!=ri and not theirresid in result:
                    result.append(theirresid)
        return result

    def interresidue_partners_of(self,i):
        """interresidue_partners_of return list of atom indices that are bonded partners of atom i and are not in the residue of atom i

        :param i: atom global index
        :type i: int
        :return: list of atom indices of atoms that are partners of i not in i's residue
        :rtype: list
        """
        result=[]
        bl=self.Topology.bondlist.partners_of(i)
        # logger.debug(f'{i} partners {bl}')
        myresid=self.Coordinates.A.iloc[i-1]['resNum']
        for j in bl:
            theirresid=self.Coordinates.A.iloc[j-1]['resNum']
            if theirresid!=myresid:
                result.append(j)
        return result

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        """minimum_distance computes the distance of closest approach between self's Coordinates that other's Coordinates

        :param other: another TopoCoord object
        :type other: TopoCoord
        :param self_excludes: list of global atom indices to ignore in self, defaults to []
        :type self_excludes: list, optional
        :param other_excludes: list of global atom indices to ignore in other, defaults to []
        :type other_excludes: list, optional
        :return: distance of closest approach: i.e., the distance between the two atoms, one from self and one from other, that are closest together
        :rtype: float
        """
        return self.Coordinates.minimum_distance(other.Coordinates,self_excludes=self_excludes,other_excludes=other_excludes)

    # def has_gro_attributes(self,attribute_list):
    #     return self.Coordinates.has_atom_attributes(attribute_list)

    def are_bonded(self,i,j):
        """are_bonded checks to see if atoms with indices i and j are bonded to each other

        :param i: an atom global index
        :type i: int
        :param j: another atom global index
        :type j: int
        :return: True if atoms are bonded, False otherwise
        :rtype: bool
        """
        return self.Topology.bondlist.are_bonded(i,j)

    # def decrement_z(self,pairs):
    #     self.Coordinates.decrement_z(pairs)

    def adjust_charges(self,atoms=[],overcharge_threshhold=0.1,netcharge=0.0,msg=''):
        """adjust_charges adjust the partial charges on atoms in list 'atoms' if the absolute net charge exceeds 'netcharge' by the 'overcharge_threshhold' 

        :param atoms: list of atom indexes to consider, defaults to []
        :type atoms: list, optional
        :param overcharge_threshhold: absolute deviation from netcharge that triggers adjustment, defaults to 0.1
        :type overcharge_threshhold: float, optional
        :param netcharge: desired net charge, defaults to 0.0
        :type netcharge: float, optional
        :param msg: a little message to echo to console if adjustment is necessary, defaults to ''
        :type msg: str, optional
        """
        self.Topology.adjust_charges(atoms=atoms,overcharge_threshhold=overcharge_threshhold,desired_charge=netcharge,msg=msg)

    def gro_DataFrame(self,name):
        """gro_DataFrame return the appropriate Coordinates dataframe based on the directive in 'name'

        :param name: either 'atoms' or 'mol2_bonds'
        :type name: str
        :return: the desired dataframe
        :rtype: pandas.DataFrame
        """
        if name=='atoms':
            return self.Coordinates.A
        elif name=='mol2_bonds':
            return self.Coordinates.mol2_bonds
        else:
            return None

    def overwrite_coords(self,other):
        """overwrite_coords overwrite coordinates in self by those in other

        :param other: another TopoCoord object
        :type other: TopoCoord
        """
        # logger.debug(f'Overwriting {other.Coordinates.A.shape[0]} coordinates')
        C=self.Coordinates.A
        C=C.set_index('globalIdx')
        # logger.debug(f'before update:\n{C.to_string()}')
        B=other.Coordinates.A
        B=B.set_index('globalIdx')
        B=B[['posX','posY','posZ']].copy()
        # logger.debug(f'new coordinates:\n{B.to_string()}')
        C.update(B)
        self.Coordinates.A=C.reset_index()
        # logger.debug(f'after update:\n{self.Coordinates.A.to_string()}')

    def linkcell_initialize(self,cutoff,ncpu=1,force_repopulate=True):
        """linkcell_initialize Initialize the linkcell structure; a wrapper for Coordinates

        :param cutoff: minimum value of cell side-length
        :type cutoff: float
        :param ncpu: number of processors to use in populating linkcell structure in parallel, default 1
        :type ncpu: int
        """
        self.Coordinates.linkcell_initialize(cutoff,ncpu=ncpu,populate=True,force_repopulate=force_repopulate)

    def linkcell_cleanup(self):
        """linkcell_cleanup removes linkcell_idx attribute from Coordinate.A
        """
        self.Coordinates.A.drop(columns=['linkcell_idx'],inplace=True)

    def atom_count(self):
        """atom_count Check to be sure the Coordinate and Topology members contain the same number of
            atoms

        :return: the number of atoms
        :rtype: int
        """
        assert self.Coordinates.A.shape[0]==self.Topology.D['atoms'].shape[0]
        return self.Coordinates.A.shape[0]

    def total_mass(self,units='SI'):
        """total_mass Returns the total mass of the system.  Just a wrapper.

        :param units: units designation, defaults to 'SI' (other option is 'gromacs')
        :type units: str, optional
        :return: mass
        :rtype: float
        """
        return self.Topology.total_mass(units=units)

    def total_volume(self,units='SI'):
        """total_volume returns total volume represented by the system's Coordinates

        :param units: unit designation, defaults to 'SI'
        :type units: str, optional
        :return: box volume
        :rtype: float
        """
        return self.Coordinates.box_volume(units=units)

    def density(self,units='SI'):
        """density returns system density

        :param units: unit designation, defaults to 'SI'
        :type units: str, optional
        :return: density
        :rtype: float
        """
        return self.total_mass(units)/self.total_volume(units)

    def wrap_coords(self):
        """wrap_coords wrap all coordinates into center periodic box
        """
        self.Coordinates.wrap_coords()

    def inherit_grx_attributes_from_molecules(self,molecule_dict,initial_composition,globally_unique=GRX_GLOBALLY_UNIQUE,unset_defaults=GRX_UNSET_DEFAULTS,overall_default=0):
        """inherit_grx_attributes_from_molecules Copy non-Gromacs-standard atom attributes in list "attributes" from molecule templates in molecule_dict according to molecule counts in dict initial_composition.

        :param attributes: list of labels of attributes to copy
        :type attributes: list
        :param molecule_dict: dictionary of available molecules (name:Molecule)
        :type molecule_dict: dict
        :param initial_composition: dictionary of initial composition (name:count)
        :type initial_composition: dict
        :param globally_unique: boolean list indicating attributes that must be globally unique, defaults to []
        :type globally_unique: list, optional
        :param unset_defaults: list of values, one per attribute, that signify UNSET, defaults to []
        :type unset_defaults: list, optional
        :param overall_default: default UNSET value for all attributes if unset_defaults is empty, defaults to 0
        :type overall_default: int, optional
        """
        # logger.debug(f'inherit grx {attributes} {unset_defaults} {globally_unique}')
        ''' set up the globally_unique and unset_defaults list if necessary '''
        if len(globally_unique)!=len(self.grxattr):
            globally_unique=[False for _ in range(len(self.grxattr))]
        if len(unset_defaults)!=len(self.grxattr):
            unset_defaults=[overall_default for _ in range(len(self.grxattr))]

        ''' drop all attribute values from current global atom dataframe '''
        adf=self.Coordinates.A
        self.Coordinates.A=adf.drop(columns=[a for a in self.grxattr if a in adf])

        logger.debug(f'{adf.shape[0]} atoms inheriting these attributes from molecular templates:')
        logger.debug(f'    Attribute name     Default value   Globally unique?')
        for attname,defval,gu in zip(self.grxattr,unset_defaults,globally_unique):
            logger.debug(f'    {attname:<15s}    {str(defval):<13s}   {gu}')
        if len(unset_defaults)==len(self.grxattr):
            att_running_maxval={}
            for k,v in zip(self.grxattr,unset_defaults):
                if type(v)==int or type(v)==float:
                    att_running_maxval[k]=0
                else:
                    att_running_maxval[k]='0'
        else:
            att_running_maxval={k:overall_default for k in self.grxattr}

        attribute_lists={k:[] for k in self.grxattr}
        value_counts={k:0 for k in self.grxattr}
        mol_idx_counter=0
        for icdict in [cc for cc in initial_composition if 'count' in cc]:
            molecule=icdict['molecule']
            count=icdict['count']
            mol_adf=molecule_dict[molecule].TopoCoord.Coordinates.A
            logger.debug(f'Inheriting from {molecule}')
            for ln in mol_adf.head().to_string().split('\n'):
                 logger.debug(ln)
            mol_attr_df=mol_adf[self.grxattr]
            for molecule_number in range(count):
                for i,k in enumerate(self.grxattr):
                    tra=mol_attr_df[k].to_list()
                    if k=='molecule': 
                        utra=[mol_idx_counter for _ in tra]
                        mol_idx_counter+=1
                    else:
                    # nv=len(tra)-tra.count(unset_defaults[i])
                        nuv=len(list(set([x for x in tra if x != unset_defaults[i]])))
                        utra=[]
                        for x in tra:
                            if globally_unique[i] and (type(x)==int or type(x)==float):
                                xx=x+value_counts[k] if x!=unset_defaults[i] else unset_defaults[i]
                            else:
                                xx=x
                            utra.append(xx)
                        value_counts[k]+=nuv
                    attribute_lists[k].extend(utra)

        for k,L in attribute_lists.items():
            self.Coordinates.A[k]=L
        logger.debug(f'postinherit adf columns {self.Coordinates.A.columns}')

    def make_resid_graph(self,json_file=None):
        """make_resid_graph make a residue connectivity graph

        :param json_file: name of output JSON file to write, defaults to None
        :type json_file: str, optional
        """
        self.Topology.make_resid_graph(json_file=json_file)

    def maxspan(self):
        """maxspan Returns the maxspan of the Coordinates (dimensions of orthorhombic
            convex hull enclosing Coordinates). Just a wrapper.

        :return: array of x-span, y-span, z-span
        :rtype: numpy.ndarray
        """
        return self.Coordinates.maxspan()

    def minmax(self):
        """minmax returns the coordinates of the atoms at the lower-leftmost and upper-rightmost positions in the constellation of points in the atoms dataframe

        :return: tuple of two points, lower-leftmost and upper-rightmost, respectively
        :rtype: tuple(np.ndarray(3,float),np.ndarray(3,float))
        """
        return self.Coordinates.minmax()

    def checkbox(self):
        """checkbox checks that the entire constellation of points in the atoms dataframe fits within the designated box for this Configuration object

        :return: True,True if both lower-leftmost and upper-rightmost points are within the box
        :rtype: tuple(bool,bool)
        """
        return self.Coordinates.checkbox()

    def write_mol2(self,filename,molname='',element_names_as_types=False):
        """write_mol2 Writes a SYBYL MOL2-format file using Coordinates, with certain
           atom attributes borrowed from the Topology

        :param filename: name of file to write
        :type filename: str
        :param molname: name of molecule to put in mol2 file, defaults to ''
        :type molname: str, optional
        """
        if molname=='':
            molname='This Molecule has no name'
        other_attributes=pd.DataFrame()
        if element_names_as_types:
            element_names=[x[0] for x in self.Topology.D['atoms']['atom']]
            other_attributes['type']=element_names
        else:
            other_attributes['type']=self.Topology.D['atoms']['type']
        other_attributes['charge']=self.Topology.D['atoms']['charge']
        self.files['mol2']=os.path.abspath(filename)
        # logger.debug(f'write_mol2, other_attributes:\n{other_attributes.to_string()}')
        if 'mol2_bonds' in self.Topology.D:
            self.Coordinates.write_mol2(filename,molname=molname,bondsDF=self.Topology.D['mol2_bonds'],other_attributes=other_attributes)
        else:
            self.Coordinates.write_mol2(filename,molname=molname,other_attributes=other_attributes)

    def merge(self,other):
        """merge Merges the TopoCoord instance "other" to self

        :param other: another TopoCoord instance
        :type other: TopoCoord
        :return: a shift tuple (returned by Coordinates.merge())
        :rtype: tuple
        """
        self.Topology.merge(other.Topology)
        shifts=self.Coordinates.merge(other.Coordinates)
        for name,idx_lists in other.idx_lists.items():
            # logger.debug(f'TopoCoord merge: list_name {name} lists {idx_lists}')
            for a_list in idx_lists:
                self.idx_lists[name].append([x+shifts[0] for x in a_list])
            self.reset_grx_attributes_from_idx_list(name)
        return shifts

    def bondtest_df(self,df:pd.DataFrame,pbc=[1,1,1],show_piercings=True):
        """bondtest_df applies bond filters to all bonds in the dataframe;

        :param df: dataframe of possible bonds; dataframe should have columns 'ai', 'aj' and 'r'; this method adds the column 'results'
        :type df: pd.DataFrame
        :param pbc: flags indicating dimensions in which pbc are used, defaults to [1,1,1]
        :type pbc: list, optional
        :param show_piercings: toggles diagnostic output for pierced rings, defaults to True
        :type show_piercings: bool, optional
        :return: input data frame with new 'results' column
        :rtype: pandas.DataFrame
        """
        if df.empty:
            return df
        results=[]
        for r in df.itertuples():
            result,dummy=self.bondtest((r.ai,r.aj,r.r),pbc=pbc,show_piercings=show_piercings)
            results.append(result)
        df['result']=results
        for ln in df[df['result']==BTRC.passed].to_string().split('\n'):
            logger.debug(ln)
        return df

    def bondtest(self,b,pbc=[1,1,1],show_piercings=True):
        """Determine if bond b is to be allowed to form based on geometric and
            topological criteria

        :param b: bond, tuple (ai,aj,rij)
        :type b: 2-tuple
        :param pbc: periodic boundary condition flags in each direction, defaults to [1,1,1]
        :type pbc: list, optional
        :param show_piercings: flag indicating you want to write gro files showing
            pierced rings
        :type show_piercings: bool
        :return: BTRC instance
        :rtype: BTRC enum
        """
        i,j,rij=b
        # i=int(i)
        # j=int(j)
        if self.makes_shortcircuit(i,j):
            return BTRC.failed_short_circuit,0
        if self.makes_cycle(i,j):
            return BTRC.failed_bond_cycle,0
        if self.pierces_ring(i,j,pbc=pbc,show_piercings=show_piercings):
            return BTRC.failed_pierce_ring,0
        # logger.debug(f'passed: {i:>7d} {j:>7d} {rij:>6.3f} nm')
        return BTRC.passed,rij

    def pierces_ring(self,i,j,pbc=[1,1,1],show_piercings=True):
        """pierces_ring checks to see if bond i-j would pierce any covalent ring structure

        :param i: global index of an atom
        :type i: int
        :param j: global index of another
        :type j: int
        :param pbc: flags indicating which dimensions have pbc applied, defaults to [1,1,1]
        :type pbc: list, optional
        :param show_piercings: toggles diagnostic output of ring piercings, defaults to True
        :type show_piercings: bool, optional
        :return: True if a ring is pierced; False otherwise
        :rtype: bool
        """
        adf=self.Coordinates.A
        LC=self.Coordinates.linkcell
        # at the current state, a linkcell is active under Coordinates
        # with spacing *greater* than the initial length of any bond.
        # so we can visit rings with one or more atom in a cell neighboring
        # the cells of the two atoms
        assert 'linkcell_idx' in adf,f'Error: atoms have no linkcell_idx attribute - bug!'
        i_lcidx=self.get_gro_attribute_by_attributes('linkcell_idx',{'globalIdx':i})
        j_lcidx=self.get_gro_attribute_by_attributes('linkcell_idx',{'globalIdx':j})
        joint_idx=[]
        for idx in LC.neighborlists[i_lcidx]+LC.neighborlists[j_lcidx]+[i_lcidx,j_lcidx]:
            if not idx in joint_idx:
                joint_idx.append(idx)
        cycle_tags=list(set(adf[(adf['cycle']!=-1)&(adf['linkcell_idx'].isin(joint_idx))]['cycle'].to_list()))
        B=adf.iloc[[i-1,j-1]].copy()
        saveB=B.copy()
        # logger.debug(f'bond {i} {j} subframe')
        # for ln in B.to_string().split('\n'):
        #     logger.debug(ln)
        for c in cycle_tags:
            C=adf[adf['cycle']==c].copy()
            # logging.debug(f'ring {c} has {C.shape[0]} members')
            # assert 4<C.shape[0]<7
            saveC=C.copy()
            if self.Coordinates.pierces(B,C,pbc):
                if show_piercings:
                    sub=self.Coordinates.subcoords(pd.concat([B,C]))
                    sub.write_gro(f'ring-{i}-{j}'+'.gro')
                    sub=self.Coordinates.subcoords(pd.concat([saveB,saveC]))
                    sub.write_gro(f'ring-orig-{i}-{j}'+'.gro')
                    logger.debug(f'Cycle {c} pierced by bond candidate {i}-{j}')
                    for ln in sub.A.to_string().split('\n'):
                        logger.debug(ln)
                return True
        return False

    def makes_shortcircuit(self,i,j):
        """Determine whether atoms i and j, if bonded, would produce a short circuit,
           defined as an instance in which i and j belong to residues that are already
           bonded to each other

        :param i: global index of first atom
        :type i: int
        :param j: global index of second atom
        :type j: int
        :return: True if a short circuit would happen, False otherwise
        :rtype: bool
        """
        i_resName,i_resNum,i_atomName,i_molNum=self.get_gro_attribute_by_attributes(['resName','resNum','atomName','molecule'],{'globalIdx':i})
        j_resName,j_resNum,j_atomName,j_molNum=self.get_gro_attribute_by_attributes(['resName','resNum','atomName','molecule'],{'globalIdx':j})
        '''
        In a cure reaction, atoms that react should be in different residues
        '''
        assert i_resNum!=j_resNum,f'shortcircuit test error {i}-{j} both in residue {i_resNum}?'
        i_neighbors=self.partners_of(i)
        j_neighbors=self.partners_of(j)
        # if i in j_neighbors or j in i_neighbors:
        #     # logger.debug(f'atoms {i} and {j} already on each other\'s list of bonded partners')
        #     return True
        '''
        Should never test a proposed bond that already exists
        '''
        assert not j in i_neighbors
        assert not i in j_neighbors
        '''
        Set up a DataFrame for heavy atoms in each molecule
        '''
        ADF=self.Coordinates.A
        R1DF=ADF[ADF['molecule']==i_molNum]
        R1DF=R1DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R1DF['atomName']]]
        R2DF=ADF[ADF['molecule']==j_molNum]
        R2DF=R2DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R2DF['atomName']]]
        assert R1DF.shape[0]>0
        assert R2DF.shape[0]>0
        '''
        Test to see if there exists any bond between these two molecules
        '''
        for i,j in product(R1DF['globalIdx'].to_list(),R2DF['globalIdx'].to_list()):
            if self.are_bonded(i,j):
                return True
        return False

    def reset_grx_attributes_from_idx_list(self,list_name):
        """reset_grx_attributes_from_idx_list uses information in the "index lists" to repopulate appropriate GRX attributes.  There are two index lists:  one for cycles and the other for chains.  Each index list is a list of lists; each element corresponds to a unique structure (cycle or chain) and is a list of global atom indices for atoms that make up that structural instance.

        :param list_name: either 'cycle' or 'chain'
        :type list_name: str
        """
        self.set_gro_attribute(list_name,-1)
        self.set_gro_attribute(f'{list_name}_idx',-1)
        for i,c in enumerate(self.idx_lists[list_name]):
            for j,x in enumerate(c):
                self.set_gro_attribute_by_attributes(list_name,i,{'globalIdx':x})
                self.set_gro_attribute_by_attributes(f'{list_name}_idx',j,{'globalIdx':x})

    def reset_idx_list_from_grx_attributes(self,list_name):
        """reset_idx_list_from_grx_attributes is the inverse of reset_grx_attributes_from_idx_list: it uses the GRX attributes to rebuild index lists.

        :param list_name: either 'cycle' or 'chain'
        :type list_name: str
        """
        adf=self.Coordinates.A
        logger.debug(f'reset: columns {adf.columns}')
        tmp_dict={}
        for i,r in adf.iterrows():
            gix=r['globalIdx']
            cid=r[list_name]
            cix=r[f'{list_name}_idx']
            if cid!=-1:
                if not cid in tmp_dict:
                    tmp_dict[cid]={}
                tmp_dict[cid][cix]=gix
        # logger.debug(f'{list_name} tmp_dict item count: {len(tmp_dict)}')
        # logger.debug(f'{tmp_dict}')
        if tmp_dict:
            consec_test=[a in tmp_dict for a in range(len(tmp_dict))]
            # logger.debug(f'{consec_test}')
            assert all(consec_test),f'{list_name} reset_idx_list for group attribute {list_name} has non-consecutive integer keys -- bug\n{[a in tmp_dict for a in range(len(tmp_dict))]}'
            ngroups=len(tmp_dict)
            self.idx_lists[list_name]=[[] for _ in range(ngroups)]
            for i in range(ngroups):
                for j in range(len(tmp_dict[i])):
                    self.idx_lists[list_name][i].append(tmp_dict[i][j])
        # logger.debug(f'-> idx_lists[{list_name}]: {self.idx_lists[list_name]}')

    def chainlist_update(self,new_bond_recs,msg=''):
        """chainlist_update updates the chain index lists due to generation of new bonds

        :param new_bond_recs: list of bond records, each a tuple of two ints corresponding to atom indices
        :type new_bond_recs: list
        :param msg: a nice message (unused), defaults to ''
        :type msg: str, optional
        """
        chainlists=self.idx_lists['chain']
        if len(chainlists)==0: return
        # logger.debug(f'pre {msg} chainlists')
        # for i,c in enumerate(chainlists):
        #     logger.debug(f'  {i} {c}')
        for b in new_bond_recs:
            aidx,bidx=b[0],b[1]
            ar=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':aidx})
            br=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bidx})
            if ar==br: continue # ignore intramolecular bonds
            # logger.debug(f'chainlist_update pair {aidx} {bidx}')
            ac=self.get_gro_attribute_by_attributes('chain',{'globalIdx':aidx})
            bc=self.get_gro_attribute_by_attributes('chain',{'globalIdx':bidx})
            logger.debug(f'ac {ac} bc {bc}')
            if ac==-1 or bc==-1:
                # neither of these newly bonded atoms is already in a chain, so
                # there is no possibility that this new bond can join two chains.
                continue
            # logger.debug(f'chain of bidx {bidx}: {chainlists[bc]}')
            # logger.debug(f'chain of aidx {aidx}: {chainlists[ac]}')
            aci=self.get_gro_attribute_by_attributes('chain_idx',{'globalIdx':aidx})
            bci=self.get_gro_attribute_by_attributes('chain_idx',{'globalIdx':bidx})
            logger.debug(f' -> {aidx}-{bidx}: ac {ac} #{len(chainlists[ac])} aci {aci} bc {bc} #{len(chainlists[bc])} bci {bci}')
            # one must be a head and the other a tail
            if aci==0: # a is a head
                if len(chainlists[bc])-1!=bci:
                    logger.warning(f'Atom {bidx} has index {bci} in chain {bc} but that chain has {len(chainlists[bc])} elements, so it cannot be the tail of the chain. The tail appears to be atom {chainlists[bc][-1]}. So something is wrong, and I am not merging these chains.')
                    logger.debug(f'chain of {aidx} {chainlists[ac]}')
                    logger.debug(f'chain of {bidx} {chainlists[bc]}')
                    continue
                c1=ac
                c2=bc
            elif bci==0: # b is a head
                if len(chainlists[ac])-1!=aci:
                    logger.warning(f'Atom {aidx} has index {aci} in chain {ac} but that chain has {len(chainlists[ac])} elements, so it cannot be the tail of the chain. The tail appears to be atom {chainlists[ac][-1]}. So something is wrong, and I am not merging these chains.')
                    logger.debug(f'chain of {aidx} {chainlists[ac]}')
                    logger.debug(f'chain of {bidx} {chainlists[bc]}')
                    continue
                c1=bc
                c2=ac
            else:
                logger.debug(f'Neither of {aidx} or {bidx} are heads of their respective chains.')
                logger.debug(f'This is most likely happening because of branching, which is permitted but')
                logger.debug(f'does not allow for merging of 1-D chains.')
                continue
            chainlists[c2].extend(chainlists[c1])
            for aidx in chainlists[c1]:
                self.set_gro_attribute_by_attributes('chain',c2,{'globalIdx':aidx})
                self.set_gro_attribute_by_attributes('chain_idx',chainlists[c2].index(aidx),{'globalIdx':aidx})
            # logger.debug(f'removing chain {c1}')
            chainlists.remove(chainlists[c1])
            # since we remove c1, all indices greater than c1 must decrement
            dec_us=np.array(self.Coordinates.A['chain'])
            bad_chain_idx=np.where(dec_us>c1)
            # logger.debug(f'bad_chain_idx: {bad_chain_idx}')
            # logger.debug(f'{dec_us[bad_chain_idx]}')
            dec_us[bad_chain_idx]-=1
            self.Coordinates.A['chain']=dec_us
        cnms=[]
        for c in self.idx_lists['chain']:
            cnms.append([self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in c])
        # logger.debug(f'post {msg} chains {self.idx_lists["chain"]} {cnms}')

    def makes_cycle(self,aidx,bidx):
        """makes_cycle checks the current chain index lists to see if a bond between aidx and bidx (global atom indices) would generate a C-C' cycle

        :param aidx: global index of an atom
        :type aidx: int
        :param bidx: global index of another atom
        :type bidx: int
        :return: True if a bond betwen aidx and bidx would form a C-C' cycle
        :rtype: bool
        """
        # is there a chain with aidx as head and bidx as tail, or vice versa?
        for c in self.idx_lists['chain']:
            if (aidx==c[0] and bidx==c[-1]) or (aidx==c[-1] and bidx==c[0]):
                return True
        return False

    def cycle_collective(self,bdf:pd.DataFrame):
        """cycle_collective Check to see if, when considered as a collective, this
        set of bondrecs leads to one or more cyclic chain; if so, longest bonds that 
        break cycles are removed from bondrecs list and resulting list is returned

        :param bdf: bonds data frame
        :type bdf: pandas.DataFrame
        :return: list of bond records that results in no new cycles
        :rtype: list
        """
        def addlink(chainlist,i,j):
            i_cidx=(-1,-1)
            j_cidx=(-1,-1)
            for cidx,c in enumerate(chainlist):
                if i in c:
                    i_cidx=(cidx,c.index(i))
                if j in c:
                    j_cidx=(cidx,c.index(j))
            logger.debug(f'i {i_cidx} j {j_cidx}')
            if not (i_cidx==(-1,-1) or (j_cidx==(-1,-1))):
                # logger.debug(f'cycle_collective: {i} is in chain {i_cidx[0]} at {i_cidx[1]}')
                # logger.debug(f'cycle_collective: {j} is in chain {j_cidx[0]} at {j_cidx[1]}')
                if i_cidx[0]==j_cidx[0]:
                    # cyclized!
                    logger.debug(f'chain {i_cidx[0]} is cyclized!')
                    return i_cidx[0],True
                old_head,old_tail=(i_cidx,j_cidx) if i_cidx[1]==0 else (j_cidx,i_cidx)
                chainlist[old_tail[0]].extend(chainlist[old_head[0]])
                # logger.debug(f'new chain {old_tail[0]} {chainlist[old_tail[0]]}')
                chainlist.remove(chainlist[old_head[0]])
                return old_tail[0],False
            else: # this bond involves one atom in a chain and another that is not; no way this can be a cyclization, but we have to return somthing
                return i_cidx[0] if i_cidx[0]!=-1 else j_cidx[0],False

        # kb=bondrecs.copy()
        # logger.debug(f'checking set of {bdf.shape[0]} bonds for nascent cycles')
        new_bdf=bdf.copy()
        new_bdf['remove-to-uncyclize']=[False for _ in range(new_bdf.shape[0])]
        chains=deepcopy(self.idx_lists['chain'])
        assert id(chains) != id(self.idx_lists['chain'])
        if not chains:
            return new_bdf
        cyclized_chains={}
        chains_of_bonds=[]
        for i,r in bdf.iterrows():
            chain_idx,cyclized=addlink(chains,r.ai,r.aj)
            # logger.debug(f'bond {i} {r.ai}-{r.aj} ({r.reactantName}) adds to chain {chain_idx}')
            chains_of_bonds.append(chain_idx)
            if cyclized and not chain_idx in cyclized_chains:
                cyclized_chains[chain_idx]=[]
        assert len(chains_of_bonds)==bdf.shape[0]
        for i,r in bdf.iterrows():
            # logger.debug(f'compiling bonds of cyclized chains: bond {i}')
            cidx=chains_of_bonds[i]
            # logger.debug(f'compiling bonds of cyclized chains: cidx {chains_of_bonds[i]}')
            if cidx in cyclized_chains:
                cyclized_chains[cidx].append(i)
        for k,v in cyclized_chains.items():
            lb=v[-1]  # last bondrec in list
            new_bdf.loc[lb,'remove-to-uncyclize']=True
        return new_bdf

    def get_bystanders(self,atom_idx):
        """get_bystanders identify and return bystanders at a particular proposed bond specified by atom_idx

        :param atom_idx: container of two global atom indices specifying a proposed bond
        :type atom_idx: container
        :return: 4-tuple of special bystander info containers
        :rtype: tuple
        """
        bystander_atomidx=[[int],[int]] # elements correspond to atoms in atom_idx
        bystander_atomnames=[[str],[str]]
        bystander_resids=[[int],[int]]
        bystander_resnames=[[str],[str]]
        for x in [0,1]:
            bystander_atomidx[x]=self.interresidue_partners_of(atom_idx[x])
            if atom_idx[1-x] in bystander_atomidx[x]: # if new bond has not yet formed
                bystander_atomidx[x].remove(atom_idx[1-x])
            bystander_atomnames[x]=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':y}) for y in bystander_atomidx[x]]
            bystander_resids[x]=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':y}) for y in bystander_atomidx[x]]
            bystander_resnames[x]=[self.get_gro_attribute_by_attributes('resName',{'globalIdx':y}) for y in bystander_atomidx[x]]
        return bystander_resids,bystander_resnames,bystander_atomidx,bystander_atomnames
    
    def get_oneaways(self,atom_idx):
        """get_oneaways identify and return one-aways for a particular proposed bond specified by atom_idx

        :param atom_idx: two-element container of atom indices specifying proposed bond
        :type atom_idx: list
        :return: tuple of oneaways-lists
        :rtype: tuple
        """
        resids=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in atom_idx]
        chains=[self.get_gro_attribute_by_attributes('chain',{'globalIdx':x}) for x in atom_idx]
        chain_idx=[self.get_gro_attribute_by_attributes('chain_idx',{'globalIdx':x}) for x in atom_idx]
        # logger.debug(f'chains {chains} chain_idx {chain_idx}')
        oneaway_resids=[None,None]
        oneaway_resnames=[None,None]
        oneaway_atomidx=[None,None]
        oneaway_atomnames=[None,None]
        # assert chain_idx[0]==0 or chain_idx[1]==0 # one must be a head
        if chains!=[-1,-1]:
            chainlists_idx=[self.idx_lists['chain'][chains[x]] for x in [0,1]]
            chainlists_atomnames=[]
            chainlists_resids=[]
            chainlists_resnames=[]
            for x in [0,1]:
                atomnamelist=[]
                residlist=[]
                resnamelist=[]
                for y in chainlists_idx[x]:
                    rnum,rname,aname=self.get_gro_attribute_by_attributes(['resNum','resName','atomName'],{'globalIdx':y})
                    residlist.append(rnum)
                    resnamelist.append(rname)
                    atomnamelist.append(aname)
                chainlists_atomnames.append(atomnamelist)
                chainlists_resids.append(residlist)
                chainlists_resnames.append(resnamelist)
            # logger.debug(f'chainlists_idx: {chainlists_idx}')
            # logger.debug(f'chainlists_atomnames: {chainlists_atomnames}')
            # logger.debug(f'chainlists_resids: {chainlists_resids}')
            # logger.debug(f'chainlists_resnames: {chainlists_resnames}')
            oa_chain_idx=[]
            if chain_idx[0]==0 or chain_idx[0]>chain_idx[1]:
                oa_chain_idx=[chain_idx[0]+2,chain_idx[1]-2]
            elif chain_idx[1]==0 or chain_idx[1]>chain_idx[0]:
                oa_chain_idx=[chain_idx[0]-2,chain_idx[1]+2]
            # logger.debug(f'oa_chain_idx {oa_chain_idx}')
            for i in range(len(oa_chain_idx)):
                if oa_chain_idx[i]<0:
                    oa_chain_idx[i]=None
            oa_idx=[None,None]
            for i in [0,1]:
                try:
                    oa_idx[i]=chainlists_idx[i][oa_chain_idx[i]]
                except:
                    oa_idx[i]=None
                if oa_idx[i]:
                    oneaway_atomidx[i]=chainlists_idx[i][oa_chain_idx[i]]
                    oneaway_atomnames[i]=chainlists_atomnames[i][oa_chain_idx[i]]
                    oneaway_resids[i]=chainlists_resids[i][oa_chain_idx[i]]
                    oneaway_resnames[i]=chainlists_resnames[i][oa_chain_idx[i]]
        # logger.debug(f'result oneaway_resids {oneaway_resids} oneaway_resnames {oneaway_resnames}')
        # logger.debug(f'oneaway_atomidx {oneaway_atomidx} oneaway_atomnames {oneaway_atomnames}')
        return oneaway_resids,oneaway_resnames,oneaway_atomidx,oneaway_atomnames

    def grab_files(self):
        """grab_files using absolute pathname information, grab the most up-to-date gromacs files for this system and deposit them into the cwd
        """
        cwd=os.getcwd()
        for ext in ['top','gro','grx']:
            filename=self.files[ext]
            logger.debug(f'{ext} grabbing {pfs.proj_abspath(filename)} into {pfs.cwd()}')  
            if os.path.commonprefix([cwd,filename])!=cwd:
                shutil.copy(filename,cwd)
            self.files[ext]=os.path.abspath(os.path.basename(filename))

    def grompp_and_mdrun(self,out,mdp,mylogger=logger.debug,**kwargs):
        """grompp_and_mdrun manages invoking a single Gromacs run using the current TopoCoord

        :param out: output filename basename
        :type out: str
        :param mdp: name of mdp file
        :type mdp: str
        :param mylogger: a logger, defaults to logger.debug
        :type mylogger: logging.logger, optional
        :return: any message generated by gromacs.grompp_and_mdrun
        :rtype: str
        """
        wrap_coords=kwargs.get('wrap_coords',False)
        self.grab_files()
        # make sure required files exist in this directory
        top=os.path.basename(self.files['top']).replace('.top','')
        gro=os.path.basename(self.files['gro']).replace('.gro','')
        assert os.path.exists(f'{top}.top')
        assert os.path.exists(f'{gro}.gro')
        assert os.path.exists(f'{mdp}.mdp')
        mylogger(f'TopoCoord: running Gromacs {top}.top, {gro}.gro, {mdp}.mdp')
        msg=grompp_and_mdrun(gro=gro,top=top,out=out,mdp=mdp,**kwargs)
        self.copy_coords(TopoCoord(grofilename=f'{out}.gro',wrap_coords=wrap_coords))
        mylogger(f'after grompp_and_run: gro {self.files["gro"]}')
        return msg

    def load_files(self,filenames:dict):
        """load_files load all gromacs files into TopoCoord

        :param filenames: dictionary of extension:filename
        :type filenames: dict
        """
        for e,n in filenames.items():
            if not e in ['gro','top','grx','mol2']: continue
            bn,ext=os.path.splitext(n)
            logger.debug(f'bn {bn} ext {ext}')
            if   ext=='.gro':  self.read_gro(n)
            elif ext=='.top':  self.read_top(n)
            elif ext=='.grx':  self.read_gro_attributes(n)
            elif ext=='.mol2': self.read_mol2(n)
            else:
                logger.debug(f'Warning: file {n} has unknown file extension.  Skipped.')

    def center_coords(self,new_boxsize:np.ndarray=None):
        """center_coords center all coordinates in box

        :param new_boxsize: new boxsize if desired, defaults to None
        :type new_boxsize: numpy.ndarray, optional
        """
        if type(new_boxsize)==np.ndarray:
            if new_boxsize.shape==(3,):
                box_vectors=new_boxsize*np.identity(3,dtype=float)
            elif new_boxsize.shape==(3,3):
                box_vectors=new_boxsize
            self.Coordinates.box=box_vectors
        center=self.Coordinates.box.diagonal()/2.0
        gc=self.Coordinates.geometric_center()
        addme=center-gc
        self.Coordinates.translate(addme)

    def vacuum_minimize(self,outname='minimized',**kwargs):
        """vacuum_minimize the minimize analog to grompp_and_mdrun; performs an energy minimization using mdrun

        :param outname: output file basename, defaults to 'minimized'
        :type outname: str, optional
        """
        pad=kwargs.get('pad',5)
        boxsize=np.array(self.maxspan())+pad*np.ones(3)
        logger.debug(f'{self.maxspan()} -> {boxsize}')
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix='single-molecule-min'
        pfs.checkout(f'mdp/{mdp_prefix}.mdp')
        # gromacs_dict={'nt':1,'nb':'cpu','pme':'cpu','pmefft':'cpu','bonded':'cpu','update':'cpu'}
        self.grompp_and_mdrun(out=f'{outname}',
            mdp=mdp_prefix,boxSize=boxsize,single_molecule=True,wrap_coords=False) #,**gromacs_dict)

    def vacuum_simulate(self,outname='simulated',**kwargs):
        """vacuum_simulate peform a vacuum MD simulation using mdrun

        :param outname: output file basename, defaults to 'simulated'
        :type outname: str, optional
        """
        pad=kwargs.get('pad',5)
        boxsize=np.array(self.maxspan())+pad*np.ones(3)
        logger.debug(f'{self.maxspan()} -> {boxsize}')
        self.center_coords(new_boxsize=boxsize)
        mdp_prefix='single-molecule-nvt'
        pfs.checkout(f'mdp/{mdp_prefix}.mdp')
        nsamples=kwargs.get('nsamples',10)
        params=kwargs.get('params',{})
        T=params.get('temperature',300.0)
        ps=params.get('ps',0.0)
        nsteps=params.get('nsteps',-1)
        begin_at=params.get('begin_at',0.0)
        assert ps!=0.0 or nsteps!=-1
        dt=float(mdp_get(f'{mdp_prefix}.mdp','dt'))
        if ps!=0.0:
            nsteps=int(float(ps)/dt)
        else:
            ps=nsteps*dt
        sampling_duration_ps=ps-begin_at
        sampling_duration_nsteps=int(sampling_duration_ps/dt)
        sample_interval=sampling_duration_nsteps//nsamples
        logger.debug(f'nsteps {nsteps} begin_at {begin_at} ps {ps} sampling_duration_ps {sampling_duration_ps} ({sampling_duration_nsteps} steps) sample_interval {sample_interval}')
        # nsteps=nsamples*(sample_interval+1)
        # nsteps=mdp_get(f'{mdp_prefix}.mdp','nsteps')
        # nstxout=mdp_get(f'{mdp_prefix}.mdp','nstxout')
        mdp_modify(f'{mdp_prefix}.mdp',{'nsteps':nsteps,'nstxout':sample_interval,'ref_t':T})
        # gromacs_dict={'nt':1,'nb':'cpu','pme':'cpu','pmefft':'cpu','bonded':'cpu','update':'cpu'}
        self.grompp_and_mdrun(out=f'{outname}',
            mdp=mdp_prefix,boxSize=boxsize,single_molecule=True,wrap_coords=False) #,**gromacs_dict)

    def equilibrate(self,deffnm='equilibrate',edict={},gromacs_dict={}):
        """equilibrate perform an MD simulation using mdrun

        :param deffnm: output file basename, defaults to 'equilibrate'
        :type deffnm: str, optional
        :param edict: dictionary of simulation directives, defaults to {}
        :type edict: dict, optional
        :param gromacs_dict: dictionary of gromacs directives, defaults to {}
        :type gromacs_dict: dict, optional
        :return: list of edr files this equilibration generates
        :rtype: list
        """
        mod_dict={}
        edr_list=[]
        ens=edict['ensemble']
        repeat=edict.get('repeat',0)
        assert ens in ['min','npt','nvt'],f'Bad ensemble: {ens}'
        pfs.checkout(f'mdp/{ens}.mdp') # plain jane?
        if ens=='min':
            msg='minimization'
        else:
            msg=f'{ens} ensemble'
        if ens in ['nvt','npt']:
            dt=float(mdp_get(f'{ens}.mdp','dt'))
            mod_dict['ref_t']=edict['temperature']
            mod_dict['gen-temp']=edict['temperature']
            mod_dict['gen-vel']='yes'
            nsteps=edict.get('nsteps',0)
            if not nsteps:
                ps=edict.get('ps',0.0)
                if not ps: return
                nsteps=int(float(ps)/dt)
            else:
                ps=nsteps*dt
            mod_dict['nsteps']=nsteps
            msg+=f'; {ps:6.2f} ps, {float(edict["temperature"]):7.2f} K'
            if ens=='npt': 
                mod_dict['ref_p']=edict['pressure']
                msg+=f', {float(edict["pressure"]):6.2f} bar'
        new_mdp=f'{deffnm}-{ens}'
        mdp_modify(f'{ens}.mdp',mod_dict,new_filename=f'{new_mdp}.mdp')
        logger.info(f'Running Gromacs: {msg}')
        self.grompp_and_mdrun(out=f'{deffnm}-{ens}',mdp=new_mdp,quiet=False,**gromacs_dict)
        edr_list=[f'{deffnm}-{ens}']
        if ens=='npt':
            box=self.Coordinates.box.diagonal()
            logger.info(f'Current box side lengths: {box[0]:.3f} nm x {box[1]:.3f} nm x {box[2]:.3f} nm')
            gmx_energy_trace(f'{deffnm}-{ens}',['Density'],report_averages=True,**gromacs_dict)
        for rep in range(repeat):
            logger.info(f'Repeat {rep+1} out of {repeat}')
            this_deffnm=f'{deffnm}-repeat-{rep+1}-{ens}'
            self.grompp_and_mdrun(out=this_deffnm,mdp=new_mdp,quiet=False,**gromacs_dict)
            edr_list.append(this_deffnm)
            if ens=='npt':
                box=self.Coordinates.box.diagonal()
                logger.info(f'Current box side lengths: {box[0]:.3f} nm x {box[1]:.3f} nm x {box[2]:.3f} nm')
            gmx_energy_trace(this_deffnm,['Density'],report_averages=True,**gromacs_dict)
        return edr_list

    def get_resid_sets(self,atom_pair):
        """get_resid_sets identifies individual sets of separate resids owned by unbonded atoms i and j 

        :param atom_pair: tuple of i,j atom indices
        :type atom_pair: tuple
        :return: tuple containing the two apposing residue sets
        :rtype: tuple
        """
        # assertion: i and j are not bonded and they represent two separate sets of residues in this topocoord.
        i,j=atom_pair
        ri=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':i})
        rj=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':j})
        ci=self.Topology.local_resid_cluster(ri)
        cj=self.Topology.local_resid_cluster(rj)
        if set(ci)==set(cj):
            logger.debug(f'resid sets overlap: {i} {ri} : {ci} and {j} {rj} : {cj}')
            return []
        return [ci,cj]

    def check_your_topology(self):
        """check_your_topology checks topology for duplicate 1-4 pair interactions and deletes them
        """
        T=self.Topology
        C=self.Coordinates
        aT={}
        pdf=T.D['pairs']
        aT['dx']=[]
        aT['dy']=[]
        aT['dz']=[]
        checked=[]
        drops=[]
        for i,r in pdf.iterrows():
            ai,aj=min([r.ai,r.aj]),max([r.ai,r.aj])
            assert ai!=aj
            if not (ai,aj) in checked: 
                checked.append((ai,aj))
            else:
                # print(f'repeated pair {ai} {aj}')
                drops.append(i)
            ri=self.get_R(ai)
            rj=self.get_R(aj)
            D=C.mic(ri-rj,pbc=[1,1,1])
            aT['dx'],aT['dy'],aT['dz']=D
        print(f'{len(drops)} duplicate pairs detected')
        T.D['pairs']=pdf.drop(drops)
        pdf['dx']=aT['dx']
        pdf['dy']=aT['dy']
        pdf['dz']=aT['dz']
        print(pdf.sort_values(by='dx').head(3).to_string())
        print(pdf.sort_values(by='dy').head(3).to_string())
        print(pdf.sort_values(by='dz').head(3).to_string())
        self.write_top('checked.top')

def find_template(BT:BondTemplate,moldict):
    """find_template searches the dictionary of available molecules to identify a bond template that matches the passed-in template, returning the corresponding template molecule and reaction-bond

    :param BT: bond template to search for
    :type BT: BondTemplate
    :param moldict: dictionary of available molecule
    :type moldict: MoleculeDict
    :raises Exception: if no matching template is found
    :return: template Molecule object, corresponding ReactionBond object from that template, and a boolean flag indicating whether or not the match required a symmetric-reversal of the template
    :rtype: tuple(Molecule,ReactionBond,bool)
    """
    use_T=None
    b_idx=-1
    reverse_bond=False
    for template_name,T in moldict.items():
        # if BT in T.bond_templates:
        #     use_T=T
        #     break
        use_T=None
        for b_idx in range(len(T.bond_templates)):
            bt=T.bond_templates[b_idx]
            # logger.debug(f'comparing to {T.name} {str(bt)}')
            if bt==BT:
                # logger.debug(f'TRUE')
                use_T=T
                break
            if bt.is_reverse_of(BT):
                # logger.debug(f'TRUE, but...')
                use_T=T
                reverse_bond=True
                # logger.debug(f'reversed bond instance {str(RB)}')
        if use_T!=None:
            break
    else:
        logger.error(f'No template is found for {str(BT)}')
        raise Exception('you have a bond for which I cannot find a template')
    logger.debug(f'Using template {use_T.name} and bond index {b_idx}')
    rb=use_T.reaction_bonds[b_idx]
    return use_T,rb,reverse_bond

