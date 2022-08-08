"""

.. module:: topocoords
   :synopsis: Class for jointly handling Topology and Coordinate instances

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

#
#  class for methods that need to work with both Topology and Coordinates
from itertools import product
import pandas as pd
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topology import Topology
from HTPolyNet.bondtemplate import BondTemplate,ReactionBond
from HTPolyNet.gromacs import grompp_and_mdrun
# from HTPolyNet.molecule import MoleculeDict,ReactionList
#from HTPolyNet.plot import network_graph
import HTPolyNet.projectfilesystem as pfs
import logging
import numpy as np
from enum import Enum
import os
import shutil
from copy import deepcopy

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
        wrap_coords=kwargs.get('wrap_coords',True)
        self.files={}
        self.files['gro']=grofilename
        self.files['top']=topfilename
        self.files['grx']=grxfilename
        self.files['mol2']=mol2filename
        self.grxattr=[]
        self.idx_lists={}
        self.idx_lists['chain']=[]
        self.idx_lists['cycle']=[]
        if grofilename!='':
            self.read_gro(grofilename,wrap_coords=wrap_coords)
            if grxfilename!='':
                self.grxattr=self.read_gro_attributes(grxfilename)
        else:
            self.Coordinates=Coordinates()  # empty
        if topfilename!='':
            self.read_top(topfilename)
        else:
            self.Topology=Topology(system_name=system_name) # empty
        if mol2filename!='':
            self.read_mol2(mol2filename)
        if grxfilename!='':
            self.read_gro_attributes(grxfilename)

    def make_bonds(self,pairs,skip_H=[]):
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
        idx_to_ignore=self.Coordinates.find_sacrificial_H(pairs,self.Topology,skip_pairs=skip_H)
        # logger.debug(f'idx_to_ignore {idx_to_ignore}')
        self.Topology.add_bonds(pairs)
        self.chainlist_update(pairs,msg='TopoCoord.make_bonds')
        self.Topology.null_check(msg='add_bonds')
        rename=True if len(skip_H)>0 else False
        idx_to_delete=self.Coordinates.find_sacrificial_H(pairs,self.Topology,skip_pairs=skip_H,rename=rename)
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
            resids=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in bb]
            # this is the product name of the reaction used to identify this bond
            product_name=b.reactantName
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
            use_T=None
            b_idx=-1
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
                        RB.reverse()
                        # logger.debug(f'reversed bond instance {str(RB)}')
                if use_T!=None:
                    break
            else:
                logger.error(f'No template is found for {bb}: {str(BT)}')
                raise Exception('you have a bond for which I cannot find a template')
            
            T=use_T
            logger.debug(f'Using template {T.name} and bond index {b_idx}')
            rb=T.reaction_bonds[b_idx]
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
            logger.debug(f'inst2temp {inst2temp}')
            logger.debug(f'temp2inst {temp2inst}')
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
                    logger.debug(f'changing charge of inst atom {inst_atom} ({inst_resn} {inst_rnam} {inst_name}) from {inst_charge} to {temp_charge}')
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
            logger.debug(f'Template angles:')
            for ln in temp_angles.to_string().split('\n'):
                logger.debug(ln)
            inst_angles=temp_angles.copy()
            inst_angles.ai=temp_angles.ai.map(temp2inst)
            inst_angles.aj=temp_angles.aj.map(temp2inst)
            inst_angles.ak=temp_angles.ak.map(temp2inst)
            logger.debug(f'Mapped instance angles:')
            for ln in inst_angles.to_string().split('\n'):
                logger.debug(ln)
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

            logger.debug(f'Template dihedrals:')
            for ln in temp_dihedrals.to_string().split('\n'):
                logger.debug(ln)
            # map from template atom indicies to system atom indicies in dihedrals
            inst_dihedrals=temp_dihedrals.copy()
            inst_dihedrals.ai=temp_dihedrals.ai.map(temp2inst)
            inst_dihedrals.aj=temp_dihedrals.aj.map(temp2inst)
            inst_dihedrals.ak=temp_dihedrals.ak.map(temp2inst)
            inst_dihedrals.al=temp_dihedrals.al.map(temp2inst)
            logger.debug(f'Mapped instance dihedrals:')
            for ln in inst_dihedrals.to_string().split('\n'):
                logger.debug(ln)
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
            logger.debug(f'Concatenating this pairs to global pairs')
            for ln in temp_pairs.to_string().split('\n'):
                logger.debug(ln)
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

    def update_topology_and_coordinates(self,bdf,template_dict={},write_mapper_to=None,**kwargs):
        """update_topology_and_coordinates updates global topology and necessary atom attributes in the configuration to reflect formation of all bonds listed in "keepbonds"

        :param bdf: bonds dataframe, columns 'ai', 'aj', 'reactantName'
        :type bdf: pandas.DataFrame
        :param template_dict: dictionary of molecule templates keyed on molecule name
        :type dict of molecules
        :return: 3-tuple: new topology file name, new coordinate file name, list of bonds with atom indices updated to reflect any atom deletions
        :rtype: 3-tuple
        """
        overcharge_threshhold=kwargs.get('overcharge_threshhold',0.1)
        logger.debug(f'begins.')
        if bdf.shape[0]>0:
            assert bdf['ai'].dtype==int
            assert bdf['aj'].dtype==int
            # pull out just the atom index pairs (first element of each tuple)
            at_idx=[(int(x.ai),int(x.aj),x.order) for x in bdf.itertuples()]
            logger.debug(f'Making {len(at_idx)} bonds.')
            idx_to_delete=self.make_bonds(at_idx)
            logger.debug(f'Deleting {len(idx_to_delete)} atoms.')
            idx_mapper=self.delete_atoms(idx_to_delete) # will result in full reindexing
            logger.debug(f'null check')
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
            logger.debug(f'calling map_from_templates')
            self.map_from_templates(ri_bdf,template_dict,overcharge_threshhold=overcharge_threshhold)
            logger.debug(f'1-4 pair update')
            # enumerate ALL pairs involving either or both of the bonded atoms
            pdf=self.Topology.D['pairs']
            # each of these bonds results in 1-4 pair interactions
            bl=self.Topology.bondlist
            pai=[]
            paj=[]
            pri=[]
            prj=[]
            prsource=[]
            for p in at_idx:
                j,k=p
                # jresnum,jresname,jname=self.get_gro_attribute_by_attributes(['resNum','resName','atomName'],{'globalIdx':j})
                # kresnum,kresname,kname=self.get_gro_attribute_by_attributes(['resNum','resName','atomName'],{'globalIdx':k})
                nj=bl.partners_of(j)
                nj.remove(k)
                nk=bl.partners_of(k)
                nk.remove(j)
                this_pairs=list(product(nj,nk))
                pai.extend([x[0] for x in this_pairs])
                paj.extend([x[1] for x in this_pairs])
                pri.extend([j for x in this_pairs])
                prj.extend([k for x in this_pairs])
                prsource.extend(['c' for x in this_pairs])
                jpdf=pdf[(pdf['ai']==j)|(pdf['aj']==j)].copy()
                # logger.debug(f'j: pairs with member {j} {jresnum} {jresname} {jname}')
                # for ln in jpdf.to_string().split('\n'):
                #     logger.debug(ln)
                pai.extend(jpdf['ai'].to_list())
                paj.extend(jpdf['aj'].to_list())
                pri.extend([j for x in jpdf['ai'].to_list()])
                prj.extend([k for x in jpdf['aj'].to_list()])
                prsource.extend(['l' for _ in range(jpdf.shape[0])])
                kpdf=pdf[(pdf['ai']==k)|(pdf['aj']==k)].copy()
                # logger.debug(f'k: pairs with member {k} {kresnum} {kresname} {kname}')
                # for ln in kpdf.to_string().split('\n'):
                #     logger.debug(ln)
                pai.extend(kpdf['ai'].to_list())
                paj.extend(kpdf['aj'].to_list())
                pri.extend([j for x in kpdf['ai'].to_list()])
                prj.extend([k for x in kpdf['aj'].to_list()])
                prsource.extend(['r' for _ in range(kpdf.shape[0])])
            pi_df=pd.DataFrame({'ai':pai,'aj':paj,'source':prsource,'bi':pri,'bj':prj})
            pi_df.drop_duplicates(inplace=True,ignore_index=True)
            self.Topology.null_check(msg='update_topology_and_coordinates')
            if write_mapper_to:
                tdf=pd.DataFrame({'old':list(idx_mapper.keys()),'new':list(idx_mapper.values())})
                tdf.to_csv(write_mapper_to,sep=' ',index=False)
            logger.debug('finished')
            return ri_bdf,pi_df

    def read_top(self,topfilename):
        """Creates a new Topology member by reading from a Gromacs-style top file.
            Just a wrapper for the read_gro method of Topology

        :param topfilename: name of topology file
        :type topfilename: str
        """
        self.files['top']=os.path.abspath(topfilename)
        self.Topology=Topology.read_gro(topfilename)

    def read_gro(self,grofilename,preserve_box=False,wrap_coords=True):
        """Creates a new Coordinates member by reading from a Gromacs-style coordinates
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

    def read_mol2(self,mol2filename,ignore_bonds=False,overwrite_coordinates=False):
        """Creates a new Coordinates member by reading from a SYBYL-style MOL2 file.
            A wrapper for read_mol2 from Coordinates, but also sets the 'mol2_bonds'
            dataframe in the Topology if the parameter ignore_bonds is False.  If
            the mol2_bonds dataframe is created, and the Topology already has a 'bonds' dataframe, a consistency check is peformed.

        :param mol2filename: name of mol2 file
        :type mol2filename: str
        :param ignore_bonds: flag indicating bond records in the mol2 file should be
            ignored, defaults to False
        :type ignore_bonds: bool, optional
        """
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
        """Swaps the names of the two atoms with global indicies ai and aj.  This is used when
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
        """Wrapper for read_top and read_gro; generates new Topology and Coordinates members

        :param topfilename: name of topology file
        :type topfilename: str
        :param grofilename: name of coordinates file (Gromacs format)
        :type grofilename: str
        """
        self.read_top(topfilename)
        self.read_gro(grofilename)

    def write_top(self,topfilename):
        """Write a Gromacs-format topology file; this will only write an in-line version,
            no itp; wrapper for Topology.to_file()

        :param topfilename: name of file to write
        :type topfilename: str
        """
        self.Topology.to_file(topfilename)
        self.files['top']=os.path.abspath(topfilename)

    def write_gro(self,grofilename):
        """Write a Gromacs-format coordinate file; wrapper for Coordinates.write_gro()

        :param grofilename: name of file to write
        :type grofilename: str
        """
        self.Coordinates.write_gro(grofilename)
        self.files['gro']=os.path.abspath(grofilename)

    def write_top_gro(self,topfilename,grofilename):
        """Writes both a Gromacs top file and Gromacs coordinate file

        :param topfilename: name of topology file to write
        :type topfilename: str
        :param grofilename: name of coordinate file to write
        :type grofilename: str
        """
        self.write_top(topfilename)
        self.write_gro(grofilename)

    def return_bond_lengths(self,bdf):
        """Return the length of all bonds in list bonds

        :param bdf: bonds dataframe, 'ai','aj','reactantName'
        :type bonds: pandas.DataFrame
        :return: list of lengths parallel to bonds
        :rtype: list of floats
        """
        return self.Coordinates.return_bond_lengths(bdf)

    def add_length_attribute(self,bdf:pd.DataFrame,attr_name='length'):
        self.Coordinates.add_length_attribute(bdf,attr_name=attr_name)

    # def return_pair_lengths(self):
    #     """Return the length of all 1-4 pairs in topology

    #     :param bdf: bonds dataframe, 'ai','aj','reactantName'
    #     :type bonds: pandas.DataFrame
    #     :return: list of lengths parallel to bonds
    #     :rtype: list of floats
    #     """
    #     return self.Coordinates.return_pair_lengths(self.Topology.D['pairs'])

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
        self.files['gro']=other.files['gro']

    def restore_bond_parameters(self,saved):
        """Retores saved bond parameters in df saved by overwriting

        :param saved: [ bonds ] dataframe
        :type saved: pandas.DataFrame
        """
        self.Topology.restore_bond_parameters(saved)

    def set_grx_attributes(self,attributes):
        self.grxattr=attributes
    
    def write_gro_attributes(self,attributes_list,grxfilename):
        """Writes atomic attributes to a file

        :param attributes_list: list of attributes to write
        :type attributes_list: list
        :param grxfilename: name of output file
        :type grxfilename: str
        """
        self.Coordinates.write_atomset_attributes(attributes_list,grxfilename)
        self.files['grx']=os.path.abspath(grxfilename)

    def write_grx_attributes(self,grxfilename):
        self.write_gro_attributes(self.grxattr,grxfilename)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
        """Read attributes from file into self.Coordinates.A

        :param grxfilename: name of input file
        :type grxfilename: str
        :param attribute_list: list of attributes to take, defaults to [] (take all)
        :type attribute_list: list, optional
        """
        self.files['grx']=os.path.abspath(grxfilename)
        attributes_read=self.Coordinates.read_atomset_attributes(grxfilename,attributes=attribute_list)
        if 'chain' in attributes_read and 'chain_idx' in attributes_read:
            self.reset_idx_list_from_grx_attributes('chain')
        if 'cycle' in attributes_read and 'cycle_idx' in attributes_read:
            self.reset_idx_list_from_grx_attributes('cycle')
        if attributes_read!=self.grxattr:
            self.grxattr=attributes_read

    def set_gro_attribute(self,attribute,srs):
        self.Coordinates.set_atomset_attribute(attribute,srs)

    def set_gro_attribute_by_attributes(self,att_name,att_value,attribute_dict):
        self.Coordinates.set_atom_attribute(att_name,att_value,attribute_dict)

    def get_gro_attribute_by_attributes(self,att_name,attribute_dict):
        return self.Coordinates.get_atom_attribute(att_name,attribute_dict)

    def increment_gro_attribute_by_attributes(self,att_name,attribute_dict):
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        # logger.debug(f'increment {att_name} {attribute_dict} {val}')
        val+=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def decrement_gro_attribute_by_attributes(self,att_name,attribute_dict):
        val=self.get_gro_attribute_by_attributes(att_name,attribute_dict)
        # logger.debug(f'decrement {att_name} {attribute_dict} {val}')
        val-=1
        self.set_gro_attribute_by_attributes(att_name,val,attribute_dict)

    def get_gro_attributelist_by_attributes(self,attribute_list,attribute_dict):
        return self.Coordinates.get_atoms_w_attribute(attribute_list,attribute_dict)

    def get_R(self,idx):
        return self.Coordinates.get_R(idx)

    def rotate(self,R):
        self.Coordinates.rotate(R)

    def translate(self,L):
        self.Coordinates.translate(L)

    def partners_of(self,i):
        return self.Topology.bondlist.partners_of(i)

    def resid_partners_of(self,ri):
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
        result=[]
        bl=self.Topology.bondlist.partners_of(i)
        logger.debug(f'{i} partners {bl}')
        myresid=self.Coordinates.A.iloc[i-1]['resNum']
        for j in bl:
            theirresid=self.Coordinates.A.iloc[j-1]['resNum']
            if theirresid!=myresid:
                result.append(j)
        return result

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        return self.Coordinates.minimum_distance(other.Coordinates,self_excludes=self_excludes,other_excludes=other_excludes)

    # def ring_detector(self):
    #     return self.Topology.ring_detector()

    def has_gro_attributes(self,attribute_list):
        return self.Coordinates.has_atom_attributes(attribute_list)

    def are_bonded(self,i,j):
        return self.Topology.bondlist.are_bonded(i,j)

    def decrement_z(self,pairs):
        self.Coordinates.decrement_z(pairs)

    def show_z_report(self):
        self.Coordinates.show_z_report()

    def make_ringlist(self):
        self.Coordinates.make_ringlist()

    def adjust_charges(self,atoms=[],overcharge_threshhold=0.1,netcharge=0.0,msg=''):
        self.Topology.adjust_charges(atoms=atoms,overcharge_threshhold=overcharge_threshhold,desired_charge=netcharge,msg=msg)

    def gro_DataFrame(self,name):
        if name=='atoms':
            return self.Coordinates.A
        elif name=='mol2_bonds':
            return self.Coordinates.mol2_bonds
        else:
            return None

    def overwrite_coords(self,other):
        logger.debug(f'Overwriting {other.Coordinates.A.shape[0]} coordinates')
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

    def label_cycle_atoms(self):
        adf=self.Coordinates.A
        cycle_idx=list(sorted(list(set(adf['cycle'].to_list()))))
        if -1 in cycle_idx:
            cycle_idx.remove(-1)
        for c in cycle_idx:
            cnames=adf[adf['cycle']==c]['atomName'].to_list()
            logger.debug(f'cycle {c}: {cnames}')
        # logger.debug(f'label_ring_atoms for {self.name}:\n{adf.to_string()}')

    def linkcell_initialize(self,cutoff,ncpu=1,force_repopulate=False):
        """Initialize the linkcell structure; a wrapper for Coordinates

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
        """Check to be sure the Coordinate and Topology members contain the same number of
            atoms

        :return: the number of atoms
        :rtype: int
        """
        assert self.Coordinates.A.shape[0]==self.Topology.D['atoms'].shape[0]
        return self.Coordinates.A.shape[0]

    def total_mass(self,units='SI'):
        """Returns the total mass of the system.  Just a wrapper.

        :param units: units designation, defaults to 'SI' (other option is 'gromacs')
        :type units: str, optional
        :return: mass
        :rtype: float
        """
        return self.Topology.total_mass(units=units)

    def wrap_coords(self):
        self.Coordinates.wrap_coords()

    def inherit_grx_attributes_from_molecules(self,molecule_dict,initial_composition,globally_unique=[],unset_defaults=[],overall_default=0):
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
        for icdict in initial_composition:
            molecule=icdict['molecule']
            count=icdict['count']
            mol_adf=molecule_dict[molecule].TopoCoord.Coordinates.A
            for ln in mol_adf.to_string().split('\n'):
                logger.debug(ln)
            mol_attr_df=mol_adf[self.grxattr]
            for i in range(count):
                for i,k in enumerate(self.grxattr):
                    tra=mol_attr_df[k].to_list()
                    # nv=len(tra)-tra.count(unset_defaults[i])
                    nuv=len(list(set([x for x in tra if x != unset_defaults[i]])))
                    utra=[]
                    for x in tra:
                        if globally_unique[i] and (type(x)==int or type(x)==float):
                            xx=x+value_counts[k] if x!=unset_defaults[i] else unset_defaults[i]
                        else:
                            xx=x
                        utra.append(xx)
                    attribute_lists[k].extend(utra)
                    value_counts[k]+=nuv

        for k,L in attribute_lists.items():
            self.Coordinates.A[k]=L
        # logger.debug(f'postinherit adf columns {self.Coordinates.A.columns}')

    # def make_ringlist(self):
    #     self.Coordinates.make_ringlist()

    def make_resid_graph(self,json_file=None,draw=None):
        self.Topology.make_resid_graph(json_file=json_file,draw=draw)

    def maxspan(self):
        """Returns the maxspan of the Coordinates (dimensions of orthorhombic
            convex hull enclosing Coordinates). Just a wrapper.

        :return: array of x-span, y-span, z-span
        :rtype: numpy.ndarray
        """
        return self.Coordinates.maxspan()


    def minmax(self):
        return self.Coordinates.minmax()

    def checkbox(self):
        return self.Coordinates.checkbox()

    def write_mol2(self,filename,molname=''):
        """Writes a SYBYL MOL2-format file using Coordinates, with certain
           atom attributes borrowed from the Topology

        :param filename: name of file to write
        :type filename: str
        :param molname: name of molecule to put in mol2 file, defaults to ''
        :type molname: str, optional
        """
        if molname=='':
            molname='This Molecule has no name'
        other_attributes=pd.DataFrame()
        other_attributes['type']=self.Topology.D['atoms']['type']
        other_attributes['charge']=self.Topology.D['atoms']['charge']
        self.files['mol2']=os.path.abspath(filename)
        # logger.debug(f'write_mol2, other_attributes:\n{other_attributes.to_string()}')
        if 'mol2_bonds' in self.Topology.D:
            self.Coordinates.write_mol2(filename,molname=molname,bondsDF=self.Topology.D['mol2_bonds'],other_attributes=other_attributes)
        else:
            self.Coordinates.write_mol2(filename,molname=molname,other_attributes=other_attributes)

    def merge(self,other):
        """Merges the TopoCoord instance "other" to self

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

    # def bondtest_par(self,B,pbc=[1,1,1],show_piercings=True):
    #     """Parallelization of bondtest

    #     :param B: list of bonds (generated by a split to be mapped onto a Pool)
    #     :type B: list of 2-tuples
    #     :param pbc: periodic boundary condition flags in each direction, defaults to [1,1,1]
    #     :type pbc: list, optional
    #     :param show_piercings: flag indicating you want to write gro files showing
    #         pierced rings
    #     :type show_piercings: bool
    #     :return: list of booleans, each is True if that bond is permitted
    #     :rtype: list
    #     """
    #     L=[]
    #     for b in B:
    #         L.append(self.bondtest(b,pbc=pbc,show_piercings=show_piercings))
    #     return L

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
        adf=self.Coordinates.A
        LC=self.Coordinates.linkcell
        # Ri=self.get_R(i)
        # Rj=self.get_R(j)
        # Rij=self.Coordinates.mic(Ri-Rj,pbc)
        # chk_rij=np.sqrt(Rij.dot(Rij))
        # assert np.isclose(rij,chk_rij,atol=1.e-3),f'{i} {j} {rij} {chk_rij}'
        # generate the nearest periodic image of Rj to Ri
        # Rjp=Ri-Rij
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
        i_resName,i_resNum,i_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i})
        j_resName,j_resNum,j_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j})
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
        Set up a DataFrame for heavy atoms in each residue
        '''
        ADF=self.Coordinates.A
        R1DF=ADF[ADF['resNum']==i_resNum]
        R1DF=R1DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R1DF['atomName']]]
        R2DF=ADF[ADF['resNum']==j_resNum]
        R2DF=R2DF[[(not (b.startswith('H') or b.startswith('h'))) for b in R2DF['atomName']]]
        assert R1DF.shape[0]>0
        assert R2DF.shape[0]>0
        '''
        Test to see if there exists any bond between these two residues
        '''
        for i,j in product(R1DF['globalIdx'].to_list(),R2DF['globalIdx'].to_list()):
            if self.are_bonded(i,j):
                return True
        return False

    def reset_grx_attributes_from_idx_list(self,list_name):
        self.set_gro_attribute(list_name,-1)
        self.set_gro_attribute(f'{list_name}_idx',-1)
        for i,c in enumerate(self.idx_lists[list_name]):
            for j,x in enumerate(c):
                self.set_gro_attribute_by_attributes(list_name,i,{'globalIdx':x})
                self.set_gro_attribute_by_attributes(f'{list_name}_idx',j,{'globalIdx':x})

    def reset_idx_list_from_grx_attributes(self,list_name):
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
        logger.debug(f'{list_name} tmp_dict item count: {len(tmp_dict)}')
        logger.debug(f'{tmp_dict}')
        if tmp_dict:
            consec_test=[a in tmp_dict for a in range(len(tmp_dict))]
            logger.debug(f'{consec_test}')
            assert all(consec_test),f'{list_name} reset_idx_list for group attribute {list_name} has non-consecutive integer keys -- bug\n{[a in tmp_dict for a in range(len(tmp_dict))]}'
            ngroups=len(tmp_dict)
            self.idx_lists[list_name]=[[] for _ in range(ngroups)]
            for i in range(ngroups):
                for j in range(len(tmp_dict[i])):
                    self.idx_lists[list_name][i].append(tmp_dict[i][j])
        logger.debug(f'-> idx_lists[{list_name}]: {self.idx_lists[list_name]}')

    # def remap_idx_list(self,list_name,mapper):
    #     logger.debug(f'{list_name}')
    #     remapped_groups=[]
    #     for c in self.idx_lists[list_name]:
    #         remapped_groups.append([mapper[x] for x in c])
    #     self.idx_lists[list_name]=remapped_groups
    #     self.reset_grx_attributes_from_idx_list(list_name)
    #     logger.debug(f'finished.')

    def chainlist_update(self,new_bond_recs,msg=''):
        chainlists=self.idx_lists['chain']
        logger.debug(f'pre {msg} chains')
        for i,c in enumerate(chainlists):
            logger.debug(f'  {i} {c}')
        for b in new_bond_recs:
            aidx,bidx=b[0],b[1]
            logger.debug(f'chainlist_update pair {aidx} {bidx}')
            ac=self.get_gro_attribute_by_attributes('chain',{'globalIdx':aidx})
            bc=self.get_gro_attribute_by_attributes('chain',{'globalIdx':bidx})
            if ac in chainlists:
                logger.debug(f'chain of aidx {aidx}: {chainlists[ac]}')
            if bc in chainlists:
                logger.debug(f'chain of bidx {bidx}: {chainlists[bc]}')
            if ac==-1 or bc==-1:
                # neither of these newly bonded atoms is already in a chain, so
                # there is no possibility that this new bond can join two chains.
                continue
            aci=self.get_gro_attribute_by_attributes('chain_idx',{'globalIdx':aidx})
            bci=self.get_gro_attribute_by_attributes('chain_idx',{'globalIdx':bidx})
            logger.debug(f' -> {aidx}-{bidx}: ac {ac} bc {bc} aci {aci} bci {bci}')
            # one must be a head and the other a tail
            if aci==0: # a is a head
                assert len(chainlists[bc])-1==bci,f'incorrect tail'
                c1=ac
                c2=bc
            elif bci==0: # b is a head
                assert len(chainlists[ac])-1==aci,f'incorrect tail'
                c1=bc
                c2=ac
            else:
                raise Exception(f'chain parse error')
            chainlists[c2].extend(chainlists[c1])
            for aidx in chainlists[c1]:
                self.set_gro_attribute_by_attributes('chain',c2,{'globalIdx':aidx})
                self.set_gro_attribute_by_attributes('chain_idx',chainlists[c2].index(aidx),{'globalIdx':aidx})
            chainlists.remove(chainlists[c1])
            # since we remove c1, all indices greater than c1 must decrement
            dec_us=np.array(self.Coordinates.A['chain'])
            bad_chain_idx=np.where(dec_us>c1)
            dec_us[bad_chain_idx]-=1
            self.Coordinates.A['chain']=dec_us

        cnms=[]
        for c in self.idx_lists['chain']:
            cnms.append([self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in c])
        # logger.debug(f'post {msg} chains {self.idx_lists["chain"]} {cnms}')

    def makes_cycle(self,aidx,bidx):
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
        bystander_atomidx=[[int],[int]]
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
        # resids=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in atom_idx]
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
        cwd=os.getcwd()
        for ext in ['top','gro']:
            filename=self.files[ext]
            logger.debug(f'{ext} grabbing {filename} for {cwd}')  
            if os.path.commonprefix([cwd,filename])!=cwd:
                shutil.copy(filename,cwd)
            self.files[ext]=os.path.abspath(os.path.basename(filename))

    def grompp_and_mdrun(self,out,mdp,**kwargs):
        self.grab_files()
        top=os.path.basename(self.files['top']).replace('.top','')
        gro=os.path.basename(self.files['gro']).replace('.gro','')
        assert os.path.exists(f'{top}.top')
        assert os.path.exists(f'{gro}.gro')
        assert os.path.exists(f'{mdp}.mdp')
        logger.debug(f'{os.getcwd()} {pfs.cwd()} {top}, {gro}, {mdp}')
        msg=grompp_and_mdrun(gro=gro,top=top,out=out,mdp=mdp,**kwargs)
        self.copy_coords(TopoCoord(grofilename=f'{out}.gro'))
        logger.debug(f'after grommp_and_run: gro {self.files["gro"]}')
        return msg

    def load_files(self):
        if self.files['gro']:
            self.read_gro(self.files['gro'])
        if self.files['top']:
            self.read_top(self.files['top'])
        if self.files['grx']:
            self.read_gro_attributes(self.files['grx'])
        if self.files['mol2']:
            self.read_mol2(self.files['mol2'])
