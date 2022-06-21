"""

.. module:: topocoords
   :synopsis: Class for jointly handling Topology and Coordinate instances
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

#
#  class for methods that need to work with both Topology and Coordinates
from itertools import product
from operator import is_
from tokenize import maybe
import pandas as pd
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topology import Topology
import logging
import numpy as np
import networkx as nx
from enum import Enum

class BTRC(Enum):
    """Bond test return codes: bond tests are applied to those bond-candidates that are within search radius of each other

    :param Enum: inherits from Enum class
    :type Enum: class
    """
    passed = 0
    failed_pierce_ring = 1        # does this bond-candidate pierce a ring?
    failed_short_circuit = 2      # does this bond-candidate result in short-circuit?
    failed_polyethylene_cycle = 3 # does this bond-candidate create a polyethylene cycle?

class TopoCoord:
    """Container for Topology and Coordinates, along with methods that 
        use either or both of them
    """
    def __init__(self,topfilename='',grofilename='',mol2filename='',system_name='htpolynet'):
        """Constructor method for TopoCoord.

        :param topfilename: name of Gromacs-format topology file (top), defaults to ''
        :type topfilename: str, optional
        :param grofilename: name of Gromacs-format coordinate file (gro), defaults to ''
        :type grofilename: str, optional
        :param mol2filename: name of SYBYL MOL2-format coordinate/bonds file, defaults to ''
        :type mol2filename: str, optional
        """
        if grofilename!='':
            self.grofilename=grofilename
            self.read_gro(grofilename)
        else:
            self.Coordinates=Coordinates()  # empty
        if topfilename!='':
            self.topfilename=topfilename
            self.read_top(topfilename)
        else:
            self.Topology=Topology(system_name=system_name) # empty
        if mol2filename!='':
            self.mol2filename=mol2filename
            self.read_mol2(mol2filename) 
            # will overwrite coords and add 'mol2_bonds' section to topology

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
        # logging.debug(f'idx_to_ignore {idx_to_ignore}')
        self.Topology.add_bonds(pairs)
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

    # def add_pairs(self,pairdf,kb=300000.0):
    #     """Adds a pair for each pair in the pairdf (['ai'],['aj'],['initial-distance'])
        
    #     :param pairdf: dataframe of pairs
    #     :type pardf: pandas.DataFrames
    #     """
    #     self.Topology.add_pairs(pairdf,kb=kb)

    def delete_atoms(self,atomlist):
        """Deletes atoms from both the Topology and Coordinates instances

        :param atomlist: list of global indexes of atoms to delete
        :type atomlist: list
        :return: old-to-new atom global index mapper dictionary resulting from reindexing
            remaining atoms to make sure global indexes are sequential
        :rtype: dict
        """
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper=self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper)==dict
        return idx_mapper

    def select_template(self,bondrec,moldict,reaction_list):
        def is_product(a,reaction_list):
            return a in [r.product for r in reaction_list]
        def is_reactant(a,reaction_list):
            reactants=[]
            for r in reaction_list:
                for v in r.reactants.values():
                    if not v in reactants:
                        reactants.append(v)
            return a in v
        def find_reaction(reaction_list,bondrec=[]):
            (ian,irn),(jan,jrn)=bondrec
            for R in reaction_list:
                # find the bond in this reaction that corresponds to the bond
                # implied by bondrec
                iar,jar=None,None
                for B in R.bonds:
                    a,b=B['atoms']
                    an,bn=R.atoms[a]['atom'],R.atoms[b]['atom']
                    if (ian,jan)==(an,bn):
                        iar,jar=a,b
                        break
                    elif (ian,jan)==(bn,an):
                        iar,jar=b,a
                        break
                if iar and jar:
                    # find the reactants for these two atoms
                    irk=R.atoms[iar]['reactant']
                    jrk=R.atoms[jar]['reactant']
                    # are these two reactants the reactants these atoms brought in?
                    if R.reactants[irk]==irn and R.reactants[jrk]==jrn:
                        return R
            
            return None
        
        ai,aj,prodname=bondrec
        # itn=self.get_gro_attribute_by_attributes('old_reactantName',{'globalIdx':ai})
        # jtn=self.get_gro_attribute_by_attributes('old_reactantName',{'globalIdx':aj})
        ian=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':ai})
        jan=self.get_gro_attribute_by_attributes('atomName',{'globalIdx':aj})

#        assert itn==jtn # both atoms will have the same reactantName

        n_of_ai=self.partners_of(ai)
        # logging.debug(f'ai {ai} n {n_of_ai}')
        n_of_ai.remove(aj)
        n_of_aj=self.partners_of(aj)
        # logging.debug(f'aj {aj} n {n_of_aj}')
        n_of_aj.remove(ai)
        an_of_n_ai=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in n_of_ai]
        an_of_n_aj=[self.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in n_of_aj]
        rn_of_iresid=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in n_of_ai]
        rn_of_jresid=[self.get_gro_attribute_by_attributes('resNum',{'globalIdx':x}) for x in n_of_aj]
        nitn=[self.get_gro_attribute_by_attributes('reactantName',{'globalIdx':x}) for x in n_of_ai]
        njtn=[self.get_gro_attribute_by_attributes('reactantName',{'globalIdx':x}) for x in n_of_aj]
        # b/c of removals above, don't expect jresid to be on iresid-neighbors
        # if jresid in rn_of_iresid:
        #     rn_of_iresid.remove(jresid)
        # if iresid in rn_of_jresid:
        #     rn_of_jresid.remove(iresid)

        in_prods=[]
        # logging.debug(f'ai {ai} {self.get_gro_attribute_by_attributes("atomName",{"globalIdx":ai})} template {self.get_gro_attribute_by_attributes("reactantName",{"globalIdx":ai})}:')
        for n,an,r,rn in zip(n_of_ai,an_of_n_ai,rn_of_iresid,nitn):
            # logging.debug(f'  {n} {an} {r} {rn} product? {is_product(rn,reaction_list)}')
            if is_product(rn,reaction_list):
                in_prods.append((n,an,r,rn))

        jn_prods=[]
        # logging.debug(f'aj {aj} {self.get_gro_attribute_by_attributes("atomName",{"globalIdx":aj})} template {self.get_gro_attribute_by_attributes("reactantName",{"globalIdx":aj})}:')
        for n,an,r,rn in zip(n_of_aj,an_of_n_aj,rn_of_jresid,njtn):
            # logging.debug(f'  {n} {an} {r} {rn} product? {is_product(rn,reaction_list)}')
            if is_product(rn,reaction_list):
                jn_prods.append((n,an,r,rn))

        # logging.debug(f'ai {ai} {self.get_gro_attribute_by_attributes("atomName",{"globalIdx":ai})} template {prodname} in_prods {in_prods}:')
        # logging.debug(f'aj {aj} {self.get_gro_attribute_by_attributes("atomName",{"globalIdx":aj})} template {prodname} jn_prods {jn_prods}:')

        if len(in_prods)==0 and len(jn_prods)==0:
            self.set_gro_attribute_by_attributes('reactantName',prodname,{'globalIdx':ai})
            self.set_gro_attribute_by_attributes('reactantName',prodname,{'globalIdx':aj})
            return moldict[prodname]
        elif len(in_prods)==0 and len(jn_prods)==1:
            irn=self.get_gro_attribute_by_attributes('old_reactantName',{'globalIdx':ai})
            jrn=jn_prods[0][3]
        elif len(in_prods)==1 and len(jn_prods)==0:
            jrn=self.get_gro_attribute_by_attributes('old_reactantName',{'globalIdx':aj})
            irn=in_prods[0][3]
        elif len(in_prods)==1 and len(jn_prods)==1:
            irn=in_prods[0][3]
            jrn=jn_prods[0][3]
        else:
            logging.error(f'ai {ai} in_prods {in_prods} aj {aj} {jn_prods}: too many reactant names')
            raise Exception('Cannot identify template product for this bond formation reaction')
        R=find_reaction(reaction_list,bondrec=[(ian,irn),(jan,jrn)])
        # logging.debug(f'template name {R.product} from reaction {R.name}')
        prodname=R.product
        # you should only do this if prodname is a reactant
        if is_reactant(prodname,reaction_list):
            self.set_gro_attribute_by_attributes('reactantName',prodname,{'globalIdx':ai})
            self.set_gro_attribute_by_attributes('reactantName',prodname,{'globalIdx':aj})
        return moldict[prodname]

    def map_from_templates(self,bdf,moldict,reaction_list):
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
        grodf=self.Coordinates.A
        grodf['old_reactantName']=grodf['reactantName'].copy()
        # first pass -- set all reactant names
        for i,b in bdf.iterrows():
            bb=[b['ai'],b['aj']]
            template_name=b['reactantName']
            self.set_gro_attribute_by_attributes('reactantName',template_name,{'globalIdx':bb[0]})
            self.set_gro_attribute_by_attributes('reactantName',template_name,{'globalIdx':bb[1]})
            i_resid=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bb[0]})
            j_resid=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bb[1]})
            logging.debug(f'Bond {bb} (resids {i_resid} and {j_resid}) template {template_name}')
        for i,b in bdf.iterrows():
            bb=[b['ai'],b['aj']]
            bbb=[b['ai'],b['aj'],b['reactantName']]
            T=self.select_template(bbb,moldict,reaction_list)  # select template based on neighbors
            # T=moldict[template_name]
            # i_resid=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bb[0]})
            # j_resid=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':bb[1]})
            # logging.debug(f'Bond {bb} (resids {i_resid} and {j_resid}) using template {T.name}')
            temp_atdf=T.TopoCoord.Topology.D['atoms']
            # get the bidirectional instance<->template mapping dictionaries
            inst2temp,temp2inst=T.idx_mappers(self,bb)
            # some hard checks on compatibility of the dicts
            assert len(inst2temp)==len(temp2inst)
            check=True
            for k,v in inst2temp.items():
                check = check and (k == temp2inst[v])
            for k,v in temp2inst.items():
                check = check and (k == inst2temp[v])
            assert check,f'Error: bidirectional dicts are incompatible; bug\n{inst2temp}\b{temp2inst}'
            # logging.debug(f'map_from_templates:inst2temp {inst2temp}')
            # logging.debug(f'map_from_templates:temp2inst {temp2inst}')
            i_idx,j_idx=bb # the actual atoms that just bonded
            # get their indicies in the template
            temp_i_idx=inst2temp[i_idx]
            temp_j_idx=inst2temp[j_idx]
            d=T.TopoCoord.Topology.D['bonds']
            # copy all bond records matching these two bonds; should be only one!!
            d=d[((d.ai==temp_i_idx)&(d.aj==temp_j_idx))|
                ((d.ai==temp_j_idx)&(d.aj==temp_i_idx))].copy()
            if d.shape[0]!=1:
                logging.error(f'map_from_templates using {T.name} is sent inst-bond {i_idx}-{j_idx} which is claimed to map to {temp_i_idx}-{temp_j_idx}, but no such unique bond is found:\n{T.TopoCoord.Topology.D["bonds"].to_string()}')
            assert d.shape[0]==1,f'Error: see log'

            # get all angle, dihedrals, and pairs from template that result from the existence of the specified bond
            temp_angles,temp_dihedrals,temp_pairs=T.get_angles_dihedrals((temp_i_idx,temp_j_idx))
            # logging.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            # map from template atom indicies to system atom indicies in angles
            inst_angles=temp_angles.copy()
            # logging.debug(f'temp_angles:\n{temp_angles.to_string()}')
            inst_angles.ai=temp_angles.ai.map(temp2inst)
            inst_angles.aj=temp_angles.aj.map(temp2inst)
            inst_angles.ak=temp_angles.ak.map(temp2inst)
            # logging.debug(f'Mapped inst_angles:\n{inst_angles.to_string()}')
            # add new angles to the system topology
            d=self.Topology.D['angles']
            self.Topology.D['angles']=pd.concat((d,inst_angles),ignore_index=True)
            # for all atoms identified in these angles, update types and charges if
            # required
            for a in ['ai','aj','ak']:
                for inst_atom,temp_atom in zip(inst_angles[a],temp_angles[a]):
                    assert inst_atom in atdf['nr'].values,f'Error: mapped atom {inst_atom} not found in [ atoms ]'
                    inst_type,inst_charge=atdf[atdf['nr']==inst_atom][['type','charge']].values[0]
                    temp_type,temp_charge=temp_atdf[temp_atdf['nr']==temp_atom][['type','charge']].values[0]
                    # logging.debug(f'ang temp {temp_atom} {temp_type} {temp_charge}')
                    # logging.debug(f'ang inst {inst_atom} {inst_type} {inst_charge}')
                    if inst_type!=temp_type:
                        # logging.debug(f'(angles){a} changing type of inst atom {inst_atom} from {inst_type} to {temp_type}')
                        atdf.loc[atdf['nr']==inst_atom,'type']=temp_type
                    if inst_charge!=temp_charge:
                        # logging.debug(f'(angles){a} changing charge of inst atom {inst_atom} from {inst_charge} to {temp_charge}')
                        atdf.loc[atdf['nr']==inst_atom,'charge']=temp_charge
            # hard check for any nan's in any atom index attribute in any angle
            d=self.Topology.D['angles']
            check=True
            for a in ['ai','aj','ak']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in angles')
                raise Exception

            # map from template atom indicies to system atom indicies in dihedrals            
            inst_dihedrals=temp_dihedrals.copy()
            inst_dihedrals.ai=temp_dihedrals.ai.map(temp2inst)
            inst_dihedrals.aj=temp_dihedrals.aj.map(temp2inst)
            inst_dihedrals.ak=temp_dihedrals.ak.map(temp2inst)
            inst_dihedrals.al=temp_dihedrals.al.map(temp2inst)
            d=inst_dihedrals
            check=False
            for a in ['ai','aj','ak','al']:
                check=check or d[a].isnull().values.any()
            if check:
                logging.error(f'a {a} NAN in dihedrals\n{inst_dihedrals.to_string()}')
                raise Exception
            # add new dihedrals to global topology
            d=self.Topology.D['dihedrals']
            self.Topology.D['dihedrals']=pd.concat((d,inst_dihedrals),ignore_index=True)
            # update any necessary atom types and charges
            for a in ['ai','aj','ak','al']:
                for inst_atom,temp_atom in zip(inst_dihedrals[a],temp_dihedrals[a]):
                    # logging.debug(f'a {a} inst_atom {inst_atom} temp_atom {temp_atom}')
                    inst_type,inst_charge=atdf[atdf['nr']==inst_atom][['type','charge']].values[0]
                    temp_type,temp_charge=temp_atdf[temp_atdf['nr']==temp_atom][['type','charge']].values[0]
                    # logging.debug(f'dih temp {temp_atom} {temp_type} {temp_charge}')
                    # logging.debug(f'dih inst {inst_atom} {inst_type} {inst_charge}')
                    if inst_type!=temp_type:
                        # logging.debug(f'(dihedrals){a} changing type of inst atom {inst_atom} from {inst_type} to {temp_type}')
                        atdf.loc[atdf['nr']==inst_atom,'type']=temp_type
                    if inst_charge!=temp_charge:
                        # logging.debug(f'(dihedrals){a} changing charge of inst atom {inst_atom} from {inst_charge} to {temp_charge}')
                        atdf.loc[atdf['nr']==inst_atom,'charge']=temp_charge
            # hard check for no nans
            d=self.Topology.D['dihedrals']
            check=True
            for a in ['ai','aj','ak','al']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in dihedrals')
                raise Exception
            d=self.Topology.D['pairs']
            check=True
            for a in ['ai','aj']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in pairs premapping')
                raise Exception

            # double-hard check to make sure pairs can be mapped
            k=np.array(list(temp2inst.keys()))
            v=np.array(list(temp2inst.values()))
            if any(np.isnan(k)):
                logging.error('null in temp2inst keys')
            if any(np.isnan(v)):
                logging.error('null in temp2inst values')
            # logging.debug(f'temp_pairs:\n{temp_pairs.to_string()}')
            isin=[not x in temp2inst for x in temp_pairs.ai]
            if any(isin):
                for ii,jj in enumerate(isin):
                    if jj:
                        logging.error(f'atom ai {temp_pairs.ai.iloc[ii]} not in temp2inst')
            isin=[not x in temp2inst for x in temp_pairs.aj]
            if any(isin):
                for ii,jj in enumerate(isin):
                    if jj:
                        logging.error(f'atom aj {temp_pairs.aj.iloc[ii]} not in temp2inst')

            # map all ai attributes of all template pairs to global ai
            temp_pairs.ai=temp_pairs.ai.map(temp2inst)
            # check AGAIN for nans (I am afraid of nans)
            if temp_pairs.ai.isnull().values.any():
                logging.error('NAN in pairs ai')
                raise Exception
            # map all aj attributes of all template pairs to global ai
            temp_pairs.aj=temp_pairs.aj.map(temp2inst)
            # check YET AGAIN for nans (eek!)
            if temp_pairs.aj.isnull().values.any():
                logging.error('NAN in pairs aj')
                raise Exception
            # add these pairs to the topology
            self.Topology.D['pairs']=pd.concat((d,temp_pairs),ignore_index=True)
            d=self.Topology.D['pairs']
            # check AGAIN for nans
            check=True
            for a in ['ai','aj']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in pairs post mapping')
                raise Exception

    def update_topology_and_coordinates(self,bdf,template_dict={},write_mapper_to=None,reaction_list=[]):    
        """update_topology_and_coordinates updates global topology and necessary atom attributes in the configuration to reflect formation of all bonds listed in "keepbonds"

        :param bdf: bonds dataframe, columns 'ai', 'aj', 'reactantName'
        :type bdf: pandas.DataFrame
        :param template_dict: dictionary of molecule templates keyed on molecule name
        :type dict of molecules
        :return: 3-tuple: new topology file name, new coordinate file name, list of bonds with atom indices updated to reflect any atom deletions
        :rtype: 3-tuple
        """
        # logging.debug(f'update_topology_and_coordinates begins.')
        if bdf.shape[0]>0:
            # pull out just the atom index pairs (first element of each tuple)
            at_idx=[(x['ai'],x['aj'],x['order']) for i,x in bdf.iterrows()]
            logging.debug(f'Making {len(at_idx)} bonds.')
            idx_to_delete=self.make_bonds(at_idx)
            # logging.debug(f'Deleting {len(idx_to_delete)} atoms.')
            idx_mapper=self.delete_atoms(idx_to_delete) # will result in full reindexing
            self.Topology.null_check(msg='delete_atoms')
            self.write_gro('tmp.gro')
            # reindex all atoms in the list of bonds sent in, and write it out
            ri_bdf=bdf.copy()
            ri_bdf.ai=ri_bdf.ai.map(idx_mapper)
            ri_bdf.aj=ri_bdf.aj.map(idx_mapper)
            self.Topology.update_polyethylenes(ri_bdf,idx_mapper)
            at_idx=[(x['ai'],x['aj']) for i,x in ri_bdf.iterrows()]

            self.decrement_z(at_idx)
            self.make_ringlist()
            self.map_from_templates(ri_bdf,template_dict,reaction_list)
            # TODO: enumerate ALL pairs involving either or both of the bonded atoms
            pdf=self.Topology.D['pairs']
            # each of these bonds results in 1-4 pair interactions 
            bl=self.Topology.bondlist
            pai=[]
            paj=[]
            for p in at_idx:
                j,k=p
                nj=bl.partners_of(j)
                nj.remove(k)
                nk=bl.partners_of(k)
                nk.remove(j)
                this_pairs=list(product(nj,nk))
                pai.extend([x[0] for x in this_pairs])
                paj.extend([x[1] for x in this_pairs])
                jpdf=pdf[(pdf['ai']==j)|(pdf['aj']==j)].copy()
                pai.extend(jpdf['ai'].to_list())
                paj.extend(jpdf['aj'].to_list())
                kpdf=pdf[(pdf['ai']==k)|(pdf['aj']==k)].copy()
                pai.extend(kpdf['ai'].to_list())
                paj.extend(kpdf['aj'].to_list())
            pi_df=pd.DataFrame({'ai':pai,'aj':paj})
            pi_df.drop_duplicates(inplace=True,ignore_index=True)
            self.Topology.null_check(msg='map_from_templates')
            self.adjust_charges(msg='')
            if write_mapper_to:
                with open(write_mapper_to,'w') as f:
                    for k,v in idx_mapper.items():
                        f.write(f'{k} {v}\n')
            return ri_bdf,pi_df

    def read_top(self,topfilename):
        """Creates a new Topology member by reading from a Gromacs-style top file.
            Just a wrapper for the read_gro method of Topology

        :param topfilename: name of topology file
        :type topfilename: str
        """
        self.Topology=Topology.read_gro(topfilename)

    def read_gro(self,grofilename):
        """Creates a new Coordinates member by reading from a Gromacs-style coordinates
            file.  Just a wrapper for the read_gro method of Coordinates

        :param grofilename: name of gro file
        :type grofilename: str
        """
        self.Coordinates=Coordinates.read_gro(grofilename)

    def read_mol2(self,mol2filename,ignore_bonds=False):
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
        self.Coordinates=Coordinates.read_mol2(mol2filename)
        if not ignore_bonds:
            self.Topology.D['mol2_bonds']=self.Coordinates.mol2_bonds.copy()
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
        # logging.debug(f'Swapping names of atoms {ai}({iname}) and {aj}({jname})')
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

    def write_gro(self,grofilename):
        """Write a Gromacs-format coordinate file; wrapper for Coordinates.write_gro()

        :param grofilename: name of file to write
        :type grofilename: str
        """
        self.Coordinates.write_gro(grofilename)

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

    def attenuate_bond_parameters(self,bonds,i,n,minimum_distance=0.0,init_colname='initial-distance'):
        """Alter the kb and b0 parameters for new crosslink bonds according to the values prior to 
            relaxation (stored in lengths), their equilibrium values, and the ratio stage/max_stages.
            Let stage/max_stages be x, and 1/max_stages <= x <= 1.  The spring constant for each
            bond is multiplied by x and the distance is 1 xth of the way from its maximum value 
            to its equilibrium value.

        :param bonds: bonds dataframe, 'ai', 'aj', 'initial-distance'
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
            to the values prior to dragging (stored in pairdf['initial-distances']), 
            the desired lower limit of interatomic distance 'draglimit_nm', 
            and the ratio stage/max_stages.
            
        :param pairdf: pairs dataframe (['ai'],['aj'],['initial-distance'])
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

    def restore_bond_parameters(self,saved):
        """Retores saved bond parameters in df saved by overwriting

        :param saved: [ bonds ] dataframe
        :type saved: pandas.DataFrame
        """
        self.Topology.restore_bond_parameters(saved)

    def write_gro_attributes(self,attributes_list,grxfilename):
        """Writes atomic attributes to a file

        :param attributes_list: list of attributes to write
        :type attributes_list: list
        :param grxfilename: name of output file
        :type grxfilename: str
        """
        self.Coordinates.write_atomset_attributes(attributes_list,grxfilename)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
        """Read attributes from file into self.Coordinates.A

        :param grxfilename: name of input file
        :type grxfilename: str
        :param attribute_list: list of attributes to take, defaults to [] (take all)
        :type attribute_list: list, optional
        """
        self.Coordinates.read_atomset_attributes(grxfilename,attributes=attribute_list)

    def set_gro_attribute(self,attribute,srs):
        self.Coordinates.set_atomset_attribute(attribute,srs)

    def set_gro_attribute_by_attributes(self,att_name,att_value,attribute_dict):
        self.Coordinates.set_atom_attribute(att_name,att_value,attribute_dict)

    def get_gro_attribute_by_attributes(self,att_name,attribute_dict):
        return self.Coordinates.get_atom_attribute(att_name,attribute_dict)

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

    def interresidue_parters_of(self,i):
        result=[]
        bl=self.Topology.bondlist.partners_of(i)
        myresid=self.Coordinates.A.iloc[i-1]['resNum']
        for j in bl:
            theirresid=self.Coordinates.A.iloc[j-1]['resNum']
            if theirresid!=myresid:
                result.append(j)
        return result

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        return self.Coordinates.minimum_distance(other.Coordinates,self_excludes=self_excludes,other_excludes=other_excludes)

    def ring_detector(self):
        return self.Topology.ring_detector()

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

    def adjust_charges(self,netcharge=0.0,msg=''):
        self.Topology.adjust_charges(desired_charge=netcharge,msg=msg)

    def gro_DataFrame(self,name):
        if name=='atoms':
            return self.Coordinates.A
        elif name=='mol2_bonds':
            return self.Coordinates.mol2_bonds
        else:
            return None

    def overwrite_coords(self,other):
        logging.debug(f'Overwriting {other.Coordinates.A.shape[0]} coordinates')
        C=self.Coordinates.A
        C=C.set_index('globalIdx')
        # logging.debug(f'before update:\n{C.to_string()}')
        B=other.Coordinates.A
        B=B.set_index('globalIdx')
        B=B[['posX','posY','posZ']].copy()
        # logging.debug(f'new coordinates:\n{B.to_string()}')
        C.update(B)
        self.Coordinates.A=C.reset_index()
        # logging.debug(f'after update:\n{self.Coordinates.A.to_string()}')

    def set_z(self,reaction_list,moldict):
        """Sets the z attribute of each atom in self using information in the list of 
        reactions and the dictionary of molecules.  'z' is an integer showing the
        number of available interresidue bonds an atom can participate in.

        :param reaction_list: List of Reactions read in from cfg file or autogenerated via 
            symmetry operations
        :type reaction_list: list of Reactions
        :param moldict:  Dictionary of all molecular templates
        :type moldict:  Dictionary
        """
        razdict={}  # [resname][atomname]=[list of z-values detected]
        for R in reaction_list:
            # logging.debug(f'Set z: scanning reaction {R.name} for raz')
            raz=R.get_raz(moldict)  # [resname][atomname]=[list of z-values detected in this reaction]
            # logging.debug(f'-> raz {raz}')
            for rn in raz:
                if not rn in razdict:
                    razdict[rn]={}
                for an in raz[rn]:
                    if not an in razdict[rn]:
                        razdict[rn][an]=[]
                    razdict[rn][an].extend(raz[rn][an])
        for rn in razdict:
            for an in razdict[rn]:
                razdict[rn][an]=max(razdict[rn][an])
        # logging.debug(f'razdict {razdict}')
        zsrs=[]
        for i,r in self.Coordinates.A.iterrows():
            rn=r['resName']
            an=r['atomName']
            idx=r['globalIdx']
            z=0
            if rn in razdict:
                if an in razdict[rn]:
                    irb=self.interresidue_parters_of(idx)
                    z=razdict[rn][an]-len(irb)
            zsrs.append(z)
        self.set_gro_attribute('z',zsrs)
        self.set_gro_attribute('nreactions',[0]*len(zsrs))
        # for i,r in self.Coordinates.A.iterrows():
        #     z=r['z']
        #     if z>0:
        #         idx=r['globalIdx']
        #         irnum=r['resNum']
        #         n=self.Topology.bondlist.partners_of(idx)
        #         for jdx in n:
        #             jz=self.get_gro_attribute_by_attributes('z',{'globalIdx':jdx})
        #             jrnum=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':jdx})
        #             if irnum==jrnum and jz>0:
        #                 # two reactive atoms in the same residue bound to each other are
        #                 # part of a polyethylene chain
        #                 self.Topology.update_polyethylenes(idx,jdx)

    def label_ring_atoms(self):
        cycles=self.ring_detector()
        adf=self.gro_DataFrame('atoms')
        self.set_gro_attribute('cycle-idx',np.zeros(adf.shape[0]).astype(int))
        cidx=1
        for l,cl in cycles.items():
            for c in cl:
                for idx in c:
                    self.set_gro_attribute_by_attributes('cycle-idx',cidx,{'globalIdx':idx})
                cidx+=1
        # logging.debug(f'label_ring_atoms for {self.name}:\n{adf.to_string()}')


    def analyze_sea_topology(self):
        """Checks for consistency of atom type, charge, and mass for all atoms in each      
            symmetry class.  If consistency is lacking, logs a warning.
        """
        tadf=self.Topology.D['atoms']
        cadf=self.Coordinates.A
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
                    logging.warning(f'Warning: atoms in symmetry class {i} have different values of {attr}\n{sea_cls.to_string()}')

    def linkcell_initialize(self,cutoff,ncpu=1,force_repopulate=False):
        """Initialize the linkcell structure; a wrapper for Coordinates

        :param cutoff: minimum value of cell side-length
        :type cutoff: float
        :param ncpu: number of processors to use in populating linkcell structure in parallel, default 1
        :type ncpu: int
        """
        self.Coordinates.linkcell_initialize(cutoff,ncpu=ncpu,populate=True,force_repopulate=force_repopulate)

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

    def inherit_attributes_from_molecules(self,attribute_list,molecule_dict):
        self.Coordinates.inherit_attributes_from_molecules(attribute_list,molecule_dict)

    def make_ringlist(self):
        self.Coordinates.make_ringlist()

    def make_resid_graph(self,json_file=None,draw=None):
        self.Topology.make_resid_graph(json_file=json_file,draw=draw)

    def maxspan(self):
        """Returns the maxspan of the Coordinates (dimensions of orthorhombic 
            convex hull enclosing Coordinates). Just a wrapper.

        :return: array of x-span, y-span, z-span
        :rtype: numpy.ndarray
        """
        return self.Coordinates.maxspan()

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
        logging.debug(f'write_mol2, other_attributes:\n{other_attributes.to_string()}')
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
        return self.Coordinates.merge(other.Coordinates)

    def bondtest_par(self,B,pbc=[1,1,1],show_piercings=True):
        """Parallelization of bondtest

        :param B: list of bonds (generated by a split to be mapped onto a Pool)
        :type B: list of 2-tuples
        :param pbc: periodic boundary condition flags in each direction, defaults to [1,1,1]
        :type pbc: list, optional
        :param show_piercings: flag indicating you want to write gro files showing 
            pierced rings
        :type show_piercings: bool
        :return: list of booleans, each is True if that bond is permitted
        :rtype: list
        """
        L=[]
        for b in B:
            L.append(self.bondtest(b,pbc=pbc,show_piercings=show_piercings))
        return L

    def bondtest(self,b,pbc=[1,1,1],show_piercings=True):
        """Determine if bond b is to be allowed to form based on geometric and 
            topological criteria

        :param b: bond, tuple of two global atom indicies
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
        i=int(i)
        j=int(j)
        # check for short-circuits, defined as a residue attempting to bond to another
        # residue to which it was already previously bonded
        if self.shortcircuit(i,j):
            return BTRC.failed_short_circuit,0
        Ri=self.get_R(i)
        Rj=self.get_R(j)
        Rij=self.Coordinates.mic(Ri-Rj,pbc)
        # chk_rij=np.sqrt(Rij.dot(Rij))
        # assert np.isclose(rij,chk_rij,atol=1.e-3),f'{i} {j} {rij} {chk_rij}'
        # generate the nearest periodic image of Rj to Ri
        Rjp=Ri-Rij
        # return array of atom coordinates of ring pierced by this bond, if any
        C=self.Coordinates.ringpierce(Ri,Rjp,pbc)
        if type(C)==np.ndarray:  # this is a ring
            # all this generate a special output file for inspection
            if show_piercings:
                idx=[i,j]
                idx.extend(C['globalIdx'].to_list()) # list of globalIdx's for this output
                sub=self.Coordinates.subcoords(self.Coordinates.A[self.Coordinates.A['globalIdx'].isin(idx)].copy())
                sub.write_gro(f'ring-{i}-{j}'+'.gro')
                logging.debug(f'Ring pierced by bond ({i}){Ri} --- ({j}){Rj} : {rij}\n{C.to_string()}')
            return BTRC.failed_pierce_ring,0
        if self.polyethylene_cycle(i,j):
            return BTRC.failed_polyethylene_cycle,0
        logging.debug(f'passed bondtest: {i:>7d} {j:>7d} {rij:>6.3f} nm')
        return BTRC.passed,rij

    def shortcircuit(self,i,j):
        """Determine whether atoms i and j, if bonded, would produce a short circuit, 
           for now defined as an instance in which i and j belong to residues that are already 
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
        i_neighbors=self.partners_of(i)
        j_neighbors=self.partners_of(j)
        if i in j_neighbors or j in i_neighbors:
            # logging.debug(f'atoms {i} and {j} already on each other\'s list of bonded partners')
            return True

        assert not j in i_neighbors # haven't made the bond yet...
        assert not i in j_neighbors # haven't made the bond yet...

        # logging.debug(f'shortcircuit: i_idx {i} i_resnum {i_resNum} i_neighbors {i_neighbors}')
        for ix in i_neighbors:
            ix_resName,ix_resNum,ix_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':ix})
            # logging.debug(f'-> ix_idx {ix} ix_resnum {ix_resNum} ix_atomName {ix_atomName}')
            if ix_resNum==j_resNum:
                # logging.debug(f'resid {i_resNum} is already bound to an atom in {j_resNum}')
                return True
        # logging.debug(f'shortcircuit: j_idx {j} j_resNum {j_resNum} j_atomName {j_atomName} j_neighbors {j_neighbors}')
        for jx in j_neighbors:
            jx_resName,jx_resNum,jx_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':jx})
            # logging.debug(f'-> jx_idx {jx} jx_resnum {jx_resNum} jx_atomName {jx_atomName}')
            if jx_resNum==i_resNum:
                # logging.debug(f'resid {j_resNum} is already bound to an atom in {i_resNum}')
                return True

        return False

    def polyethylene_cycle(self,a,b):
        return self.Topology.polyethylene_cycle(a,b)

    def set_polyethylenes(self):
        """set_polyethylenes regenerate the polyethylene network using z/nreactions
        """
        CDF=self.Coordinates.A
        T=self.Topology
        T.polyethylenes=nx.DiGraph()
        BL=T.bondlist
        for i,r in CDF.iterrows():
            idx=r['globalIdx']
            irn=r['resNum']
            iz=r['z']
            inr=r['nreactions']
            rtv=iz>0 or inr>0
            n=BL.partners_of(idx)
            for jdx in n:
                jrn=self.get_gro_attribute_by_attributes('resNum',{'globalIdx':jdx})
                if irn!=jrn:
                    T.polyethylenes.add_edge(idx,jdx) # bonded, different residues
                else:
                    jz=self.get_gro_attribute_by_attributes('z',{'globalIdx':jdx})
                    jnr=self.get_gro_attribute_by_attributes('nreactions',{'globalIdx':jdx})
                    jtv=jz>0 or jnr>0
                    if rtv and jtv: # same residue, bonded, both reactive
                        T.polyethylenes.add_edge(idx,jdx)
        cycles=list(nx.simple_cycles(T.polyethylenes))
        clens={}
        for c in cycles:
            l=len(c)
            if not l in clens:
                clens[l]=0
            clens[l]+=1
        logging.debug(f'polyethyelene cycle lengths:counts {clens}')
                

