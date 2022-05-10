# class for methods that need to work with both Topology and Coordinates
import pandas as pd
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topology import Topology
import logging
import numpy as np
from enum import Enum
class BTRC(Enum):
    passed = 0
    fail_linkcell = 1
    fail_beyond_cutoff = 2
    fail_pierce_ring = 3
    fail_short_circuit = 4


class TopoCoord:
    def __init__(self,topfilename='',grofilename='',mol2filename=''):
        if grofilename!='':
            self.read_gro(grofilename)
        else:
            self.Coordinates=Coordinates()
        if topfilename!='':
            self.read_top(topfilename)
        else:
            self.Topology=Topology()
        if mol2filename!='':
            self.read_mol2(mol2filename)

    def make_bonds(self,pairs,skip_H=[]):
        idx_to_ignore=self.Coordinates.find_sacrificial_H(pairs,self.Topology,skip_pairs=skip_H)
        self.Topology.add_bonds(pairs,ignores=idx_to_ignore)
        rename=True if len(skip_H)>0 else False
        idx_to_delete=self.Coordinates.find_sacrificial_H(pairs,self.Topology,skip_pairs=skip_H,rename=rename)
        assert type(idx_to_delete)==list
        return idx_to_delete

    def delete_atoms(self,atomlist):
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper=self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper)==dict
        return idx_mapper

    def map_from_templates(self,bonds,moldict):
        atdf=self.Topology.D['atoms']
        for b in bonds:
            bb,template_name=b
            T=moldict[template_name]
            self.set_gro_attribute_by_attributes('reactantName',template_name,{'globalIdx':bb[0]})
            self.set_gro_attribute_by_attributes('reactantName',template_name,{'globalIdx':bb[1]})
            temp_atdf=T.TopoCoord.Topology.D['atoms']
            inst2temp,temp2inst=T.idx_mappers(self,bb)
            assert len(inst2temp)==len(temp2inst)
            check=True
            for k,v in inst2temp.items():
                check = check and (k == temp2inst[v])
            assert check,f'Error: bidirectional dicts are incompatible; bug'
            # logging.debug(f'map_from_templates:inst2temp {inst2temp}')
            # logging.debug(f'map_from_templates:temp2inst {temp2inst}')
            i_idx,j_idx=bb
            temp_i_idx=inst2temp[i_idx]
            temp_j_idx=inst2temp[j_idx]
            d=T.TopoCoord.Topology.D['bonds']
            d=d[((d.ai==temp_i_idx)&(d.aj==temp_j_idx))|
                ((d.ai==temp_j_idx)&(d.aj==temp_i_idx))].copy()
            assert d.shape[0]==1,f'Error: map_from_templates using {T.name} is sent inst-bond {i_idx}-{j_idx} which is claimed to map to {temp_i_idx}-{temp_j_idx}, but no such bond is found:\n{self.Topology.D["bonds"].to_string()}'
            temp_angles,temp_dihedrals,temp_pairs=T.get_angles_dihedrals((inst2temp[i_idx],inst2temp[j_idx]))
            # logging.debug(f'Mapping {temp_angles.shape[0]} angles, {temp_dihedrals.shape[0]} dihedrals, and {temp_pairs.shape[0]} pairs from template {T.name}')
            inst_angles=temp_angles.copy()
            inst_angles.ai=temp_angles.ai.map(temp2inst)
            inst_angles.aj=temp_angles.aj.map(temp2inst)
            inst_angles.ak=temp_angles.ak.map(temp2inst)
            d=self.Topology.D['angles']
            self.Topology.D['angles']=pd.concat((d,inst_angles),ignore_index=True)
            for a in ['ai','aj','ak']:
                for inst_atom,temp_atom in zip(inst_angles[a],temp_angles[a]):
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
            d=self.Topology.D['angles']
            check=True
            for a in ['ai','aj','ak']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in angles')
                raise Exception
            inst_dihedrals=temp_dihedrals.copy()
            inst_dihedrals.ai=temp_dihedrals.ai.map(temp2inst)
            inst_dihedrals.aj=temp_dihedrals.aj.map(temp2inst)
            inst_dihedrals.ak=temp_dihedrals.ak.map(temp2inst)
            inst_dihedrals.al=temp_dihedrals.al.map(temp2inst)
            d=self.Topology.D['dihedrals']
            self.Topology.D['dihedrals']=pd.concat((d,inst_dihedrals),ignore_index=True)
            for a in ['ai','aj','ak','al']:
                for inst_atom,temp_atom in zip(inst_dihedrals[a],temp_dihedrals[a]):
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

            temp_pairs.ai=temp_pairs.ai.map(temp2inst)
            
            if temp_pairs.ai.isnull().values.any():
                logging.error('NAN in pairs ai')
                raise Exception
            temp_pairs.aj=temp_pairs.aj.map(temp2inst)
            if temp_pairs.aj.isnull().values.any():
                logging.error('NAN in pairs aj')
                raise Exception
            self.Topology.D['pairs']=pd.concat((d,temp_pairs),ignore_index=True)
            d=self.Topology.D['pairs']
            check=True
            for a in ['ai','aj']:
                check=check and d[a].isnull().values.any()
            if check:
                logging.error('NAN in pairs post mapping')
                raise Exception

            # bondtree=self.Topology.bondtree_as_list((i_idx,j_idx),depth=3)
            # mappables=[]
            # for B in bondtree:
            #     for i in range(2):
            #         if not B[i] in mappables:
            #             mappables.append(B[i])
            # madf=atdf[atdf['nr'].isin(mappables)]
            # logging.debug(f'Bond {b} mappable atoms:\n{madf.to_string()}')
            # Tatdf=T.TopoCoord.Topology.D['atoms']
            # Tai=inst2temp[i_idx]
            # Taj=inst2temp[j_idx]
            # T_bondtree=T.TopoCoord.Topology.bondtree_as_list((Tai,Taj),depth=3)
            # frommables=[]
            # for B in T_bondtree:
            #     for i in range(2):
            #         if not B[i] in frommables:
            #             frommables.append(B[i])
            # T_madf=Tatdf[Tatdf['nr'].isin(frommables)]
            # logging.debug(f'Template {T.name} bond {Tai}-{Taj} mappable atoms:\n{T_madf.to_string()}')
            # logging.debug(f'Same size? {madf.shape[0]==T_madf.shape[0]}')
            # for i,r in T_madf.iterrows():
            #     an,rn,ty,ch=r['atom'],r['residue'],r['type'],r['charge']
            #     atdf.loc[(atdf['nr'].isin(mappables))&(madf['atom']==an)&(madf['residue']==rn),'type']=ty
            #     atdf.loc[(atdf['nr'].isin(mappables))&(madf['atom']==an)&(madf['residue']==rn),'charge']=ch
            # logging.debug(f'Bond {b} mapped atoms:\n{atdf[atdf["nr"].isin(mappables)].to_string()}')

    def read_top(self,topfilename):
        self.Topology=Topology.read_gro(topfilename)

    def num_atoms(self):
        return self.Topology.D['atoms'].shape[0]
    
    def read_gro(self,grofilename):
        self.Coordinates=Coordinates.read_gro(grofilename)

    def read_mol2(self,mol2filename,ignore_bonds=False):
        self.Coordinates=Coordinates.read_mol2(mol2filename)
        if not ignore_bonds:
            self.Topology.D['mol2_bonds']=self.Coordinates.mol2_bonds.copy()
            if 'bonds' in self.Topology.D:
                self.Topology.bond_source_check()

    def swap_atom_names(self,ai,aj):
        T=self.Topology.D['atoms']
        C=self.Coordinates.A
        l1=T.columns=='atom'
        l2=C.columns=='atomName'
        iname=T.iloc[ai-1,l1].values[0]
        jname=T.iloc[aj-1,l1].values[0]
        logging.debug(f'Swapping names of atoms {ai}({iname}) and {aj}({jname})')
        tmpNm=T.iloc[ai-1,l1].values[0]
        T.iloc[ai-1,l1]=T.iloc[aj-1,l1]
        T.iloc[aj-1,l1]=tmpNm
        C.iloc[ai-1,l2]=C.iloc[aj-1,l2]
        C.iloc[aj-1,l2]=tmpNm

    def read_top_gro(self,topfilename,grofilename):
        self.read_top(topfilename)
        self.read_gro(grofilename)

    def write_top(self,topfilename):
        self.Topology.to_file(topfilename)

    def write_gro(self,grofilename):
        self.Coordinates.write_gro(grofilename)

    def write_top_gro(self,topfilename,grofilename):
        self.write_top(topfilename)
        self.write_gro(grofilename)

    def return_bond_lengths(self,bonds):
        return self.Coordinates.return_bond_lengths(bonds)

    def copy_bond_parameters(self,bonds):
        return self.Topology.copy_bond_parameters(bonds)

    def attenuate_bond_parameters(self,bonds,i,n,lengths):
        self.Topology.attenuate_bond_parameters(bonds,i,n,lengths)

    def copy_coords(self,other):
        self.Coordinates.copy_coords(other.Coordinates)
        self.Coordinates.box=other.Coordinates.box.copy()

    def restore_bond_parameters(self,saved):
        self.Topology.restore_bond_parameters(saved)

    def write_gro_attributes(self,attributes_list,grxfilename):
        self.Coordinates.write_atomset_attributes(attributes_list,grxfilename)

    def read_gro_attributes(self,grxfilename,attribute_list=[]):
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
        # logging.debug(f'Calling Coordinates.translate with arg {L}')
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
        razdict={}  # [resname][atomname]=[list of z-values detected]
        for R in reaction_list:
            logging.debug(f'Set z: scanning reaction {R.name} for raz')
            raz=R.get_raz(moldict)  # [resname][atomname]=[list of z-values detected in this reaction]
            logging.debug(f'-> raz {raz}')
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
        logging.debug(f'razdict {razdict}')
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

    def analyze_sea_topology(self):
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
                    logging.debug(f'Error: atoms in symmetry class {i} have different values of {attr}')

    def linkcell_initialize(self,cutoff):
        self.Coordinates.linkcell_initialize(cutoff)

    def atom_count(self):
        assert self.Coordinates.A.shape[0]==self.Topology.D['atoms'].shape[0]
        return self.Coordinates.A.shape[0]

    def total_mass(self,units='SI'):
        return self.Topology.total_mass(units=units)

    def inherit_attributes_from_molecules(self,attribute_list,molecule_dict):
        self.Coordinates.inherit_attributes_from_molecules(attribute_list,molecule_dict)

    def make_ringlist(self):
        self.Coordinates.make_ringlist()

    def make_resid_graph(self):
        self.Topology.make_resid_graph()

    def maxspan(self):
        return self.Coordinates.maxspan()

    def write_mol2(self,filename,molname=''):
        if molname=='':
            molname='This Molecule has no name'
        other_attributes=pd.DataFrame()
        other_attributes['type']=self.Topology.D['atoms']['type']
        other_attributes['charge']=self.Topology.D['atoms']['charge']
        if 'mol2_bonds' in self.Topology.D:
            self.Coordinates.write_mol2(filename,molname=molname,bondsDF=self.Topology.D['mol2_bonds'],other_attributes=other_attributes)
        else:
            self.Coordinates.write_mol2(filename,molname=molname,other_attributes=other_attributes)
    
    def merge(self,other):
        self.Topology.merge(other.Topology)
        return self.Coordinates.merge(other.Coordinates)

    def bondtest_par(self,B,radius,pbc=[1,1,1]):
        L=[]
        for b in B:
            L.append(self.bondtest(b,radius,pbc=pbc))
        return L

    def bondtest(self,b,radius,pbc=[1,1,1]):
        i,j=b
        if not self.Coordinates.linkcelltest(i,j):
            return BTRC.fail_linkcell,0
        Ri=self.get_R(i)
        Rj=self.get_R(j)
        Rij=self.Coordinates.mic(Ri-Rj,pbc)
        rij=np.sqrt(Rij.dot(Rij))
        if rij>radius:
            return BTRC.fail_beyond_cutoff,0
        Rjp=Ri-Rij # generate the nearest periodic image of Rj to Ri
        C=self.Coordinates.ringpierce(Ri,Rjp,pbc)
        if type(C)==np.ndarray:
            cidx=C[:,0].astype(int) # get globalIdx's
            idx=[i,j]
            idx.extend(cidx) # list of globalIdx's for this output
            sub=self.Coordinates.subcoords(self.Coordinates.A[self.Coordinates.A['globalIdx'].isin(idx)].copy())
            sub.write_gro(f'ring-{i}-{j}='+'-'.join([f'{x}' for x in cidx])+'.gro')
            logging.debug(f'Ring pierced by bond ({i}){Ri} --- ({j}){Rj} : {rij}')
            logging.debug('-'.join([f'{x}' for x in cidx]))
            logging.debug(f'\n+{C[:,1:]}')
            return BTRC.fail_pierce_ring,0
        # TODO: check for short-circuits or loops
        if self.shortcircuit(i,j):
            return BTRC.fail_short_circuit,0
        logging.debug(f'bondtest {b} {rij:.3f} ({radius})')
        return BTRC.passed,rij

    def shortcircuit(self,i,j):
        i_resName,i_resNum,i_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':i})
        j_resName,j_resNum,j_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':j})
        i_neighbors=self.partners_of(i)
        j_neighbors=self.partners_of(j)
        assert not j in i_neighbors # haven't made the bond yet...
        assert not i in j_neighbors # haven't made the bond yet...

        for ix in i_neighbors:
            ix_resName,ix_resNum,ix_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':ix})
            if ix_resNum==j_resNum:
                # logging.debug(f'resid {i_resNum} is already bound to an atom in {j_resNum}')
                return True
        for jx in j_neighbors:
            jx_resName,jx_resNum,jx_atomName=self.get_gro_attribute_by_attributes(['resName','resNum','atomName'],{'globalIdx':jx})
            if jx_resNum==i_resNum:
                # logging.debug(f'resid {j_resNum} is already bound to an atom in {i_resNum}')
                return True

        return False


