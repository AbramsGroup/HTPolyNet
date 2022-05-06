# class for methods that need to work with both Topology and Coordinates
import pandas as pd
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topology import Topology
import logging

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

    def make_bonds(self,pairs):
        idx_to_ignore=self.Coordinates.find_sacrificial_H(pairs,self.Topology)
        self.Topology.add_bonds(pairs,ignores=idx_to_ignore)
        idx_to_delete=self.Coordinates.find_sacrificial_H(pairs,self.Topology,rename=True)
        assert type(idx_to_delete)==list
        return idx_to_delete

    def delete_atoms(self,atomlist):
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper=self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper)==dict
        return idx_mapper

    def map_atomtypes_and_charges_from_templates(self,bonds):
        atdf=self.Topology.D['atoms']
        # TODO: map angles, dihedrals, and impropers!!
        for b in bonds:
            bb,template_name=b
            ai,aj=bb
            irn=self.Coordinates.get_atom_attribute('resName',{'globalIdx':ai})
            jrn=self.Coordinates.get_atom_attribute('resName',{'globalIdx':aj})
            ian=self.Coordinates.get_atom_attribute('atomName',{'globalIdx':ai})
            jan=self.Coordinates.get_atom_attribute('atomName',{'globalIdx':aj})
            bondtree=self.Topology.bondtree_as_list((ai,aj),depth=4)
            mappables=[]
            for B in bondtree:
                for i in range(2):
                    if not B[i] in mappables:
                        mappables.append(B[i])
            madf=atdf[atdf['nr'].isin(mappables)]
            # logging.debug(f'Bond {b} mappable atoms:\n{madf.to_string()}')
            T=self.molecules[template_name]
            Tatdf=T.Topology.D['atoms']
            Tair=T.Coords.get_atoms_w_attribute('resNum',{'resName':irn,'atomName':ian})
            Tairx=max(Tair)
            Tajr=T.Coords.get_atoms_w_attribute('resNum',{'resName':jrn,'atomName':jan})
            Tajrx=max(Tajr)
            Tai=T.Coords.get_atom_attribute('globalIdx',{'resNum':Tairx,'resName':irn,'atomName':ian})
            Taj=T.Coords.get_atom_attribute('globalIdx',{'resNum':Tajrx,'resName':jrn,'atomName':jan})
            T_bondtree=T.Topology.bondtree_as_list((Tai,Taj),depth=4)
            frommables=[]
            for B in T_bondtree:
                for i in range(2):
                    if not B[i] in frommables:
                        frommables.append(B[i])
            T_madf=Tatdf[Tatdf['nr'].isin(frommables)]
            # logging.debug(f'Template {T.name} bond {Tai}-{Taj} mappable atoms:\n{T_madf.to_string()}')
            # logging.debug(f'Same size? {madf.shape[0]==T_madf.shape[0]}')
            for i,r in T_madf.iterrows():
                an,rn,ty,ch=r['atom'],r['residue'],r['type'],r['charge']
                atdf.loc[(atdf['nr'].isin(mappables))&(madf['atom']==an)&(madf['residue']==rn),'type']=ty
                atdf.loc[(atdf['nr'].isin(mappables))&(madf['atom']==an)&(madf['residue']==rn),'charge']=ch
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