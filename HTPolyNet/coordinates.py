'''
coordinates.py -- simple class for handling atom coordinates from gro and mol2 files
'''

import pandas as pd
import numpy as np
from io import StringIO
from copy import deepcopy
#import hashlib
import logging

from HTPolyNet.bondlist import Bondlist

def get_atom_attribute(df,name,attributes):
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find atom with attributes {attributes}'
    if type(name)==list:
        name_in_df=all([n in df for n in name])
    else:
        name_in_df= name in df
    assert name_in_df,f'Attribute(s) {name} not found'
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
#        print(name,attributes,df[list(l)][name].values)
    return df[list(l)][name].values[0]

def get_atoms_w_attribute(df,name,attributes):
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find atom with attributes {attributes}'
    if type(name)==list:
        name_in_df=all([n in df for n in name])
    else:
        name_in_df= name in df
    assert name_in_df,f'Attribute(s) {name} not found'
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
#        print(name,attributes,df[list(l)][name].values)
    return df[list(l)][name].values

def set_atom_attribute(df,name,value,attributes):
    ga={k:v for k,v in attributes.items() if k in df}
    exla={k:v for k,v in attributes.items() if not k in df}
    if len(exla)>0:
        logging.warning(f'using unknown attributes to refer to atom: {exla}')
    if name in df and len(ga)>0:
        c=[df[k] for k in ga]
        V=list(ga.values())
        l=[True]*df.shape[0]
        for i in range(len(c)):
            l = (l) & (c[i]==V[i])
        cidx=[c==name for c in df.columns]
        df.loc[list(l),cidx]=value
        logging.debug(f'Result: {df.loc[list(l),cidx]}')

_ANGSTROM_='Ångström'

class Coordinates:
    ''' Handle coordinates from mol2 or gro files 
        A Coordinates instance has a dictionary of dataframes D with keys 'atoms' and 'bonds' 
        ** The bonds entry is only set by a mol2 file and only needed to write mol2 output ** '''
    gro_colnames = ['resNum', 'resName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    gro_colunits = ['*','*','*','*','nm','nm','nm','nm/ps','nm/ps','nm/ps']
    mol2_atom_colnames = ['globalIdx','atomName','posX','posY','posZ','type','resNum','resName','charge']
    mol2_atom_colunits = ['*','*',_ANGSTROM_,_ANGSTROM_,_ANGSTROM_,'*','*']
    mol2_bond_colnames = ['bondIdx','ai','aj','type']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_colnames, [int, int, int, str])}
    coord_aux_attributes = ['sea-idx','z']

    def __init__(self,name=''):
        self.name=name
#        self.format=''
#        self.units={}
        self.metadat={}
        self.N=0
        self.A=pd.DataFrame()
        self.mol2_bonds=pd.DataFrame()
        self.mol2_bondlist=Bondlist()
        self.empty=True
        self.box=np.zeros((3,3))
        
    @classmethod
    def read_gro(cls,filename=''):
        inst=cls(filename)
#        inst.format='gro'
#        inst.units['length']='nm'
#        inst.units['velocity']='nm/ps'
        assert filename!=''
        logging.debug(f'coordinates:read_gro {filename}')
        if filename!='':
            with open(filename,'r') as f:
                data=f.read().split('\n')
                while '' in data:
                    data.remove('')
                inst.name=data[0]
                inst.N=int(data[1])
                inst.metadat['N']=inst.N
                series={k:[] for k in cls.gro_colnames}
                for x in data[2:-1]:
                    series['resNum'].append(int(x[0:5].strip()))
                    series['resName'].append(x[5:10].strip())
                    series['atomName'].append(x[10:15].strip())
                    ''' if formatted correctly, globalIdx is row index + 1 always! '''
                    series['globalIdx'].append(int(x[15:20].strip()))
                    numbers=list(map(float,[y.strip() for y in x[20:].split()]))
                    series['posX'].append(numbers[0])
                    series['posY'].append(numbers[1])
                    series['posZ'].append(numbers[2])
                    if len(numbers)==6:
                        series['velX'].append(numbers[3])
                        series['velY'].append(numbers[4])
                        series['velZ'].append(numbers[5])
                if len(series['velX'])==0:
                    del series['velX']
                    del series['velY']
                    del series['velZ']
                assert inst.N==len(series['globalIdx']), f'Atom count mismatch inside {filename}'
                # for k,v in series.items():
                #     logging.debug(f'in coordinates.read_gro: {k} has {len(v)} items.')
                inst.A=pd.DataFrame(series)
                boxdataline=data[-1]
                n=10
                boxdata=list(map(float,[boxdataline[i:i+n].strip() for i in range(0,len(boxdataline),n)]))
                inst.box[0][0],inst.box[1][1],inst.box[2][2]=boxdata[0:3]
                if len(boxdata)==9:
                    inst.box[0][1],inst.box[0][2],inst.box[1][0],inst.box[1][2],inst.box[2][0],inst.box[2][1]=boxdata[3:]
        inst.empty=False
        return inst

    @classmethod
    def read_mol2(cls,filename=''):
        ''' Reads in a Sybyl MOL2 file into a Coordinates instance. 
            Note that this method only reads in
            MOLECULE, ATOM, and BOND sections.  
            ***ALL LENGTHS CONVERT FROM ANGSTROMS TO NM'''
        inst=cls(name=filename)
#        inst.format='mol2'
        ''' Length units in MOL2 are always Ångström '''
#        inst.units['length']=_ANGSTROM_
#        inst.units['velocity']=_ANGSTROM_+'/ps'
        if filename=='':
            return inst
        with open(filename,'r') as f:
            rawsections=f.read().split('@<TRIPOS>')[1:]
            sections={}
            for rs in rawsections:
                s=rs.split('\n')
                key=s[0].strip().lower()
                val=[a.strip() for a in s[1:] if len(a)>0]
                if key=='atom' or key=='bond':
                    val=StringIO('\n'.join(val))
                sections[key]=val
            inst.name=sections['molecule'][0]
            imetadat=list(map(int,sections['molecule'][1].strip().split()))
            inst.metadat['N']=inst.N=imetadat[0]
            inst.metadat['nBonds']=imetadat[1]
            inst.metadat['nSubs']=imetadat[2]
            inst.metadat['nFeatures']=imetadat[3]
            inst.metadat['nSets']=imetadat[4]
            inst.metadat['mol2type']=sections['molecule'][2]
            inst.metadat['mol2chargetype']=sections['molecule'][3]
            inst.A=pd.read_csv(sections['atom'],sep='\s+',names=Coordinates.mol2_atom_colnames)
            inst.A[['posX','posY','posZ']]*=[0.1,0.1,0.1]
            # print(inst.A.to_string())
            inst.N=inst.A.shape[0]
            inst.mol2_bonds=pd.read_csv(sections['bond'],sep='\s+',names=Coordinates.mol2_bond_colnames,dtype=Coordinates.mol2_bond_types)
            inst.mol2_bondlist=Bondlist.fromDataFrame(inst.mol2_bonds)
        inst.empty=False
        return inst

    def set_box(self,box):
        if box.shape==(3,1):
            for i in range(3):
                self.box[i,i]=box[i]
        elif box.shape==(3,3):
            self.box=np.copy(box)

    def copy_coords(self,other):
        assert self.A.shape[0]==other.A.shape[0],f'Cannot copy -- atom count mismatch {self.A.shape[0]} vs {other.A.shape[0]}'
        # otherfac=1.0
        # if self.units['length']=='nm' and other.units['length']==_ANGSTROM_:
        #     otherfac=0.1
        # elif self.units['length']==_ANGSTROM_ and other.units['length']=='nm':
        #     otherfac=10.0
        for c in ['posX','posY','posZ']:
            otherpos=other.A[c].copy()
            # otherpos*=otherfac
            self.A[c]=otherpos
        self.box=np.copy(other.box)

    def geometric_center(self):
        a=self.A
        return np.array([a.posX.mean(),a.posY.mean(),a.posZ.mean()])

    def rij(self,i,j,pbc=[1,1,1]):
        ''' compute distance between atoms i and j
            We assume that the DF row index is the
            globalIdx-1! '''
        if np.any(pbc) and not np.any(self.box):
            logging.warning('Interatomic distance calculation using PBC with no boxsize set.')
        ri=self.A.iloc[i-1][['posX','posY','posZ']].values
        rj=self.A.iloc[j-1][['posX','posY','posZ']].values
        rij=ri-rj
        for c in range(3):
            if pbc[c]:
                hbx=self.box[c][c]/2
                if rij[c]<-hbx:
                    rij[c]+=self.box[c][c]
                elif rij[c]>hbx:
                    rij[c]-=self.box[c][c]
        # logging.info(f'i {i} j {j} rij {rij} mag {np.sqrt(rij.dot(rij))}')
        # logging.info(f'pbc {pbc} box {self.box[0][0]} {self.box[1][1]} {self.box[2][2]}')
        return np.sqrt(rij.dot(rij))

    def calc_distance_matrix(self):
        M=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                d=self.rij(i,j)
                M[i][j]=M[j][i]=d
        self.distance_matrix=M

    def merge(self,other):
        # if self.units['length']!=other.units['length']:
        #     raise Exception('Cannot merge coordinates in different units')
        ''' get atom index, bond index, and resnum index shifts '''
        idxshift=self.A.shape[0]
        bdxshift=self.mol2_bonds.shape[0]
        rdxshift=0 if self.A.empty else self.A.iloc[-1]['resNum']
        nOtherBonds=0
        if not other.A.empty:
            oa=other.A.copy()
            ''' shift residue indices in other before merging '''
            oa['globalIdx']+=idxshift
            oa['resNum']+=rdxshift
            self.A=pd.concat((self.A,oa),ignore_index=True)
            self.N+=oa.shape[0]
        if not other.mol2_bonds.empty:
            ob=other.mol2_bonds.copy()
            ''' count number of mol2_bonds in other '''
            nOtherBonds=ob.shape[0]
            ''' shift bond indices in other '''
            ob['bondIdx']+=bdxshift
            for i in ['ai','aj']:
                ob[i]+=idxshift
            self.mol2_bonds=pd.concat((self.mol2_bonds,ob),ignore_index=True)
            self.mol2_bondlist.update(ob)
        self.metadat['N']=self.N
        if 'nBonds' in self.metadat:
            self.metadat['nBonds']+=nOtherBonds
        else:
            self.metadat['nBonds']=nOtherBonds
        
        return (idxshift,bdxshift,rdxshift)
            
    def write_atomset_attributes(self,attributes=[],filename='',formatters=[]):
        for a in attributes:
            if not a in self.A.columns:
                raise Exception(f'There is no column "{a}" in this atoms dataframe')
        if filename=='':
            raise Exception('Please provide a file name to write atom attribute data')
        with open(filename,'w') as f:
            if len(formatters)>0:
                f.write(self.A[['globalIdx']+attributes].to_string(header=True,index=False,formatters=formatters))
            else:
                f.write(self.A[['globalIdx']+attributes].to_string(header=True,index=False))

    def read_atomset_attributes(self,filename='',attributes=[]):
        if filename=='':
            raise Exception('Please provide a file name from which you want to read atom attributes')
        df=pd.read_csv(filename,sep='\s+',names=['globalIdx']+attributes)
        self.A=self.A.merge(df,how='outer',on='globalIdx')
        logging.debug('Atomset attributes read from {filename}; new Coords\n'+self.A.to_string())

    def set_atomset_attribute(self,attribute='',srs=[]):
        if attribute!='':
           self.A[attribute]=srs

    def atomcount(self):
        return self.N

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        ''' computes the minimum distance between two configurations '''
        sp=self.A[~self.A['globalIdx'].isin(self_excludes)][['posX','posY','posZ']]
        op=other.A[~other.A['globalIdx'].isin(other_excludes)][['posX','posY','posZ']]
        minD=1.e9
        for i,srow in sp.iterrows():
            ri=srow.values
            for j,orow in op.iterrows():
                rj=orow.values
                rij=ri-rj
                # print(i,j,rij)
                D=np.sqrt(np.dot(rij,rij))
                if D<minD:
                    minD=D
        return minD

    def rotate(self,R):
        ''' multiplies rotation matrix R by position of each atom '''
        # print('rotate R.shape',R.shape)
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            # print(i,ri,ri.shape)
            newri=np.matmul(R,ri)
            # print('rot ri',ri,'newri',newri)
            self.A.loc[i,'posX':'posZ']=newri
            # print(sp.loc[i,'posX':'posZ'])

    def translate(self,L):
        ''' translates all atom positions by L '''
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            newri=ri+L
            # print('tra ri',ri,'newri',newri)
            self.A.loc[i,'posX':'posZ']=newri

    def maxspan(self):
        sp=self.A[['posX','posY','posZ']]
        return np.array(
            [
                sp.posX.max()-sp.posX.min(),
                sp.posY.max()-sp.posY.min(),
                sp.posZ.max()-sp.posZ.min()
            ]
        )

    def get_idx(self,attributes):
        df=self.A
        return get_atom_attribute(df,'globalIdx',attributes)
    
    def get_R(self,idx):
        df=self.A
        return get_atom_attribute(df,['posX','posY','posZ'],{'globalIdx':idx})
    
    def get_atom_attribute(self,name,attributes):
        df=self.A
        return get_atom_attribute(df,name,attributes)

    def get_atoms_w_attribute(self,name,attributes):
        df=self.A
        return get_atoms_w_attribute(df,name,attributes)

    def set_atom_attribute(self,name,value,attributes):
        df=self.A
        set_atom_attribute(df,name,value,attributes)

    def has_atom_attributes(self,attributes):
        df=self.A
        return all([name in df for name in attributes])

    def delete_atoms(self,idx=[],reindex=True):
        '''
        Deletes atoms whose global indices appear in the list idx.
        If parameter 'reindex' is true, then the global indices 
        are recalculated so that they are sequential starting at 1 with no
        gaps, and two new columns are added to self.DF:
          - 'oldGlobalIdx' contains the global index values before 
             the deletion.
          - 'globalIdxShift' is the change from the old to the new
             global index for each atom.
        '''
        # print('delete_atoms',idx)
        adf=self.A
        indexes_to_drop=adf[adf.globalIdx.isin(idx)].index
        indexes_to_keep=set(range(adf.shape[0]))-set(indexes_to_drop)
        self.A=adf.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            adf=self.A
            oldGI=adf['globalIdx'].copy()
            adf['globalIdx']=adf.index+1
            mapper={k:v for k,v in zip(oldGI,adf['globalIdx'])}
        self.N-=len(idx)
        # print('mapper',mapper)
        ''' delete bonds '''
        # print('delete bonds containing',idx)
        if not self.mol2_bonds.empty:
            d=self.mol2_bonds
            # print('bonds before deletion:\n',d.to_string(index=False))
#            print(d.ai.isin(idx),d.aj.isin(idx))
            indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
            # print(self.D['bonds'].iloc[indexes_to_drop].to_string(index=False))
            indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
            self.mol2_bonds=d.take(list(indexes_to_keep)).reset_index(drop=True)
            # print('bonds after deletion:\n',self.D['bonds'].to_string(index=False))
            if reindex:
                d=self.mol2_bonds
                d.ai=d.ai.map(mapper)
                d.aj=d.aj.map(mapper)
                d.bondIdx=d.index+1
                # print('bonds after reindexing:\n',self.D['bonds'].to_string(index=False))
            if 'nBonds' in self.metadat:
                self.metadat['nBonds']=len(self.mol2_bonds)
            self.bondlist=Bondlist.fromDataFrame(self.mol2_bonds)

    def write_gro(self,filename=''):
        ''' write coordinates in Gromacs format '''
        if filename=='':
            raise Exception('write_gro needs a filename')
        has_vel='velX' in self.A.columns
        with open(filename,'w') as f:
            f.write(self.name+'\n')
            f.write(f'{self.N:>5d}\n')
            # C-format: “%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f”
            atomformatters = [
                lambda x: f'{x:>5d}',
                lambda x: f'{x:<5s}',
                lambda x: f'{x:>5s}',
                lambda x: f'{x:5d}']+[lambda x: f'{x:8.3f}']*3 + [lambda x: f'{x:8.4f}']*3
            # unfortunately, DataFrame.to_string() can't write fields with zero whitespace
            for i,r in self.A.iterrows():
                if has_vel:
                    f.write(''.join([atomformatters[i](v) for i,v in enumerate(list(r[self.gro_colnames]))])+'\n')
                else:
                    f.write(''.join([atomformatters[i](v) for i,v in enumerate(list(r[self.gro_colnames[:-3]]))])+'\n')
            if not np.any(self.box):
                logging.warning('Writing Gromacs coordinates file but boxsize is not set.')
            f.write(f'{self.box[0][0]:10.5f}{self.box[1][1]:10.5f}{self.box[2][2]:10.5f}')
            # output off-diagonals only if at least one of them is non-zero
            x,y=self.box.nonzero()
            if not all(x==y):
                f.write(f'{self.box[0][1]:10.5f}{self.box[0][2]:10.5f}')
                f.write(f'{self.box[1][0]:10.5f}{self.box[1][2]:10.5f}')
                f.write(f'{self.box[2][0]:10.5f}{self.box[2][1]:10.5f}')
            f.write('\n')

    def write_mol2(self,filename='',bondsDF=pd.DataFrame(),molname=''):
        if bondsDF.empty and self.mol2_bonds.empty:
            logging.warning(f'Cannot write any bonds to MOL2 file {filename}')
        for i in self.mol2_atom_colnames:
            if not i in self.A.columns:
                logging.warning(f'No attribute "{i}" found.')
                self.A[i]=[pd.NA]*self.A.shape[0]
        # assert self.format=='mol2','This config instance is not in mol2 format'
        # posscale=1.0
        # if self.units['length']=='nm':
        #     posscale=10.0
        com=self.geometric_center()
        if filename!='':
            atomformatters = [
                lambda x: f'{x:>7d}',
                lambda x: f'{x:<8s}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:<5s}',
                lambda x: f'{x:>3d}',
                lambda x: f' {x:<7s}',
                lambda x: f'{x:>9.4f}'
            ]
            bondformatters = [
                lambda x: f'{x:>6d}',
                lambda x: f'{x:>5d}',
                lambda x: f'{x:>5d}',
                lambda x: f'{str(x):>4s}'
            ]
            substructureformatters = [
                lambda x: f'{x:>6d}',
                lambda x: f'{x:<7s}',
                lambda x: f'{x:>6d}',
                lambda x: f'{x:<7s}'
            ]
            with open(filename,'w') as f:
                f.write('@<TRIPOS>MOLECULE\n')
                if molname=='':
                    f.write(f'{self.name}\n')
                else:
                    f.write(f'{molname}\n')
                N=self.N
                # Infer the residue names and resids from the atom records
                rdf=self.A[['resNum','resName']].copy().drop_duplicates()
                rdf['rootatom']=[1]*len(rdf)
                rdf['residue']=['RESIDUE']*len(rdf)
                nBonds=self.metadat.get('nBonds',0)
                if not nBonds:
                    nBonds=bondsDF.shape[0]
                nSubs=len(rdf)
                nFeatures=self.metadat.get('nFeatures',0)
                nSets=self.metadat.get('nSets',0)
                f.write('{:>6d}{:>6d}{:>3d}{:>3d}{:>3d}\n'.format(N,nBonds,nSubs, nFeatures,nSets))
                f.write(f"{self.metadat.get('mol2type','SMALL')}\n")
                f.write(f"{self.metadat.get('mol2chargetype','GASTEIGER')}\n")
                f.write('\n')
                f.write('@<TRIPOS>ATOM\n')
                # if posscale!=1.0:
                sdf=self.A.copy()
                pos=(sdf.loc[:,['posX','posY','posZ']]-com)*10.0
                sdf.loc[:,['posX','posY','posZ']]=pos
                f.write(sdf.to_string(columns=self.mol2_atom_colnames,header=False,index=False,formatters=atomformatters))
                # else:
                # f.write(self.D['atoms'].to_string(columns=self.mol2_atom_colnames,header=False,index=False,formatters=atomformatters))
                f.write('\n')
                f.write('@<TRIPOS>BOND\n')
                if not bondsDF.empty:
                    logging.info(f'Mol2 bonds from outside')
                    bdf=bondsDF[['bondIdx','ai','aj','type']]
                    bdf['bondIdx']=bdf['bondIdx'].astype(int)
                    bdf['ai']=bdf['ai'].astype(int)
                    bdf['aj']=bdf['aj'].astype(int)
                    f.write(bdf.to_string(columns=self.mol2_bond_colnames,header=False,index=False,formatters=bondformatters))
                elif not self.mol2_bonds.empty:
                    logging.info(f'Mol2 bonds from mol2_bonds')
                    f.write(self.mol2_bonds.to_string(columns=self.mol2_bond_colnames,header=False,index=False,formatters=bondformatters))
                f.write('\n')

                ''' write substructure section '''
                f.write('@<TRIPOS>SUBSTRUCTURE\n')
                f.write(rdf.to_string(header=False,index=False,formatters=substructureformatters))
"""
MOL2 Format

SUBSTRUCTURE

@<TRIPOS>SUBSTRUCTURE
Each data record associated with this RTI consists of a single data line. The data
line contains the substructure ID, name, root atom of the substructure,
substructure type, dictionary type, chain type, subtype, number of inter
substructure bonds, SYBYL status bits, and user defined comment.
Format:
subst_id subst_name root_atom [subst_type [dict_type
[chain [sub_type [inter_bonds [status
[comment]]]]]]]
• subst_id (integer) = the ID number of the substructure. This is provided
for reference only and is not used by the MOL2 command when reading
the file.
• subst_name (string) = the name of the substructure.
• root_atom (integer) = the ID number of the substructures root atom.
• subst_type - string) = the substructure type: temp, perm, residue, group
or domain.
• dict_type (integer) = the type of dictionary associated with the
substructure.
• chain (string) = the chain to which the substructure belongs (ð 4 chars).
• sub_type (string) = the subtype of the chain.
• inter_bonds (integer) = the number of inter substructure bonds.
• status (string) = the internal SYBYL status bits. These should never be
set by the user. Valid status bit values are LEAF, ROOT, TYPECOL, DICT,
BACKWARD and BLOCK.
• comment (remaining strings on data line) = the comment for the
substructure.
Example:
1 ALA1 1 RESIDUE 1 A ALA 1 ROOT|DICT Comment here
The substructure has 1 as ID, ALA1 as name and atom 1 as root atom. It is
of type RESIDUE and the associated dictionary type is 1 (protein). It is part
of the A chain in the molecule and it is an ALAnine. There is only one inter
substructure bonds. The SYBYL status bits indicate it is the ROOT
substructure of the chain and it came from a dictionary. The comment reads
“Comment here”.
1 ALA1 1
Minimal representation of a substructure.


"""