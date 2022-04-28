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
from HTPolyNet.linkcell import Linkcell
from HTPolyNet.ring import Ring,Segment

''' Generic dataframe functions '''
def _get_row_attribute(df,name,attributes):
    ''' returns a scalar value of attribute "name" in row
        expected to be uniquely defined by attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find row with attributes {attributes}'
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
    return df[list(l)][name].values[0]

def _get_rows_w_attribute(df,name,attributes):
    ''' returns a series of values of attribute "name" from
        all rows matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    assert len(ga)>0,f'Cannot find any rows with attributes {attributes}'
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

def _set_row_attribute(df,name,value,attributes):
    ''' set value of attribute name to value in all rows matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    exla={k:v for k,v in attributes.items() if not k in df}
    if len(exla)>0:
        logging.warning(f'Caller attempts to use unrecognized attributes to refer to row: {exla}')
    if name in df and len(ga)>0:
        c=[df[k] for k in ga]
        V=list(ga.values())
        l=[True]*df.shape[0]
        for i in range(len(c)):
            l = (l) & (c[i]==V[i])
        cidx=[c==name for c in df.columns]
        df.loc[list(l),cidx]=value

def _set_rows_attributes_from_dict(df,valdict,attributes):
    ''' set values of attributes in valdict dict of all rows
        matching attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    exla={k:v for k,v in attributes.items() if not k in df}
    if len(exla)>0:
        logging.warning(f'using unknown attributes to refer to atom: {exla}')
    if all([x in df for x in valdict]) and len(ga)>0:
        c=[df[k] for k in ga]
        V=list(ga.values())
        l=[True]*df.shape[0]
        for i in range(len(c)):
            l = (l) & (c[i]==V[i])
        for k,v in valdict.items():
            cidx=[c==k for c in df.columns]
            df.loc[list(l),cidx]=v 

_ANGSTROM_='Ångström'

class Coordinates:
    ''' A class for handling atom coordinates (and other attributes)
        Instance attributes:
        - A : a DataFrame that contains atom attributes, one atom per row.
        - metadat : a dictionary of metadata populated from a MOLECULE section of a mol2 file
        - N : number of atoms (A.shape[0])
        - mol2_bonds: a DataFrame that contains bond attributes, one bond per row; typically only
          read from a mol2 file, since gro coordinate files don't have bonds
        - mol2_bondlist: a Bondlist built from the mol2_bonds DataFrame (see bondlist.py)
        - linkcell: a Linkcell instance built from the atomic positions, box size, and desired cutoff
        - empty: a Boolean indicating whether or not this is an "empty" instance
        - box: 3x3 numpy array containing the box side vectors (for a rectilinear box, only diagonals have values)
        '''
    gro_attributes = ['resNum', 'resName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    #gro_colunits = ['*','*','*','*','nm','nm','nm','nm/ps','nm/ps','nm/ps']
    mol2_atom_attributes = ['globalIdx','atomName','posX','posY','posZ','type','resNum','resName','charge']
    #mol2_atom_colunits = ['*','*',_ANGSTROM_,_ANGSTROM_,_ANGSTROM_,'*','*']
    mol2_bond_attributes = ['bondIdx','ai','aj','type']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_attributes, [int, int, int, str])}
    #coord_aux_attributes = ['sea-idx','z','cycle-idx','linkcell-idx]

    def __init__(self,name=''):
        self.name=name
        self.metadat={}
        self.N=0
        self.A=pd.DataFrame()
        self.mol2_bonds=pd.DataFrame()
        self.mol2_bondlist=Bondlist()
        self.linkcell=Linkcell()
        self.empty=True
        self.box=np.zeros((3,3))
        
    @classmethod
    def read_gro(cls,filename=''):
        ''' read atom attributes from a Gromacs gro file '''
        inst=cls(filename)
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
                series={k:[] for k in cls.gro_attributes}
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
            ***ALL LENGTHS CONVERTED FROM ANGSTROMS TO NM'''
        inst=cls(name=filename)
        ''' Length units in MOL2 are always Ångström '''
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
            inst.A=pd.read_csv(sections['atom'],sep='\s+',names=Coordinates.mol2_atom_attributes)
            inst.A[['posX','posY','posZ']]*=[0.1,0.1,0.1]
            inst.N=inst.A.shape[0]
            inst.mol2_bonds=pd.read_csv(sections['bond'],sep='\s+',names=Coordinates.mol2_bond_attributes,dtype=Coordinates.mol2_bond_types)
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
        ''' copy the posX, posY, and posZ atom attributes, and the box size, from self.A to other.A '''
        assert self.A.shape[0]==other.A.shape[0],f'Cannot copy -- atom count mismatch {self.A.shape[0]} vs {other.A.shape[0]}'
        for c in ['posX','posY','posZ']:
            otherpos=other.A[c].copy()
            self.A[c]=otherpos
        self.box=np.copy(other.box)

    def inherit_attributes_from_molecules(self,attributes=[],molecules={}):
        ''' transfer any resname-atomname specific attributes from molecular templates to all residues in 
            self.  Attributes to transfer are in the list 'attributes' and the dictionary
            of molecular templates is in 'molecules' '''
        assert type(molecules)==dict,'Must pass a *dictionary* name:Molecule as parameter "molecules"'
        a=self.A
        logging.debug(f'{a.shape[0]} atoms inheriting values of {attributes} from molecular templates.')
        for att in attributes:
            a[att]=[0]*a.shape[0]
        for resname in a['resName'].unique():
            ''' resname is a unique residue name in the system '''
            if resname in molecules:
                adf=molecules[resname].Coords.A
                for i,ma in adf.iterrows():
                    atomName=ma['atomName']
                    z=ma[attributes].to_dict()
                    # logging.debug(f'{resname} {atomName} {z}')
                    _set_rows_attributes_from_dict(a,z,{'atomName':atomName,'resName':resname})
            else:
                logging.warning(f'Resname {resname} not found in molecular templates')

    def rings(self): # an iterator
        a=self.A
        for resid in a['resNum'].unique():
            mr=a[(a['resNum']==resid)&(a['cycle-idx']>0)]
            if not mr.empty:
                for ri in mr['cycle-idx'].unique():
                    R=mr[mr['cycle-idx']==ri][['posX','posY','posZ']].values
                    yield R

    def unwrap(self,R):
        for c in range(3):
            hbx=self.box[c][c]/2
            if R[c]<-hbx:
                R[c]+=self.box[c][c]
            elif R[c]>hbx:
                R[c]-=self.box[c][c]
        return R

    def pierces(self,i,j,C):
        # this filter is only called for atoms that are within 
        # cutoff distance of each other
        ri=self.A.iloc[i-1][['posX','posY','posZ']].values
        rj=self.A.iloc[j-1][['posX','posY','posZ']].values
        o=0.5*(ri+rj)
        ris=self.unwrap(ri-o)
        rjs=self.unwrap(rj-o)
        CS=[]
        for c in C:
            CS.append(self.unwrap(c-o))
        # we will shift all coords to this origin and then 
        # unwrap if necessary
        S=Segment(np.array([ris,rjs]))
        R=Ring(np.array(CS))
        pierces,point=R.segint(S)
        return pierces

    def linkcellrings(self,ends=[]):
        assert len(ends)==2
        c=[]
        for i in ends:
            c.append(self.get_atom_attribute('linkcell-idx',{'globalIdx':i}))
        pass

    def ringpierce(self,i,j):
        if hasattr(self,'linkcell'):
            # only sample rings that lie in cells through
            # which the i-j vector passes
            for r in self.linkcellrings([i,j]):
                if self.pierces(i,j,r):
                    return True
        else:
            for r in self.rings():
                if self.pierces(i,j,r):
                    return True
        return False

    def linkcell_initialize(self,cutoff=0.0):
        logging.debug('Initializing link-cell structure')
        self.linkcell.create(cutoff,self.box)
        self.linkcell.populate(self)

    def linkcelltest(self,i,j):
        ''' return True if atoms i and j are within potential interaction
            range based on current link-cell structure '''
        ci=self.get_atom_attribute('linkcell-idx',{'globalIdx':i})
        cj=self.get_atom_attribute('linkcell-idx',{'globalIdx':j})
        if self.linkcell.are_cellidx_neighbors(ci,cj):
            return True
        return False

    def bondtest(self,b,radius):
        i,j=b
        if not self.linkcelltest(i,j):
            return False
        rij=self.rij(i,j)
        # logging.debug(f'pbond {b} rij {rij:.3f}')
        if rij>radius:
            return False
        # if self.ringpierce(i,j):
        #     return False
        return True


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
        return np.sqrt(rij.dot(rij))

    def calc_distance_matrix(self):
        M=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                d=self.rij(i,j)
                M[i][j]=M[j][i]=d
        self.distance_matrix=M

    def merge(self,other):
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
        logging.debug(f'Atomset attributes read from {filename}; new Coords\n'+self.A.to_string())

    def set_atomset_attribute(self,attribute='',srs=[]):
        if attribute!='':
           self.A[attribute]=srs

    def atomcount(self):
        return self.N

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        ''' computes the minimum distance between two collections of atoms '''
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
        ''' premultiplies position of each atom by rotation matrix R '''
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            newri=np.matmul(R,ri)
            self.A.loc[i,'posX':'posZ']=newri

    def translate(self,L):
        ''' translates all atom positions by L '''
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            self.A.loc[i,'posX':'posZ']=srow.values+L

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
        return _get_row_attribute(df,'globalIdx',attributes)
    
    def get_R(self,idx):
        df=self.A
        return _get_row_attribute(df,['posX','posY','posZ'],{'globalIdx':idx})
    
    def get_atom_attribute(self,name,attributes):
        df=self.A
        return _get_row_attribute(df,name,attributes)

    def get_atoms_w_attribute(self,name,attributes):
        df=self.A
        return _get_rows_w_attribute(df,name,attributes)

    def set_atom_attribute(self,name,value,attributes):
        df=self.A
        _set_row_attribute(df,name,value,attributes)

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
        logging.debug(f'Coordinates:delete_atoms {idx}')
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
        ''' delete appropriate bonds '''
        if not self.mol2_bonds.empty:
            d=self.mol2_bonds
            indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
            indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
            self.mol2_bonds=d.take(list(indexes_to_keep)).reset_index(drop=True)
            if reindex:
                d=self.mol2_bonds
                d.ai=d.ai.map(mapper)
                d.aj=d.aj.map(mapper)
                d.bondIdx=d.index+1
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
                    f.write(''.join([atomformatters[i](v) for i,v in enumerate(list(r[self.gro_attributes]))])+'\n')
                else:
                    f.write(''.join([atomformatters[i](v) for i,v in enumerate(list(r[self.gro_attributes[:-3]]))])+'\n')
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
        ''' write a mol2-format file from coordinates, and optionally, a bonds DataFrame
            provided externally and passed in as "bondsDF" (typically this would be
            from a Topology instance). '''
        if bondsDF.empty and self.mol2_bonds.empty:
            logging.warning(f'Cannot write any bonds to MOL2 file {filename}')
        for i in self.mol2_atom_attributes:
            if not i in self.A.columns:
                logging.warning(f'No attribute "{i}" found.')
                self.A[i]=[pd.NA]*self.A.shape[0]
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
                sdf=self.A.copy()
                # remember to convert to Angstroms
                pos=(sdf.loc[:,['posX','posY','posZ']]-com)*10.0
                sdf.loc[:,['posX','posY','posZ']]=pos
                f.write(sdf.to_string(columns=self.mol2_atom_attributes,header=False,index=False,formatters=atomformatters))
                f.write('\n')
                f.write('@<TRIPOS>BOND\n')
                if not bondsDF.empty:
                    logging.info(f'Mol2 bonds from outside')
                    bdf=bondsDF[['bondIdx','ai','aj','type']]
                    bdf['bondIdx']=bdf['bondIdx'].astype(int)
                    bdf['ai']=bdf['ai'].astype(int)
                    bdf['aj']=bdf['aj'].astype(int)
                    f.write(bdf.to_string(columns=self.mol2_bond_attributes,header=False,index=False,formatters=bondformatters))
                elif not self.mol2_bonds.empty:
                    logging.info(f'Mol2 bonds from mol2_bonds')
                    f.write(self.mol2_bonds.to_string(columns=self.mol2_bond_attributes,header=False,index=False,formatters=bondformatters))
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