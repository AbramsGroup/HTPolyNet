'''
coordinates.py -- simple class for handling atom coordinates from gro and mol2 files
'''

from dis import dis
import pandas as pd
import numpy as np
from io import StringIO
from copy import deepcopy
import os
#import hashlib
import logging


from HTPolyNet.bondlist import Bondlist
from HTPolyNet.linkcell import Linkcell
from HTPolyNet.ring import Ring,Segment

# def my_check(A,item):
#     idx=item[0].astype(int)
#     for a in A:
#         if all(a[0].astype(int)==idx):
#             return True
#     return False

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
    # logging.debug(f'_get_row_attribute {name} {attributes} -> {df[list(l)][name].values}')
    return df[list(l)][name].values[0]

def _get_row_as_string(df,attributes):
    ''' returns a scalar value of attribute "name" in row
        expected to be uniquely defined by attributes dict '''
    ga={k:v for k,v in attributes.items() if k in df}
    c=[df[k] for k in ga]
    V=list(ga.values())
    l=[True]*df.shape[0]
    for i in range(len(c)):
        l = (l) & (c[i]==V[i])
    return df[list(l)].to_string()

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
        self.ringlist=[]
        
    @classmethod
    def read_gro(cls,filename=''):
        ''' read atom attributes from a Gromacs gro file '''
        inst=cls(filename)
        assert filename!=''
        # logging.debug(f'coordinates:read_gro {filename}')
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
                    # split won't work since sometimes there might be no spaces
                    # "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
                    #numbers=list(map(float,[y.strip() for y in x[20:].split()]))
                    numbers=list(map(float,[x[20+8*i:20+8*(i+1)] for i in range(0,3)]))
                    if len(x)>44:
                        numbers.extend(list(map(float,[x[44+8*i:44+8*(i+1)] for i in range(0,3)])))
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
        inst.wrap_coords()
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

    def subcoords(self,sub_adf):
        newC=Coordinates()
        newC.set_box(self.box)
        newC.A=sub_adf
        newC.N=sub_adf.shape[0]
        return newC

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
                adf=molecules[resname].TopoCoord.Coordinates.A
                for i,ma in adf.iterrows():
                    atomName=ma['atomName']
                    z=ma[attributes].to_dict()
                    _set_rows_attributes_from_dict(a,z,{'atomName':atomName,'resName':resname})
            else:
                logging.warning(f'Resname {resname} not found in molecular templates')

    def make_ringlist(self):
        self.ringlist=list(self.rings())

    def reindex_ringlist(self,idx_mapper):
        pass

    def rings(self): # an iterator over all rings
        a=self.A
        for resid in a['resNum'].unique():
            mr=a[(a['resNum']==resid)&(a['cycle-idx']>0)]
            if not mr.empty:
                for ri in mr['cycle-idx'].unique():
                    R=mr[mr['cycle-idx']==ri][['globalIdx','posX','posY','posZ']].values
                    # TODO: cast R as a pandas DataFrame to preserve int type of globalIdx
                    # logging.debug(f'visiting a ring ({resid}:{ri}) of length {R.shape[0]}')
                    yield R

    def unwrap(self,P,O,pbc):
        ''' shift all points in P so that all are CPI* to point O
            *CPI=closest periodic image (unwrapped) '''
        sP=[]
        for p in P:
            ROp=self.mic(O-p,pbc)
            pp=O-ROp
            sP.append(pp)
        return np.array(sP)

    def pierces(self,Ri,Rj,iC,pbc):
        # this filter is only called for atoms that are within 
        # cutoff distance of each other
        C=iC[:,1:]
        sC=self.unwrap(C,Ri,pbc)
        S=Segment(np.array([Ri,Rj]))
        R=Ring(sC)
        R.analyze()
        pierces,point=R.segint(S)
        return pierces

    def linkcellrings(self,Ri,Rj,discretization=0.2):
        ''' 
        Ri and Rj are locations such that Rj is a CPI to Ri
        '''
        nearby_rings=np.array([])
        Rij=Ri-Rj
        rij=np.sqrt(Rij.dot(Rij))
        nip=int(rij/discretization)
        if nip==1:
            nip=2
        collisions=0
        total_rings=0
        for C in self.ringlist:
            total_rings+=1
            # logging.debug(f'ring\n{C}')
            lcids=[]
            for ci in C:
                idx=int(ci[0])
                # logging.debug(f'asking for linkcell-idx of atom {idx}')
                try:
                    rci=self.get_atom_attribute('linkcell-idx',{'globalIdx':idx})      
                except:
                    logging.debug(f'asking for linkcell-idx of atom {idx} failed!!')
                    logging.debug(f'{self.spew_atom({"globalIdx":idx})}')
                    logging.debug(f'{len(self.ringlist)} rings; ring: {C}\n')
                    raise Exception('why?')
                lcids.append(rci)
            for p in np.linspace(Ri,Rj,nip):  # make a series of points along the bond
                cpi=self.linkcell.ldx_of_cellndx(self.linkcell.cellndx_of_point(self.wrap_point(p)))
                # logging.debug(f'intermediate point {p} in cell {cpi}...')
                nears=[]
                for rci in lcids:
                    nears.append(self.linkcell.are_ldx_neighbors(cpi,rci))
                    # logging.debug(f'ringc {rci} acpi {cpi} neighbors {self.linkcell.are_ldx_neighbors(cpi,rci)}')
                # logging.debug(f'any(nears) {any(nears)}')
                if any(nears):
                    # logging.debug(f'adding C(shape={C.shape}) to nearby_rings:\n{C}')
                    # logging.debug(f'nearby rings {nearby_rings.shape}\n{nearby_rings}')
                    if nearby_rings.size==0:
                        nearby_rings=np.array([C])
                    else:
                        # logging.debug(f'{np.any(C==nearby_rings,axis=0)} {np.all(np.any(C==nearby_rings,axis=0))}')
                        #is_in_list=my_check(nearby_rings,C)
                        is_in_list=np.all(np.any(C==nearby_rings,axis=0))
                        # logging.debug(f'not C in nearby_rings {not is_in_list}')
                        if not is_in_list:
                            nearby_rings=np.append(nearby_rings,np.array([C]),axis=0)
                        else:
                            collisions+=1
                    # logging.debug(f'after: nearby rings {nearby_rings.shape}\n{nearby_rings}')
        # logging.debug(f'linkcellrings(): {nearby_rings.shape[0]}/{total_rings} rings to be tested.')
        # logging.debug(f'Discretization of {discretization} of bond length {rij:.3f}')
        # logging.debug(f'into {nip} points resulted in {collisions} overcounts.')
        return nearby_rings

    def ringpierce(self,Ri,Rj,pbc):
        for C in self.linkcellrings(Ri,Rj):
            if self.pierces(Ri,Rj,C,pbc):
                # logging.debug(f'\n{C}')
                return C
        return False

    def linkcell_initialize(self,cutoff=0.0,ncpu=1,populate=True,force_repopulate=False):
        logging.debug('Initializing link-cell structure')
        self.linkcell.create(cutoff,self.box)
        if populate:
            if os.path.exists('linkcell.grx') and not force_repopulate:
                logging.debug(f'Found linkcell.grx; no need to populate.')
                self.read_atomset_attributes('linkcell.grx')
            else:
                self.linkcell.populate(self,ncpu=ncpu)
                self.write_atomset_attributes(['linkcell-idx'],'linkcell.grx')

    def linkcelltest(self,i,j):
        ''' return True if atoms i and j are within potential interaction
            range based on current link-cell structure '''
        ci=self.get_atom_attribute('linkcell-idx',{'globalIdx':i})
        cj=self.get_atom_attribute('linkcell-idx',{'globalIdx':j})
        if ci==cj:
            return True
        if self.linkcell.are_ldx_neighbors(ci,cj):
            return True
        return False

    def geometric_center(self):
        a=self.A
        return np.array([a.posX.mean(),a.posY.mean(),a.posZ.mean()])

    def rij(self,i,j,pbc=[1,1,1]):
        ''' compute distance between atoms i and j '''
        if np.any(pbc) and not np.any(self.box):
            logging.warning('Interatomic distance calculation using PBC with no boxsize set.')
        ri=self.get_R(i)
        rj=self.get_R(j)
        Rij=self.mic(ri-rj,pbc)
        return np.sqrt(Rij.dot(Rij))

    def mic(self,r,pbc):
        for c in range(3):
            if pbc[c]:
                hbx=self.box[c][c]/2
                if r[c]<-hbx:
                    r[c]+=self.box[c][c]
                elif r[c]>hbx:
                    r[c]-=self.box[c][c]
        return r

    _nwrap=0
    def wrap_point(self,ri):
        R=ri.copy()
        for i in range(3):
            if R[i]<0 or R[i]>=self.box[i][i]:
                self._nwrap+=1
            while R[i]<0:
                R[i]+=self.box[i][i]
            while R[i]>=self.box[i][i]:
                R[i]-=self.box[i][i]
        return R
        # bmults=np.ones(3).astype(int)-np.ceil(np.multiply(ri,np.reciprocal(self.box.diagonal()))).astype(int)
        # if any(bmults):
        #     self._nwrap+=1
        #     return ri+np.multiply(self.box.diagonal(),bmults)
        # return ri

    def wrap_coords(self):
        self._nwrap=0
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            self.A.loc[i,'posX':'posZ']=self.wrap_point(srow.values)
        logging.debug(f'Wrapped {self._nwrap}/{self.A.shape[0]} coordinates.')

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
        # if no particular attributes are asked for, read them all in
        if len(attributes)==0:
            df=pd.read_csv(filename,sep='\s+',header=0)
            assert 'globalIdx' in df,f'Error: {filename} does not have a \'globalIdx\' column'
        else:
            df=pd.read_csv(filename,sep='\s+',names=['globalIdx']+attributes,header=0)
        self.A=self.A.merge(df,how='outer',on='globalIdx')
        # logging.debug(f'Atomset attributes read from {filename}; new Coords\n'+self.A.to_string())

    def set_atomset_attribute(self,attribute='',srs=[]):
        if attribute!='':
           self.A[attribute]=srs

    def atomcount(self):
        return self.N

    def decrement_z(self,pairs):
        for b in pairs:
            ai,aj=b
            # ain=self.get_atom_attribute('atomName',{'globalIdx':ai})
            # ajn=self.get_atom_attribute('atomName',{'globalIdx':aj})
            iz=self.get_atom_attribute('z',{'globalIdx':ai})-1
            assert iz>=0,f'Error: decrementing z of atom {ai} gives erroneous z {iz}'
            jz=self.get_atom_attribute('z',{'globalIdx':aj})-1
            assert jz>=0,f'Error: decrementing z of atom {aj} gives erroneous z {jz}'
            # logging.debug(f'Setting z of {ain}-{ai} to {iz}')
            # logging.debug(f'Setting z of {ajn}-{aj} to {jz}')
            self.set_atom_attribute('z',iz,{'globalIdx':ai})
            self.set_atom_attribute('z',jz,{'globalIdx':aj})

    def show_z_report(self):
        zhists={}
        for i,r in self.A.iterrows():
            n=r['atomName']
            nn=r['resName']
            k=f'{nn}:{n}'
            z=r['z']
            if not k in zhists:
                zhists[k]=np.zeros(4).astype(int)
            zhists[k][z]+=1

        for n in zhists:
            if any([zhists[n][i]>0 for i in range(1,4)]):
                logging.debug(f'Z-hist for {n} atoms:')
                for i in range(4):
                    logging.debug(f'{i:>5d} ({zhists[n][i]:>6d}): '+'*'*(zhists[n][i]//10))

    def return_bond_lengths(self,bdf):
        lengths=[]
        for i,b in bdf.iterrows():
            lengths.append(self.rij(b['ai'],b['aj']))
        return lengths

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
        # logging.debug(f'Rotating {sp.shape[0]} atom positions by\n{R}')
        # logging.debug(f'before rotation:\n{self.A.to_string()}')
        for i,srow in sp.iterrows():
            ri=srow.values
            newri=np.matmul(R,ri)
            self.A.loc[i,'posX':'posZ']=newri
        # logging.debug(f'after rotation:\n{self.A.to_string()}')

    def translate(self,L):
        ''' translates all atom positions by L '''
        # logging.debug(f'Translating by {L}')
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            # logging.debug(f'{srow.values} to {srow.values+L}')
            self.A.loc[i,'posX':'posZ']=srow.values+L

    def maxspan(self):
        """Returns dimensions of orthorhombic convex hull enclosing Coordinates

        :return: array of x-span, y-span, z-span
        :rtype: numpy.ndarray
        """
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
        # assert name in self.A.columns
        if type(name)==list:
            assert all([i in self.A.columns for i in name])
        else:
            assert name in self.A.columns
        return _get_row_attribute(df,name,attributes)
    
    def spew_atom(self,attributes):
        df=self.A
        return _get_row_as_string(df,attributes)

    def get_atoms_w_attribute(self,name,attributes):
        df=self.A
        return _get_rows_w_attribute(df,name,attributes)

    def set_atom_attribute(self,name,value,attributes):
        df=self.A
        _set_row_attribute(df,name,value,attributes)

    def has_atom_attributes(self,attributes):
        df=self.A
        return all([name in df for name in attributes])

    def find_sacrificial_H(self,pairs,T,rename=False,skip_pairs=[]):
        idx_to_delete=[]
        for i,b in enumerate(pairs):
            if not i in skip_pairs:
                ai,aj=b
                idx_to_delete.extend(self.sacH(ai,aj,T,rename=rename))
        return idx_to_delete

    def sacH(self,ai,aj,T,rename=False):
        ''' find the two H's closest to each other to delete '''
        bondlist=T.bondlist
        i_partners=bondlist.partners_of(ai)
        j_partners=bondlist.partners_of(aj)
        i_Hpartners={k:v for k,v in zip(i_partners,[self.A[self.A['globalIdx']==i]['atomName'].values[0] for i in i_partners]) if v.startswith('H')}
        j_Hpartners={k:v for k,v in zip(j_partners,[self.A[self.A['globalIdx']==i]['atomName'].values[0] for i in j_partners]) if v.startswith('H')}
        assert len(i_Hpartners)>0,f'Error: atom {ai} does not have a deletable H atom!'
        assert len(j_Hpartners)>0,f'Error: atom {aj} does not have a deletable H atom!'
        minHH=(1.e9,-1,-1)
        for ih in i_Hpartners:
            RiH=self.get_R(ih)
            for jh in j_Hpartners:
                RjH=self.get_R(jh)
                RijH=RiH-RjH
                rijh=np.sqrt(RijH.dot(RijH))
                if rijh<minHH[0]:
                    minHH=(rijh,ih,jh)
        ''' rename remaining H atoms '''
        if rename:
            i_avails=list(sorted(i_Hpartners.values()))[:-1]
            j_avails=list(sorted(j_Hpartners.values()))[:-1]
            del i_Hpartners[ih]
            del j_Hpartners[jh]
            Top=T.D['atoms']
            Cor=self.A
            for h in i_Hpartners:
                i_Hpartners[h]=i_avails.pop(0)
                Top.iloc[h-1,Top.columns=='atom']=i_Hpartners[h]
                Cor.iloc[h-1,Cor.columns=='atomName']=i_Hpartners[h]
            for h in j_Hpartners:
                j_Hpartners[h]=j_avails.pop(0)
                Top.iloc[h-1,Top.columns=='atom']=j_Hpartners[h]
                Cor.iloc[h-1,Cor.columns=='atomName']=j_Hpartners[h]
        # this makes sure that it always looks like the same atom was deleted
        return [ih,jh]

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
        # logging.debug(f'Coordinates:delete_atoms {idx}')
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

    def write_mol2(self,filename='',bondsDF=pd.DataFrame(),molname='',other_attributes=pd.DataFrame()):
        ''' write a mol2-format file from coordinates, and optionally, a bonds DataFrame
            provided externally and passed in as "bondsDF" (typically this would be
            from a Topology instance). '''
        acopy=self.A.copy()
        if bondsDF.empty and self.mol2_bonds.empty:
            logging.warning(f'Cannot write any bonds to MOL2 file {filename}')
        for i in self.mol2_atom_attributes:
            # logging.debug(f'checking {i}')
            if not i in self.A.columns:
                if not i in other_attributes.columns:
                    logging.debug(f'No attribute "{i}" found.')
                    acopy[i]=[pd.NA]*self.A.shape[0]
                else:
                    acopy[i]=other_attributes[i]
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
                N=acopy.shape[0] #self.N
                # Infer the residue names and resids from the atom records
                rdf=acopy[['resNum','resName']].copy().drop_duplicates()
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
                # remember to convert to Angstroms
                pos=(acopy.loc[:,['posX','posY','posZ']]-com)*10.0
                acopy.loc[:,['posX','posY','posZ']]=pos
                f.write(acopy.to_string(columns=self.mol2_atom_attributes,header=False,index=False,formatters=atomformatters))
                f.write('\n')
                f.write('@<TRIPOS>BOND\n')
                if not bondsDF.empty:
                    # logging.info(f'Mol2 bonds from outside')
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