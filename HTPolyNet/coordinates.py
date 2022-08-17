"""

.. module:: coordinates
   :synopsis: Class for managing gromacs .gro file data
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

import pandas as pd
import numpy as np
from io import StringIO
import os
import logging
from itertools import product

from HTPolyNet.bondlist import Bondlist
from HTPolyNet.linkcell import Linkcell
from HTPolyNet.ring import Ring,Segment

logger=logging.getLogger(__name__)

# def my_check(A,item):
#     idx=item[0].astype(int)
#     for a in A:
#         if all(a[0].astype(int)==idx):
#             return True
#     return False

GRX_ATTRIBUTES=['z','nreactions','reactantName','sea_idx','cycle','cycle_idx','chain','chain_idx','molecule','molecule_name']
GRX_GLOBALLY_UNIQUE=[False,False,False,True,False,True,False,True,True],
GRX_UNSET_DEFAULTS=[0,0,'UNSET',-1,-1,-1,-1,-1,'UNSET']

def _dfrotate(df:pd.DataFrame,R):
    for i,srow in df.iterrows():
        ri=srow[['posX','posY','posZ']].values
        newri=np.matmul(R,ri)
        df.loc[i,'posX':'posZ']=newri

def _get_row_attribute(df:pd.DataFrame,name,attributes):
    """_get_row_attribute returns a scalar value of attribute "name" in row
        expected to be uniquely defined by attributes dict

    :param df: dataframe to search
    :type df: pandas.DataFrame
    :param name: name of attribute whose value you want
    :type name: str
    :param attributes: dictionary of attribute:value pairs that defines target set or row
    :type attributes: dict
    :return: value of attribute name
    :rtype: scalar
    """
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
    # logger.debug(f'_get_row_attribute {name} {attributes} -> {df[list(l)][name].values}')
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
        logger.warning(f'Caller attempts to use unrecognized attributes to refer to row: {exla}')
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
        logger.warning(f'using unknown attributes to refer to atom: {exla}')
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
    mol2_atom_attributes = ['globalIdx','atomName','posX','posY','posZ','type','resNum','resName','charge']
    mol2_bond_attributes = ['bondIdx','ai','aj','order']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_attributes, [int, int, int, str])}

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
        self.grx_attributes=GRX_ATTRIBUTES
        
    @classmethod
    def read_gro(cls,filename,wrap_coords=True):
        """read_gro Read a Gromacs gro file

        :param filename: name of gro file
        :type filename: str
        :return: a new Coordinates instance
        :rtype: Coordinates
        """
        inst=cls(filename)
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
                #     logger.debug(f'in coordinates.read_gro: {k} has {len(v)} items.')
                inst.A=pd.DataFrame(series)
                boxdataline=data[-1]
                boxdata=list(map(float,boxdataline.split()))
                # logger.debug(f'boxdata {boxdata}')
                inst.box[0][0]=boxdata[0]
                inst.box[1][1]=boxdata[1]
                inst.box[2][2]=boxdata[2]
                # logger.debug(f'box: {inst.box}')
                if len(boxdata)==9:
                    inst.box[0][1],inst.box[0][2],inst.box[1][0],inst.box[1][2],inst.box[2][0],inst.box[2][1]=boxdata[3:]
        inst.empty=False
        # logger.debug(f'{inst.checkbox()}')
        # logger.debug('Box vectors:')
        # for ln in str(inst.box).split('\n'):
        #     logger.debug(ln)
        if wrap_coords:
            inst.wrap_coords()
        return inst

    @classmethod
    def read_mol2(cls,filename):
        """read_mol2 Reads in a Sybyl MOL2 file into a Coordinates instance. 
            Note that this method only reads in
            MOLECULE, ATOM, and BOND sections. 

        :param filename: name of input mol2 file
        :type filename: str
        :return: a new Coordinates instance
        :rtype: Coordinates
        """
        '''***ALL LENGTHS CONVERTED FROM ANGSTROMS TO NM***'''
        inst=cls(name=filename)
        ''' Length units in MOL2 are always Ångström '''
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
            # sort so atom indices are increasing in each bond
            for i,r in inst.mol2_bonds.iterrows():
                ai=r['ai']
                aj=r['aj']
                if aj<ai:
                    logger.debug(f'mol2 bonds swapping {ai} and {aj}')
                    inst.mol2_bonds.iloc[i,inst.mol2_bonds.columns=='ai']=aj
                    inst.mol2_bonds.iloc[i,inst.mol2_bonds.columns=='aj']=ai
            inst.mol2_bondlist=Bondlist.fromDataFrame(inst.mol2_bonds)
        inst.empty=False
        return inst

    def set_box(self,box:np.ndarray):
        """set_box Set the box size from box

        :param box: 3x1 or 3x3 box size matrix
        :type box: numpy.ndarray
        """
        if box.shape==(3,1):
            for i in range(3):
                self.box[i,i]=box[i]
        elif box.shape==(3,3):
            self.box=np.copy(box)

    def total_volume(self,units='gromacs'):
        """Returns total volume of box.

        :param units: unit system designation; if 'SI' returns m^3, defaults to 'gromacs'
        :type units: str, optional
        :return: volume (in nm^3 if units is 'gromacs' or m^3 if units is 'SI')
        :rtype: float
        """
        nm_per_m=1.e9
        vol=np.prod(self.box.diagonal())  # nm^3
        return vol if units!='SI' else vol/(nm_per_m**3)

    def copy_coords(self,other):
        """copy_coords copy_coords copy the posX, posY, and posZ atom attributes, and the box size, 
        from other.A to self.A

        :param other: the other Coordinates instance
        :type other: Coordinates
        """
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
    
    def reconcile_subcoords(self,subc,attr):
        jdx=list(subc.A.columns).index(attr)
        for r in subc.A.itertuples(index=False):
            idx=r.globalIdx
            lc_idx=r[jdx]
            self.A.loc[idx-1,attr]=lc_idx

    def unwrap(self,P,O,pbc):
        ''' shift point P to its CPI* to point O
            *CPI=closest periodic image (unwrapped) '''
        ROP=self.mic(O-P,pbc)
        PCPI=O-ROP
        return PCPI

    def pierces(self,B:pd.DataFrame,C:pd.DataFrame,pbc):
        BC=np.array(B[['posX','posY','posZ']])
        CC=np.array(C[['posX','posY','posZ']])
        # get both atoms in bond into CPI
        # get all atoms in ring into CPI (but not necessarily wrt bond)
        CC[1:]=np.array([self.unwrap(c,CC[0],pbc) for c in CC[1:]])
        BC[0]=self.unwrap(BC[0],CC[0],pbc)
        BC[1]=self.unwrap(BC[1],BC[0],pbc)
        B[['posX','posY','posZ']]=BC
        C[['posX','posY','posZ']]=CC
        S=Segment(BC)
        R=Ring(CC)
        R.analyze()
        do_it,point=R.segint(S)
        return do_it

    def linkcell_initialize(self,cutoff=0.0,ncpu=1,populate=True,force_repopulate=False,save=True):
        logger.debug('Initializing link-cell structure')
        self.linkcell.create(cutoff,self.box)
        if populate:
            lc_file=f'linkcell-{cutoff:.2f}.grx'
            if os.path.exists(lc_file) and not force_repopulate:
                logger.debug(f'Found {lc_file}; no need to populate.')
                results=self.read_atomset_attributes(lc_file)
                logger.debug(f'Read linkcell_idx from {lc_file} {("linkcell_idx" in self.A)} {results}')
                self.linkcell.make_memberlists(self.A)
            else:
                self.set_atomset_attribute('linkcell_idx',-1*np.ones(self.A.shape[0]).astype(int))
                sc=self.subcoords(self.A[(self.A['cycle_idx']>0)|(self.A['z']>0)].copy())
                self.linkcell.populate(sc,ncpu=ncpu)
                self.reconcile_subcoords(sc,'linkcell_idx')
                if save:
                    self.write_atomset_attributes(['linkcell_idx'],lc_file)

    def linkcelltest(self,i,j):
        ''' return True if atoms i and j are within potential interaction
            range based on current link-cell structure '''
        ci=self.get_atom_attribute('linkcell_idx',{'globalIdx':i})
        cj=self.get_atom_attribute('linkcell_idx',{'globalIdx':j})
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
            logger.warning('Interatomic distance calculation using PBC with no boxsize set.')
        ri=self.get_R(i)
        rj=self.get_R(j)
        Rij=self.mic(ri-rj,pbc)
        return np.sqrt(Rij.dot(Rij))

    def mic(self,r,pbc):
        for c in range(0,3):
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

    def wrap_coords(self):
        """wrap_coords Wraps all atomic coordinates into box
        """
        assert np.any(self.box),f'Cannot wrap if boxsize is not set: {self.box}'
        self._nwrap=0
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            self.A.loc[i,'posX':'posZ']=self.wrap_point(srow.values)
        # logger.debug(f'Wrapped {self._nwrap}/{self.A.shape[0]*3} coordinates.')

    def calc_distance_matrix(self):
        M=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                d=self.rij(i,j)
                M[i][j]=M[j][i]=d
        self.distance_matrix=M

    def merge(self,other):
        """merge Merge two Coordinates instances

        :param other: the other Coordinates instance
        :type other: Coordinates
        :return: integer shifts in atom index, bond index, and residue index as a 3-tuple
        :rtype: tuple
        """
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
            
    def write_atomset_attributes(self,attributes,filename,formatters=[]):
        """write_atomset_attributes Writes atom attributes to a file

        :param attributes: List of attribute names to write
        :type attributes: list, optional
        :param filename: Name of file to write
        :type filename: str
        :param formatters: formatting methods per attribute, defaults to []
        :type formatters: list
        :raises Exception: All items in attributes must exist in the coordinates dataframe
        """
        for a in attributes:
            if not a in self.A.columns:
                raise Exception(f'There is no column "{a}" in this atoms dataframe')
        with open(filename,'w') as f:
            if len(formatters)>0:
                f.write(self.A[['globalIdx']+attributes].to_string(header=True,index=False,formatters=formatters)+'\n')
            else:
                f.write(self.A[['globalIdx']+attributes].to_string(header=True,index=False)+'\n')

    def read_atomset_attributes(self,filename,attributes=[]):
        """Reads atomic attributes from input file

        :param filename: name of file
        :type filename: str
        :param attributes: list of attributes to take, defaults to [] (take all)
        :type attributes: list, optional
        """
        assert os.path.exists(filename),f'Error: {filename} not found'
        # if no particular attributes are asked for, read them all in
        if len(attributes)==0:
            df=pd.read_csv(filename,sep='\s+',header=0)
            assert 'globalIdx' in df,f'Error: {filename} does not have a \'globalIdx\' column'
            attributes_read=list(df.columns)
            attributes_read.remove('globalIdx')
        else:
            df=pd.read_csv(filename,sep='\s+',names=['globalIdx']+attributes,header=0)
            attributes_read=attributes
            # logger.debug(f'Read from {filename}\n{df.head().to_string()}')
        # logger.debug(f'Merge:\n{self.A.head().to_string()}\nand\n{df.head().to_string()}')
        self.A=self.A.merge(df,how='outer',on='globalIdx')
        # logger.debug(f'Result:\n{self.A.head().to_string()}')
        return attributes_read

    def set_atomset_attribute(self,attribute,srs):
        """Set attribute of atoms to srs

        :param attribute: name of attribute
        :type attribute: str
        :param srs: scalar or list-like attribute values in same ordering as self.A
        :type srs: scalar or list-like
        """
        self.A[attribute]=srs

    def atomcount(self):
        return self.N

    def decrement_z(self,pairs):
        """decrement_z Decrements value of z attributes of all atoms found in pairs

        :param pairs: list of atom index pairs, interpreted as new bonds that just formed
        :type pairs: list of 2-tuples
        """
        for b in pairs:
            ai,aj=b
            # ain=self.get_atom_attribute('atomName',{'globalIdx':ai})
            # ajn=self.get_atom_attribute('atomName',{'globalIdx':aj})
            iz=self.get_atom_attribute('z',{'globalIdx':ai})-1
            assert iz>=0,f'Error: decrementing z of atom {ai} gives erroneous z {iz}'
            jz=self.get_atom_attribute('z',{'globalIdx':aj})-1
            assert jz>=0,f'Error: decrementing z of atom {aj} gives erroneous z {jz}'
            # logger.debug(f'Setting z of {ain}-{ai} to {iz}')
            # logger.debug(f'Setting z of {ajn}-{aj} to {jz}')
            self.set_atom_attribute('z',iz,{'globalIdx':ai})
            self.set_atom_attribute('z',jz,{'globalIdx':aj})
            inr=self.get_atom_attribute('nreactions',{'globalIdx':ai})+1
            jnr=self.get_atom_attribute('nreactions',{'globalIdx':aj})+1
            # logger.debug(f'Setting z of {ain}-{ai} to {iz}')
            # logger.debug(f'Setting z of {ajn}-{aj} to {jz}')
            self.set_atom_attribute('nreactions',inr,{'globalIdx':ai})
            self.set_atom_attribute('nreactions',jnr,{'globalIdx':aj})
            
    def show_z_report(self):
        zhists={}
        for r in self.A.itertuples():
            n=r.atomName
            nn=r.resName
            k=f'{nn}:{n}'
            z=r.z
            if not k in zhists:
                zhists[k]=np.zeros(4).astype(int)
            zhists[k][z]+=1

        for n in zhists:
            if any([zhists[n][i]>0 for i in range(1,4)]):
                logger.debug(f'Z-hist for {n} atoms:')
                for i in range(4):
                    logger.debug(f'{i:>5d} ({zhists[n][i]:>6d}): '+'*'*(zhists[n][i]//10))

    def return_bond_lengths(self,bdf):
        lengths=[]
        for b in bdf.itertuples():
            lengths.append(self.rij(b.ai,b.aj))
        return lengths

    def add_length_attribute(self,bdf,attr_name='length'):
        lengths=[]
        for b in bdf.itertuples():
            lengths.append(self.rij(b.ai,b.aj))
        bdf[attr_name]=lengths

    # def return_pair_lengths(self,pdf):
    #     lengths=[]
    #     for i,b in pdf.iterrows():
    #         lengths.append(self.rij(b['ai'],b['aj']))
    #     return lengths

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        """minimum_distance Computes and returns distance of closest approach between two sets of atoms

        :param other: other Coordinates instance
        :type other: Coordinates
        :param self_excludes: list of atom indexes in self to NOT consider, defaults to []
        :type self_excludes: list, optional
        :param other_excludes: list of atom indexes in other to NOT consider, defaults to []
        :type other_excludes: list, optional
        :return: distance of closest approach: i.e., the distance between the two atoms, one from self and one from other, that are closest together
        :rtype: float
        """
        sp=self.A[~self.A['globalIdx'].isin(self_excludes)][['posX','posY','posZ']]
        op=other.A[~other.A['globalIdx'].isin(other_excludes)][['posX','posY','posZ']]
        minD=1.e9
        for i,srow in sp.iterrows():
            ri=srow.values
            for j,orow in op.iterrows():
                rj=orow.values
                rij=ri-rj
                D=np.sqrt(np.dot(rij,rij))
                if D<minD:
                    minD=D
        return minD

    def rotate(self,R):
        """rotate Rotates all coordinate vectors by rotation matrix R

        :param R: rotation matrix (3x3)
        :type R: numpy.ndarray
        """
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            newri=np.matmul(R,ri)
            self.A.loc[i,'posX':'posZ']=newri

    def translate(self,L):
        """translate Translates all coordinate vectors by displacement vector L

        :param L: displacement vector (nm)
        :type L: numpy.ndarray
        """
        sp=self.A[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
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

    def minmax(self):
        sp=self.A[['posX','posY','posZ']]
        return np.array([sp.posX.min(),sp.posY.min(),sp.posZ.min()]),np.array([sp.posX.max(),sp.posY.max(),sp.posZ.max()])

    def checkbox(self):
        mm,MM=self.minmax()
        bb=self.box.diagonal()
        return mm<bb,MM>bb

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
            assert name in self.A.columns,f'{name} not found in attributes\n{df.columns}'
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
                ai,aj,o=b
                idx_to_delete.extend(self.sacH(ai,aj,T,rename=rename))
        return idx_to_delete

    def sacH(self,ai,aj,T,rename=False):
        """sacH Find the two H's, one bound to ai, the other to aj, that are closest to each other

        :param ai: index of one atom in bond
        :type ai: int
        :param aj: index of other atom in bond
        :type aj: int
        :param T: global topology
        :type T: Topology
        :param rename: whether to rename remaining H atoms bound to ai and aj so that it appears highest-sorted by name atoms are found, defaults to False
        :type rename: bool, optional
        :return: global indexes of two H atoms
        :rtype: list
        """
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
            # reverse sort names of hydrogen ligands by their number
            i_avails=list(sorted(i_Hpartners.values(),key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0')))[:-1]
            j_avails=list(sorted(j_Hpartners.values(),key=lambda x: int(x.split('H')[1] if x.split('H')[1]!='' else '0')))[:-1]
            logger.debug(f'i_avails {i_avails}')
            logger.debug(f'j_avails {j_avails}')
            # remove the globalIdx of the sacrificial H's from their atom's dictionaries of H-atoms
            del i_Hpartners[ih]
            del j_Hpartners[jh]
            Top=T.D['atoms']
            Cor=self.A
            # for all remaining H neighbor globalIdx of each atom, rename starting from lowest number
            for h in i_Hpartners:
                i_Hpartners[h]=i_avails.pop(0)
                Top.iloc[h-1,Top.columns=='atom']=i_Hpartners[h]
                Cor.iloc[h-1,Cor.columns=='atomName']=i_Hpartners[h]
                logger.debug(f'i: changed name of {h} to {i_Hpartners[h]}')
            for h in j_Hpartners:
                j_Hpartners[h]=j_avails.pop(0)
                Top.iloc[h-1,Top.columns=='atom']=j_Hpartners[h]
                Cor.iloc[h-1,Cor.columns=='atomName']=j_Hpartners[h]
                logger.debug(f'j: changed name of {h} to {j_Hpartners[h]}')
        # this makes sure that it always looks like the same atom was deleted
        return [ih,jh] # return the globalIdx's of the two sacrificial H's

    def delete_atoms(self,idx=[],reindex=True):
        """delete_atoms Deletes atoms whose global indices appear in the list idx.
        If parameter 'reindex' is true, then the global indices 
        are recalculated so that they are sequential starting at 1 with no
        gaps, and two new columns are added to self.DF:
          - 'oldGlobalIdx' contains the global index values before 
             the deletion.
          - 'globalIdxShift' is the change from the old to the new
             global index for each atom.

        :param idx: list of atom indexes to delete, defaults to []
        :type idx: list, optional
        :param reindex: reindex remaining atoms, defaults to True
        :type reindex: bool, optional
        """
        # logger.debug(f'Coordinates:delete_atoms {idx}')
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

    def write_gro(self,filename):
        """write_gro Write coordinates and if present, velocities, to a Gromacs-format coordinate file

        :param filename: name of file to write
        :type filename: str
        """
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
                logger.debug('Writing Gromacs coordinates file but boxsize is not set.')
            f.write(f'{self.box[0][0]:10.5f}{self.box[1][1]:10.5f}{self.box[2][2]:10.5f}')
            # output off-diagonals only if at least one of them is non-zero
            x,y=self.box.nonzero()
            if not all(x==y):
                f.write(f'{self.box[0][1]:10.5f}{self.box[0][2]:10.5f}')
                f.write(f'{self.box[1][0]:10.5f}{self.box[1][2]:10.5f}')
                f.write(f'{self.box[2][0]:10.5f}{self.box[2][1]:10.5f}')
            f.write('\n')

    def write_mol2(self,filename,bondsDF=pd.DataFrame(),molname='',other_attributes=pd.DataFrame()):
        """write_mol2 Write a mol2-format file from coordinates, and optionally, a bonds DataFrame
            provided externally and passed in as "bondsDF" (typically this would be
            from a Topology instance).

        :param filename: name of file name to write
        :type filename: str, optional
        :param bondsDF: dataframe of bonds ['ai','aj'], defaults to pd.DataFrame()
        :type bondsDF: pandas.DataFrame, optional
        :param molname: name of molecule, defaults to ''
        :type molname: str, optional
        :param other_attributes: auxiliary dataframe of attributes, defaults to pd.DataFrame()
        :type other_attributes: pandas.DataFrame, optional
        """
        # logger.debug(f'write_mol2 {filename}')
        acopy=self.A.copy()
        if bondsDF.empty and self.mol2_bonds.empty:
            logger.warning(f'Cannot write any bonds to MOL2 file {filename}')
        if not self.mol2_bonds.empty and bondsDF.empty:
            bdf=self.mol2_bonds
        elif not bondsDF.empty and self.mol2_bonds.empty:
            bdf=bondsDF
        else:
            logger.info('Coordinates.write_mol2 provided with both a bondsDF parameter and a mol2_bonds attribute')
            logger.info('Using the parameter')
            bdf=bondsDF
        for i in other_attributes.columns: #self.mol2_atom_attributes:
            # logger.debug(f'importing/overwriting other_attribute {i}...')
            acopy[i]=other_attributes[i]
        # logger.debug(f'Updated [ atoms ]:\n{acopy.to_string()}')
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
                nBonds=bdf.shape[0]
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
                    logger.debug(f'Mol2 bonds from outside')
                    bdf=bondsDF[['bondIdx','ai','aj','order']].copy()
                    bdf['bondIdx']=bdf['bondIdx'].astype(int)
                    bdf['ai']=bdf['ai'].astype(int)
                    bdf['aj']=bdf['aj'].astype(int)
                    f.write(bdf.to_string(columns=self.mol2_bond_attributes,header=False,index=False,formatters=bondformatters))
                elif not self.mol2_bonds.empty:
                    logger.debug(f'write_mol2 ({filename}): Mol2 bonds from mol2_bonds attribute')
                    f.write(self.mol2_bonds.to_string(columns=self.mol2_bond_attributes,header=False,index=False,formatters=bondformatters))
                f.write('\n')
                ''' write substructure section '''
                f.write('@<TRIPOS>SUBSTRUCTURE\n')
                f.write(rdf.to_string(header=False,index=False,formatters=substructureformatters))
                f.write('\n')
   