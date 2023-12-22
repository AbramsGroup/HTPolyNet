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
from HTPolyNet.dataframetools import *
from HTPolyNet.matrix4 import Matrix4

logger=logging.getLogger(__name__)

GRX_ATTRIBUTES     =[  'z','nreactions','reactantName','sea_idx','chain','chain_idx','molecule','molecule_name']
"""Extended atom attributes

    - 'z' number of sacrificial H's on atom 
    - 'nreactions' number of H's sacrificed so far to form bonds
    - 'reactantName' name of most recent reactant to which atom belonged
    - 'sea_idx' index of the group of symmetry-related atoms this atom belongs to (atoms with the same sea_idx in the same resid are considered symmetry-equivalent)
    - 'chain' index of the unique chain this atom belongs to
    - 'chain_idx' index of this atom within this chain
    - 'molecule' index of the unique molecule this atom belongs to
    - 'molecule_name' name of that molecule
"""
GRX_GLOBALLY_UNIQUE=[False,       False,         False,     True,  True,      False,      True,          False]
GRX_UNSET_DEFAULTS =[    0,           0,       'UNSET',       -1,    -1,          -1,       -1,        'UNSET']

def dfrotate(df:pd.DataFrame,R):
    """dfrotate applies rotation matrix R to coordinates in dataframe

    :param df: coordinates dataframe; must have 'posX', 'posY', and 'posZ' columns
    :type df: pd.DataFrame
    :param R: rotation matrix
    :type R: np.ndarray((3,3))
    """
    for i,srow in df.iterrows():
        ri=srow[['posX','posY','posZ']].values
        newri=np.matmul(R,ri)
        df.loc[i,'posX':'posZ']=newri

class Coordinates:
    """ Handles atom coordinates.

    The primary object is `A`, a pandas DataFrame with one row per atom.  Each atom has attributes that may be found in a `gro` file and/or a `mol2` file, along with so-called extended attributes, which are used solely  by HTPolyNet.  

    """
    gro_attributes = ['resNum', 'resName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    """GRO format atom attributes
    
    - 'resNum' unique index of residue to which atom belongs
    - 'resName' name of that residue (usually a 3-letter designation)
    - 'atomName' name of this atom, must be unique within a residue
    - 'globalIdx' global index of atom in whole system
    - 'posX', 'posY', 'posZ' cartesian coordinates
    - 'velX', 'velY', 'velZ' cartesian velocities
    """
    mol2_atom_attributes = ['globalIdx','atomName','posX','posY','posZ','type','resNum','resName','charge']
    """MOL2 format atom attributes
    
    - 'globalIdx' global index of atom in whole system
    - 'atomName' name of this atom, must be unique within a residue
    - 'posX', 'posY', 'posZ' cartesian coordinates
    - 'resNum' unique index of residue to which atom belongs
    - 'resName' name of that residue (usually a 3-letter designation)
    - 'charge' charge on atom
    """
    mol2_bond_attributes = ['bondIdx','ai','aj','order']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_attributes, [int, int, int, str])}

    def __init__(self,name=''):
        """__init__ contructs an empty Coordinates object

        :param name: a name string, defaults to ''
        :type name: str, optional
        """
        self.name=name
        self.metadat={}
        self.N=0
        self.A=pd.DataFrame()
        self.mol2_bonds=pd.DataFrame()
        self.mol2_bondlist=Bondlist()
        self.linkcell=Linkcell(pbc_wrapper=self.wrap_point)
        self.empty=True
        self.box=np.zeros((3,3))
        self.grx_attributes=GRX_ATTRIBUTES
        self.parent=None
        
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
                lc_globalIdx=1
                for x in data[2:-1]:
                    series['resNum'].append(int(x[0:5].strip()))
                    series['resName'].append(x[5:10].strip())
                    series['atomName'].append(x[10:15].strip())
                    ''' if formatted correctly, globalIdx is row index + 1 always! '''
                    series['globalIdx'].append(lc_globalIdx)
                    lc_globalIdx+=1
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

    @classmethod
    def fcc(cls,a,nc=[1,1,1]):
        """fcc generates a Coordinates object that represents an FCC crystal

        :param a: lattice parameter
        :type a: float
        :param nc: number of unit cells in the three lattice vector directions, defaults to [1,1,1]
        :type nc: list, optional
        :return: a Coordinates object
        :rtype: Coordinates
        """
        inst=cls()
        basis=np.identity(3)*a
        base_atoms=np.array([[0.0,0.0,0.0],[0.5,0.5,0.0],[0.5,0.0,0.5],[0.0,0.5,0.5]])*a
        n=0
        p=[]
        x=[list(range(n)) for n in nc]
        for i,j,k in product(*x):
            ll=np.dot(basis,np.array([i,j,k]))
            for m in range(len(base_atoms)):
                p.append(ll+base_atoms[m])
        # print(p)
        posn=np.array(p)
        # print(posn)
        N=len(posn)
        adf=inst.A
        adf['globalIdx']=list(range(1,N+1))
        adf['atomName']='AL'
        adf['resNum']=list(range(1,N+1))
        adf['posX']=posn[:,0]
        adf['posY']=posn[:,1]
        adf['posZ']=posn[:,2]
        adf['resName']='MET'
        inst.N=N
        return inst

    def claim_parent(self,parent):
        self.parent=parent

    def set_box(self,box:np.ndarray):
        """set_box Set the box size from box

        :param box: 3-by-1 or 3-by-3 box size matrix
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

    def subcoords(self,sub_adf:pd.DataFrame):
        """subcoords generates a new Coordinates object to hold the atoms dataframe in 'sub_adf' parameter

        :param sub_adf: an atom dataframe
        :type sub_adf: pd.DataFrame
        :return: a new Coordinates object
        :rtype: Coordinates
        """
        newC=Coordinates()
        newC.set_box(self.box)
        newC.A=sub_adf
        newC.N=sub_adf.shape[0]
        return newC
    
    def reconcile_subcoords(self,subc,attr):
        """moves all values of attribute name contained in attr from Coordinates object subc to self

        :param subc: A separate, independent Coordinates object
        :type subc: Coordinates
        :param attr: attribute name whose value is to be copied from subc to self
        :type attr: str
        """
        jdx=list(subc.A.columns).index(attr)
        for r in subc.A.itertuples(index=False):
            idx=r.globalIdx
            lc_idx=r[jdx]
            self.A.loc[idx-1,attr]=lc_idx

    def unwrap(self,P,O,pbc):
        """unwrap shifts point P to its unwrapped closest periodic image to point O

        :param P: a point
        :type P: np.ndarray(3,float)
        :param O: origin
        :type O: np.ndarray(3,float)
        :param pbc: directions in which pbc are applied
        :type pbc: np.ndarray(3,int)
        :return: a point
        :rtype: nd.ndarray(3,float)
        """
        ROP=self.mic(O-P,pbc)
        PCPI=O-ROP
        return PCPI

    def pierces(self,B:pd.DataFrame,C:pd.DataFrame,pbc):
        """pierces determines whether or not bond represented by two points in B pierces ring represented by N points in C

        :param B: two points defining a bond
        :type B: pd.DataFrame
        :param C: N points defining a ring
        :type C: pd.DataFrame
        :param pbc: periodic boundary condition flags, one per dimension
        :type pbc: list
        :return: True if ring is pierced
        :rtype: bool
        """
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
        """linkcell_initialize initializes link-cell structure for ring-pierce testing

        :param cutoff: cutoff radius for link-cell structure, defaults to 0.0
        :type cutoff: float, optional
        :param ncpu: number of cpus to use to populate link-cell structure, defaults to 1
        :type ncpu: int, optional
        :param populate: If True, an actual population of the link-cell structure is performed; defaults to True
        :type populate: bool, optional
        :param force_repopulate: If True and this link-cell structure is already populated, a repopulation is performed based on the current Coordinates; defaults to False
        :type force_repopulate: bool, optional
        :param save: If true, all atoms' linkcell_idx attributes are written to an output file; defaults to True
        :type save: bool, optional
        """
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
                # we only populate with atoms whose positions will be needed in interatomic
                # distance calculations; these are those (a) in rings, or (b) are reactive
                sc=self.subcoords(self.A[(self.A['globalIdx'].isin(self.parent.rings.all_atoms()))|(self.A['z']>0)].copy())
                self.linkcell.populate(sc,ncpu=ncpu)
                self.reconcile_subcoords(sc,'linkcell_idx')
                if save:
                    self.write_atomset_attributes(['linkcell_idx'],lc_file)

    def linkcelltest(self,i,j):
        """linkcelltest returns True if atoms i and j are within potential interaction
            range based on current link-cell structure

        :param i: An atom index
        :type i: int
        :param j: another atom index
        :type j: int
        :return: True if i and j are in the same cell or in neighboring cells
        :rtype: bool
        """
        ci=self.get_atom_attribute('linkcell_idx',{'globalIdx':i})
        cj=self.get_atom_attribute('linkcell_idx',{'globalIdx':j})
        if ci==cj:
            return True
        if self.linkcell.are_ldx_neighbors(ci,cj):
            return True
        return False

    def geometric_center(self):
        """geometric_center computes and returns the geometric center of the atoms in self's A dataframe

        :return: geometric center
        :rtype: np.ndarray(3,float)
        """
        a=self.A
        return np.array([a.posX.mean(),a.posY.mean(),a.posZ.mean()])

    def rij(self,i,j,pbc=[1,1,1]):
        """rij compute distance between atoms i and j

        :return: distance between i and j
        :rtype: float
        """
        if np.any(pbc) and not np.any(self.box):
            logger.warning('Interatomic distance calculation using PBC with no boxsize set.')
        ri=self.get_R(i)
        rj=self.get_R(j)
        Rij=self.mic(ri-rj,pbc)
        return np.sqrt(Rij.dot(Rij))

    def mic(self,r,pbc):
        """mic applies minimum image convention to displacement vector r

        :param r: displacement vector
        :type r: np.ndarray(3,float)
        :param pbc: periodic boundary condition 
        :type pbc: _type_
        :return: _description_
        :rtype: _type_
        """
        '''  '''
        for c in range(0,3):
            if pbc[c]:
                hbx=self.box[c][c]/2
                while r[c]<-hbx:
                    r[c]+=self.box[c][c]
                while r[c]>hbx:
                    r[c]-=self.box[c][c]
        return r

    def wrap_point(self,ri):
        """wrap_point wraps point ri into the central periodic image

        :param ri: a point
        :type ri: np.ndarray(3,float)
        :return: a tuple containing (1) the wrapped point and (2) number of box lengths required to wrap this point, per dimension
        :rtype: tuple
        """
        R=ri.copy()
        box_lengths=np.array([0,0,0],dtype=int)
        for i in range(3):
            while R[i]<0:
                R[i]+=self.box[i][i]
                box_lengths[i]+=1
            while R[i]>=self.box[i][i]:
                R[i]-=self.box[i][i]
                box_lengths[i]-=1
        return R,box_lengths

    def wrap_coords(self):
        """wrap_coords Wraps all atomic coordinates into box
        """
        assert np.any(self.box),f'Cannot wrap if boxsize is not set: {self.box}'
        sp=self.A[['posX','posY','posZ']]
        boxL=[]
        for i,srow in sp.iterrows():
            p,box_lengths=self.wrap_point(srow.values)
            self.A.loc[i,'posX':'posZ']=p
            boxL.append(box_lengths)
        boxL=np.array(boxL)
        self.A['boxLx']=boxL[:,0]
        self.A['boxLy']=boxL[:,1]
        self.A['boxLz']=boxL[:,2]
        # logger.debug(f'Wrapped {self._nwrap}/{self.A.shape[0]*3} coordinates.')

    def merge(self,other):
        """merge Merge two Coordinates objects

        :param other: the other Coordinates object
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
            logger.debug(f'Read from {filename}\n{df.head().to_string()}')
        logger.debug(f'Merge:\n{self.A.head().to_string()}\nand\n{df.head().to_string()}')
        self.A=self.A.merge(df,how='outer',on='globalIdx')
        logger.debug(f'Result:\n{self.A.head().to_string()}')
        return attributes_read

    def set_atomset_attribute(self,attribute,srs):
        """set_atomset_attribute sets attribute of atoms to srs

        :param attribute: name of attribute
        :type attribute: str
        :param srs: scalar or list-like attribute values in same ordering as self.A
        :type srs: scalar or list-like
        """
        self.A[attribute]=srs

    def atomcount(self):
        """atomcount returns the number of atoms in the Coordinates object

        :return: number of atoms
        :rtype: int
        """
        return self.N

    def decrement_z(self,pairs):
        """decrement_z decrements value of z attributes of all atoms found in pairs

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
        """show_z_report generates a little text-based histogram displaying number of atoms with each value of z between 0 and 3; atoms are keyed by resname:atomname.
        """
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
        """return_bond_lengths returns an ordered list of bond lengths computed based on bonds indicated by the parallel 'ai' and 'aj' columns of the parameter dataframe bdf

        :param bdf: a pandas dataframe with 'ai' and 'aj' columns of atom indices indicating bonds
        :type bdf: pd.DataFrame
        :return: list of atom distances
        :rtype: list
        """
        lengths=[]
        for b in bdf.itertuples():
            lengths.append(self.rij(b.ai,b.aj))
        return lengths

    def add_length_attribute(self,bdf,attr_name='length'):
        """add_length_attribute computes bond lengths based on bonds indicated by the parallel 'ai' and 'aj' columns of the parameter dataframe bdf and stores result in a new column called attr_name

        :param bdf: a pandas dataframe with 'ai' and 'aj' columns of atom indices indicating bonds
        :type bdf: pd.DataFrame
        :param attr_name: name of length attribute column, defaults to 'length'
        :type attr_name: str, optional
        """
        lengths=[]
        for b in bdf.itertuples():
            lengths.append(self.rij(b.ai,b.aj))
        bdf[attr_name]=lengths

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
    
    def homog_trans(self,M:Matrix4,indices=[]):
        """applies homogeneous transformation matrix M [4x4] to coordinates

        :param M: homogeneous transformation matrix
        :type M: Matrix4
        """
        df=self.A
        for i,srow in df.iterrows():
            if len(indices)==0 or (srow['globalIdx'] in indices):
                ri=np.array(list(srow[['posX','posY','posZ']].values))
                df.loc[i,'posX':'posZ']=M.transform(ri)

    def rotate(self,R):
        """rotate Rotates all coordinate vectors by rotation matrix R

        :param R: rotation matrix (3x3)
        :type R: numpy.ndarray
        """
        M=Matrix4(R)
        self.homog_trans(M)
    #     sp=self.A[['posX','posY','posZ']]
    #     for i,srow in sp.iterrows():
    #         ri=srow.values
    #         newri=np.matmul(R,ri)
    #         self.A.loc[i,'posX':'posZ']=newri

    def translate(self,L):
        """translate Translates all coordinate vectors by displacement vector L

        :param L: displacement vector (nm)
        :type L: numpy.ndarray
        """
        M=Matrix4(L)
        self.homog_trans(M)
    #     sp=self.A[['posX','posY','posZ']]
    #     for i,srow in sp.iterrows():
    #         self.A.loc[i,'posX':'posZ']=srow.values+L

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
        """minmax returns the coordinates of the atoms at the lower-leftmost and upper-rightmost positions in the constellation of points in the atoms dataframe

        :return: tuple of two points, lower-leftmost and upper-rightmost, respectively
        :rtype: tuple(np.ndarray(3,float),np.ndarray(3,float))
        """
        sp=self.A[['posX','posY','posZ']]
        return np.array([sp.posX.min(),sp.posY.min(),sp.posZ.min()]),np.array([sp.posX.max(),sp.posY.max(),sp.posZ.max()])

    def checkbox(self):
        """checkbox checks that the entire constellation of points in the atoms dataframe fits within the designated box for this Configuration object

        :return: True,True if both lower-leftmost and upper-rightmost points are within the box
        :rtype: tuple(bool,bool)
        """
        mm,MM=self.minmax()
        bb=self.box.diagonal()
        return mm<bb,MM>bb

    def get_idx(self,attributes):
        """get_idx returns the global atom index of the atom in the atoms dataframe which has the attribute:value pairs indicated in the parameter attributes

        :param attributes: dictionary of attribute:value pairs that identify the atom
        :type attributes: dict
        :return: global atom index
        :rtype: int
        """
        df=self.A
        return get_row_attribute(df,'globalIdx',attributes)
    
    def get_R(self,idx):
        """get_R return the cartesian position of atom with global index idx

        :param idx: global index of atom
        :type idx: int
        :return: its cartesian position
        :rtype: numpy.ndarray(3,float)
        """
        df=self.A
        return np.array(get_row_attribute(df,['posX','posY','posZ'],{'globalIdx':idx}))
    
    def get_atom_attribute(self,name,attributes):
        """get_atom_attribute return values of attributes listed in name from atoms specified by attribute:value pairs in attributes

        :param name: list of attributes whose values are to be returned
        :type name: list
        :param attributes: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attributes: dict
        :return: scalar or list of one or more return attribute values
        :rtype: list if name is a list; scalar otherwise
        """
        df=self.A
        if type(name)==list:
            assert all([i in self.A.columns for i in name])
        else:
            assert name in self.A.columns,f'{name} not found in attributes\n{df.columns}'
        return get_row_attribute(df,name,attributes)
    
    def spew_atom(self,attributes):
        """spew_atom outputs all attributes of atom identified by the attributes dict

        :param attributes: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attributes: dict  
        :return: stringified dataframe
        :rtype: string
        """
        df=self.A
        return get_row_as_string(df,attributes)

    def get_atoms_w_attribute(self,name,attributes):
        """get_atoms_w_attribute returns all rows of atoms dataframe and columns named in names of atoms identified by the attributes dict

        :param name: list of attributes to be in the rows that are returned
        :type name: list
        :param attributes: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attributes: dict
        :return: a dataframe segment
        :rtype: pd.DataFrame
        """
        df=self.A
        return get_rows_w_attribute(df,name,attributes)

    def set_atom_attribute(self,name,value,attributes):
        """set_atom_attribute set the attributes named in name to values named in values (names||values) for the set of atoms specified in the attributes dict

        :param name: list of names of attributes
        :type name: list
        :param value: list of values of attributes to be set
        :type value: list
        :param attributes: dictionary of attribute:value pairs that specify the set of atoms to be considered
        :type attributes: dict
        """
        df=self.A
        set_row_attribute(df,name,value,attributes)

    def has_atom_attributes(self,attributes):
        """has_atom_attributes returns True if all atoms in atoms dataframe have the attribute named in attributes

        :param attributes: list of attribute names to look for
        :type attributes:  list or list-like container
        :return: True if all atoms have the attributes
        :rtype: bool
        """
        df=self.A
        return all([name in df for name in attributes])

    def find_sacrificial_H(self,pairs,T,rename=False,explicit_sacH={}):
        """find_sacrificial_H identifies all appropriate sacrificial hydrogen atoms determined by the bonds indicated in the pairs list of 3-tuples (ai,aj,order)

        :param pairs: list of atom pairs/order tuples
        :type pairs: tuple(int,int,int)
        :param T: the global Topology
        :type T: Topology
        :param rename: If true, renames any remaining H atoms so that it appears as though highest-order named H atoms are the ones sacrificed, defaults to False
        :type rename: bool, optional
        :param explicit_sacH: dictionary of pre-chosen sacrificial H atoms, keyed by pair index, defaults to {}
        :type explicit_sacH: dict, optional
        :return: list of global atom indices to delete
        :rtype: list
        """
        idx_to_delete=[]
        for i,b in enumerate(pairs):
            if not i in explicit_sacH:
                ai,aj,o=b
                idx_to_delete.extend(self.sacH(ai,aj,T,rename=rename))
            else:
                idx_to_delete.extend(explicit_sacH[i])
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
        gaps, and two new columns are added to self.DF: 'oldGlobalIdx' contains the global index values before the deletion, and 'globalIdxShift' is the change from the old to the new global index for each atom.

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

    def write_gro(self,filename,grotitle=''):
        """write_gro Write coordinates and if present, velocities, to a Gromacs-format coordinate file

        :param filename: name of file to write
        :type filename: str
        """
        title=self.name if not grotitle else grotitle
        has_vel='velX' in self.A.columns
        with open(filename,'w') as f:
            f.write(title+'\n')
            f.write(f'{self.N:>5d}\n')
            # C-format: “%5i%5s%5s%5i%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f”
            # Note that the gro atom number is not used; gromacs assigns atom indicies based
            # on counting input lines!  We will wrap the index so that it only has 5 digits.
            atomformatters = [
                lambda x: f'{x:>5d}',
                lambda x: f'{x:<5s}',
                lambda x: f'{x:>5s}',
                lambda x: f'{x%10000:5d}']+[lambda x: f'{x:8.3f}']*3 + [lambda x: f'{x:8.4f}']*3
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
        bdf=pd.DataFrame()
        if bondsDF.empty and self.mol2_bonds.empty:
            logger.warning(f'Cannot write any bonds to MOL2 file {filename}')
        elif (not self.mol2_bonds.empty) and bondsDF.empty:
            bdf=self.mol2_bonds
        elif (not bondsDF.empty) and self.mol2_bonds.empty:
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
   