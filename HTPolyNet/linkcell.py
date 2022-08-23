import numpy as np
import pandas as pd
import logging
from itertools import product
from multiprocessing import Pool
from functools import partial
logger=logging.getLogger(__name__)

class Linkcell:
    """ Handles the link-cell algorithm for searching for bonding partners within a
        cutoff distance from each other
    """
    def __init__(self,box=[],cutoff=None,pbc_wrapper=None):
        """__init__ Constructor for an empy instance

        :param box: system box size, defaults to []
        :type box: list, optional
        :param cutoff: cutoff distance, defaults to None
        :type cutoff: float, optional
        """
        self.box=box
        self.cutoff=cutoff
        self.pbc_wrapper=pbc_wrapper

    def create(self,cutoff,box,origin=np.array([0.,0.,0.])):
        """create Creates the link-cell structure in a previously initialized instance

        :param cutoff: cutoff distance
        :type cutoff: float
        :param box: box size
        :type box: numpy.ndarray
        :param origin: origin, defaults to np.array([0.,0.,0.])
        :type origin: numpy.ndarry, optional
        """
        if box.shape==(3,3):
            box=np.diagonal(box)
        self.cutoff=cutoff
        self.box=box
        self.origin=origin
        # number of cells along x, y, and z directions
        self.ncells=np.floor(self.box/self.cutoff).astype(int)
        # dimensions of one cell
        self.celldim=box/self.ncells
        # 3-d array of lower left corner as a 3-space point, indexed by i,j,k
        # initialized to all zeros, calculated below
        self.cells=np.zeros((*self.ncells,3))
        # 1-d array of (i,j,k) indices indexed by linear cell index (0...ncells-1)
        self.cellndx=np.array(list(product(*[np.arange(x) for x in self.ncells])))
        # 3-d array of lower left corner as a 3-space point, indexed by i,j,k
        for t in self.cellndx:
            i,j,k=t
            self.cells[i,j,k]=self.celldim*np.array([i,j,k])+self.origin
        # set up neighbor lists using linear indices
        self.make_neighborlists()
        logger.debug(f'Linkcell structure: {len(self.cellndx)} cells ({self.ncells}) dim {self.celldim}')

    def cellndx_of_point(self,R):
        """cellndx_of_point returns the (i,j,k) cell index of point R

        :param R: a 3-space point
        :type R: numpy.ndarray
        :return: (i,j,k) index of cell
        :rtype: (int,int,int)
        """
        wrapR,bl=self.pbc_wrapper(R)
        C=np.floor(wrapR*np.reciprocal(self.celldim)).astype(int)
        lowdim=(C<np.zeros(3).astype(int)).astype(int) # will never happen if R is wrapped
        hidim=(C>=self.ncells).astype(int) # could happen if exactly there
        if (any(lowdim) or any(hidim)):
            logger.warning(f'Warning: point {R} maps to out-of-bounds-cell {C} ({self.ncells})')
            logger.warning(f'box: {self.box}')
        C+=lowdim
        C-=hidim
        return C

    def point_in_cellndx(self,R,C):
        """point_in_cellndx returns True if point R is located in cell with (i,j,k) index C

        :param R: 3-space point
        :type R: numpy.ndarray
        :param C: (i,j,k) cell index
        :type C: (int,int,int)
        :return: 
            - True: R in C
            - False: R not in C
        :rtype: boolean
        """
        LL,UU=self.corners_of_cellndx(C)
        wrapR,bl=self.pbc_wrapper(R)
        return all(wrapR<UU) and all(wrapR>=LL)

    def corners_of_cellndx(self,C):
        """corners_of_cellndx returns the lower-left and upper-right corners of cell with (i,j,k) index C, as an array of 3-space points

        :param C: (i,j,k) index
        :type C: (int,int,int)
        :return: 2x3 array of lower-left and upper-right corner coordinates
        :rtype: numpy.ndarray
        """
        LL=self.cells[C[0],C[1],C[2]]
        UU=LL+self.celldim
        return np.array([LL,UU])

    def cellndx_in_structure(self,C):
        """cellndx_in_structure test to see if (i,j,k) index given is within the established linkcell structure

        :param C: (i,j,k) index
        :type C: (int,int,int)
        :return:  - True: C is in the linkcell structure
                  - False: C is not in the linkcell structure
        :rtype: boolean
        """
        return all(np.zeros(3)<=C) and all(C<self.ncells)

    def ldx_of_cellndx(self,C):
        """ldx_of_cellndx return scalar index of cell with (i,j,k)-index C

        :param C: (i,j,k)-index
        :type C: (int,int,int)
        :return: scalar index of C
        :rtype: int
        """
        nc=self.ncells
        xc=C[0]*nc[1]*nc[2]+C[1]*nc[1]+C[2]
        return xc

    def cellndx_of_ldx(self,i):
        return self.cellndx[i]

    def populate_par(self,adf):
        """populate_par populate the linkcell structure by setting the "linkcell_idx" attribute of each atom in the coordinates dataframe adf

        :param adf: gromacs-format coordinate dataframe
        :type adf: pandas.DataFrame
        :return: modified dataframe
        :rtype: pandas.DataFrame
        """
        sp=adf[['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            C=self.cellndx_of_point(ri)
            idx=self.ldx_of_cellndx(C)
            adf.loc[i,'linkcell_idx']=idx
        return adf

    def _return_list_lens(self,idx_list,mlists):
        """_return_list_lens returns a list of lengths of lists in list mlists

        :param idx_list: indicies of mlists that will be tallied
        :type idx_list: list of ints
        :param mlists: list of lists, each a list of integer atom global indices
        :type mlists: list of list of ints
        :return: list of lengths of lists
        :rtype: list of ints
        """
        return [len(mlists[i]) for i in idx_list]

    def populate(self,Coordinates,ncpu=1):
        """populate Populates linkcell structure

        :param Coordinates: Coordinates instance from which atom coordinates are taken
        :type Coordinates: Coordinates
        :param ncpu: number of processors to split the operation over, defaults to 1
        :type ncpu: int, optional
        :raises Exception: dies if a point's assigned (i,j,k) cell is outside the cell structure (this would mean the atom's position is outside the periodic box, which is an error.  Atom coordinates are always to be held in wrapped configuration, but be careful: gro files unwrap!)
        """
        N=Coordinates.A.shape[0]
        logger.debug(f'Linkcell: assigning cell indices to {N} atoms in {self.box}...')
        Coordinates.set_atomset_attribute('linkcell_idx',-1*np.ones(N).astype(int))
        self.memberlists=[[] for _ in range(self.cellndx.shape[0])]
        ess='s' if ncpu>1 else ''
        logger.debug(f'Linkcell assignment will use {ncpu} processor{ess}')
        p=Pool(processes=ncpu)
        adf_split=np.array_split(Coordinates.A,ncpu)
        result=p.map(partial(self.populate_par),adf_split)
        p.close()
        p.join()
        Coordinates.A=pd.DataFrame()
        for a in result:
            Coordinates.A=pd.concat((Coordinates.A,a))

        idx_list=Coordinates.A['globalIdx'].to_list()
        for i in idx_list:
            lc_idx=Coordinates.get_atom_attribute('linkcell_idx',{'globalIdx':i})
            try:
                self.memberlists[lc_idx].append(i)
            except:
                logger.debug(f'Linear linkcell index {lc_idx} of atom {i} is out of range.\ncellndx.shape[0] is {self.cellndx.shape[0]}\nThis is a bug.')
                raise Exception

        # for i in range(len(self.memberlists)):
        #     logger.debug(f'{i} {len(self.memberlists[i])}')

        idx_list=list(range(len(self.memberlists)))
        p=Pool(processes=ncpu)
        idx_list_split=np.array_split(idx_list,ncpu)
        result=p.map(partial(self._return_list_lens,mlists=self.memberlists),idx_list_split)
        p.close()
        p.join()
        result=np.array([item for sublist in result for item in sublist])
        avg_cell_pop=result.mean()
        min_cell_pop=int(result.min())
        max_cell_pop=int(result.max())
        logger.debug(f'Avg/min/max cell pop: {avg_cell_pop:>8.3f}/{min_cell_pop:>8d}/{max_cell_pop:>8d}')
        # logger.debug(f'Linkcell.populate() ends.')

    def make_neighborlists(self):
        self.neighborlists=[[] for _ in range(self.cellndx.shape[0])]
        for C in self.cellndx:
            idx=self.ldx_of_cellndx(C)
            for D in self.neighbors_of_cellndx(C):
                if self.ldx_of_cellndx(D)!=idx:
                    self.neighborlists[idx].append(self.ldx_of_cellndx(D))

    def make_memberlists(self,cdf):
        self.memberlists=[[] for _ in range(self.cellndx.shape[0])]
        rdf=cdf[cdf['linkcell_idx']!=-1]
        # logger.debug(f'Generated {len(self.memberlists)} empty memberlists.')
        for i,r in rdf.iterrows():
            cidx=r['linkcell_idx']
            idx=r['globalIdx']
            self.memberlists[cidx].append(idx)
        rl=np.array([len(self.memberlists[i]) for i in range(self.cellndx.shape[0])])
        assert int(rl.sum())==rdf.shape[0] # check to make sure all atoms are counted
        avg_cell_pop=rl.mean()
        min_cell_pop=int(rl.min())
        max_cell_pop=int(rl.max())
        logger.debug(f'Avg/min/max cell pop: {avg_cell_pop:>8.3f}/{min_cell_pop:>8d}/{max_cell_pop:>8d}')

    def neighbors_of_cellndx(self,Ci):
        assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
        retlist=[]
        dd=np.array([-1,0,1])
        S=list(product(dd,dd,dd))
        S.remove((0,0,0))
        for s in S:
            nCi=Ci+np.array(s)
            for d in range(3):
                if nCi[d]==self.ncells[d]:
                    nCi[d]=0
                elif nCi[d]==-1:
                    nCi[d]=self.ncells[d]-1
            retlist.append(nCi)
        assert len(retlist)==(len(dd)**3)-1
        return retlist

    def searchlist_of_ldx(self,i):
        assert i!=-1
        retlist=[]
        C=self.cellndx[i]
        for c in self.neighbors_of_cellndx(C):
            retlist.append(self.ldx_of_cellndx(c))
        assert len(retlist)==27,f'Error: not counting enough neighbor cells'
        return retlist

    def are_cellndx_neighbors(self,Ci,Cj):
        assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
        assert self.cellndx_in_structure(Cj),f'Error: cell {Cj} outside of cell structure {self.ncells}'
        dij=Ci-Cj
        low=(dij<-self.ncells/2).astype(int)
        hi=(dij>self.ncells/2).astype(int)
        dij+=(low-hi)*self.ncells
        return all([x in [-1,0,1] for x in dij])

    def are_ldx_neighbors(self,ildx,jldx):
        # should never call this for atoms with unset lc indices
        assert ildx!=-1
        assert jldx!=-1
        return jldx in self.neighborlists[ildx]

if __name__=='__main__':
    box=np.array([4.7,4.7,4.7])
    cutoff=0.5
    L=Linkcell()
    L.create(cutoff,box)
    check1=[]
    for i,C in enumerate(L.cellndx):
        check1.append(i==L.ldx_of_cellndx(C))
    print(f'check1 {all(check1)}')
    check2=[]
    for i,C in enumerate(L.cellndx):
        n=0
        for j,dldx in enumerate(L.neighborlists[i]):
            D=L.cellndx[dldx]
            check2.append(L.are_cellndx_neighbors(C,D))
            n+=1
        
    print(f'check2 {all(check2)}')
