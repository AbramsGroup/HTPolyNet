import numpy as np
import logging
from itertools import product

class Linkcell:
    def __init__(self,box=[],cutoff=None):
        self.box=box
        self.cutoff=None

    def create(self,cutoff,box,origin=np.array([0.,0.,0.])):
        logging.debug(f'Linkcell.create() begins')
        if box.shape==(3,3):
            box=np.diagonal(box)
        self.cutoff=cutoff
        self.box=box
        self.origin=origin
        self.ncells=np.floor(self.box/self.cutoff).astype(int)
        self.celldim=box/self.ncells
        self.cells=np.zeros((*self.ncells,3))
        self.cellidx=np.array(list(product(*[np.arange(x) for x in self.ncells])))
        for t in self.cellidx:
            i,j,k=t
            self.cells[i,j,k]=self.celldim*np.array([i,j,k])+self.origin
        logging.debug(f'Linkcell.create() ends; structure has {len(self.cellidx)} cells ({self.ncells}) dim {self.celldim}')

    def point_in_box(self,R):
        LL=self.origin
        UU=self.origin+self.box
        gt=all(R>=LL)
        lt=all(R<UU)
        return lt and gt
    
    def wrap_point(self,R):
        if self.point_in_box(R):
            return R
        nR=R.copy()
        for d in range(0,3):
            if nR[d]<self.origin[d]:
                nR[d]+=self.box[d]
            elif nR[d]>self.origin[d]+self.box[d]:
                nR[d]-=self.box[d]
        return nR

    def cell_of_point(self,R):
        return np.floor(self.wrap_point(R)*np.reciprocal(self.celldim)).astype(int)

    def point_in_cell(self,R,C):
        LL,UU=self.corners_of_cell(C)
        gt=all(R>=LL)
        lt=all(R<UU)
        return lt and gt

    def corners_of_cell(self,C):
        LL=self.cells[C[0],C[1],C[2]]
        UU=LL+self.celldim
        return np.array([LL,UU])

    def cell_in_structure(self,C):
        return all(np.zeros(3)<=C) and all(C<self.ncells)

    def index_of_cell(self,C):
        return self.ncells[0]*self.ncells[1]*C[0]+self.ncells[1]*C[1]+C[2]

    def populate(self,Coordinates):
        ''' Set the linkcell-idx attribute of every atom and create
            each cell's list of members (each element is a globalIdx) '''
        N=Coordinates.A.shape[0]  # Coordinates.N should also work
        logging.debug('Linkcell.populate() assigning cell indices to {N} atoms in {self.box}')
        Coordinates.set_atomset_attribute('linkcell-idx',-1*np.ones(N).astype(int))
        # self.memberlists=[[]]*len(self.cellidx)
        for i in range(N):
            R=Coordinates.get_R(i+1)
            C=self.cell_of_point(self.wrap_point(R))
            idx=self.index_of_cell(C)
            Coordinates.set_atom_attribute('linkcell-idx',idx,{'globalIdx':i+1})
            # self.memberlists[idx].append(i+1)
        logging.debug('Linkcell.populate() ends.')

    def next_j(self,i,Coordinates):
        pass

    def neighbors_of_cell(self,Ci):
        assert self.cell_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
        retlist=[]
        for d in range(3):
            p=np.zeros(3).astype(int)
            p[d]=1
            for s in [-1,1]:
                nCi=Ci+s*p
                if nCi[d]==self.ncells[d]:
                    nCi[d]=0
                elif nCi[d]==-1:
                    nCi[d]=self.ncells[d]-1
                retlist.append(nCi)
        return retlist

    def are_cells_neighbors(self,Ci,Cj):
        assert self.cell_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
        assert self.cell_in_structure(Cj),f'Error: cell {Cj} outside of cell structure {self.ncells}'
        oneaway=np.array([False,False,False])
        for d in range(0,3):
            dd=Ci[d]-Cj[d]
            if dd>self.ncells[d]/2:
                dd-=self.ncells[d]
            elif dd<-self.ncells[d]/2:
                dd+=self.ncells[d]
            add=np.abs(dd)
            oneaway[d]=add==1
        return oneaway.astype(int).sum()==1

    def are_cellidx_neighbors(self,cix,cjx):
        Ci=self.cellidx[cix]
        Cj=self.cellidx[cjx]
        return self.are_cells_neighbors(Ci,Cj)

    def update(self,adf):
        pass

if __name__=='__main__':
    box=np.array([10.,10.,10.])
    cutoff=2.1
    L=Linkcell()
    L.create(cutoff,box)
    c2=np.array([0,0,0])
    for C in L.neighbors_of_cell(c2):
        print(c2,C,L.are_cells_neighbors(C,c2))
    
    p=np.array([3.4,5.7,7.9])
    C=L.cell_of_point(p)
    i=L.index_of_cell(C)
    print(p,C,L.point_in_cell(p,C),i,L.cellidx[i])
