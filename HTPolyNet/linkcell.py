import numpy as np
import logging
from itertools import product

class Linkcell:
    def __init__(self,box=[],cutoff=None):
        self.box=box
        self.cutoff=None

    def create(self,cutoff,box,origin=np.array([0.,0.,0.])):
        if box.shape==(3,3):
            box=np.diagonal(box)
        self.box=box
        self.ncells=np.floor(box/cutoff).astype(int)
        logging.info(f'Linkcells: ncells {self.ncells}')
        self.celldim=box/self.ncells
        logging.info(f'Linkcells: celldim {self.celldim}')
        self.cells=np.zeros((*self.ncells,3))
        self.cellidx=np.array(list(product(*[np.arange(x) for x in self.ncells])))
        for t in self.cellidx:
            i,j,k=t
            self.cells[i,j,k]=self.celldim*np.array([i,j,k])

    def cell_of(self,R):
        return np.floor(R*np.reciprocal(self.celldim)).astype(int)

    def corners_of(self,C):
        LL=self.cells[C[0],C[1],C[2]]
        UU=LL+self.celldim
        return np.array([LL,UU])

    def populate(self,adf):
        pass

    def neighbors(self,Ci,Cj):
        # print(np.abs(Ci[1]-Cj[1])%self.ncells[1])
        inx=np.abs(Ci[0]-Cj[0])%self.ncells[0]==1 and (Ci[1]==Cj[1] and Ci[2]==Cj[2])
        iny=np.abs(Ci[1]-Cj[1])%self.ncells[1]==1 and (Ci[2]==Cj[2] and Ci[0]==Cj[0])
        inz=np.abs(Ci[2]-Cj[2])%self.ncells[2]==1 and (Ci[0]==Cj[0] and Ci[1]==Cj[1])
        # print(inx,iny,inz)
        return any([inx,iny,inz])

    def update(self,adf):
        pass

if __name__=='__main__':
    box=np.array([10.,10.,10.])
    cutoff=2.1
    L=Linkcell()
    L.create(cutoff,box)
    print(L.ncells)
    print(L.celldim)
    print(L.cellidx)
    
    R=np.array([2.3,4.2,5.3])
    print(L.corners_of(L.cell_of(R)))
    c1=np.array([0,0,0])
    c2=np.array([0,4,0])
    print(L.neighbors(c1,c2))