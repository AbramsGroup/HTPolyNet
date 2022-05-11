import numpy as np
import logging
from itertools import product

class Linkcell:
    ''' A class for implementing the link-cell algorithm for pairwise searches '''
    def __init__(self,box=[],cutoff=None):
        self.box=box
        self.cutoff=cutoff

    def create(self,cutoff,box,origin=np.array([0.,0.,0.])):
        logging.debug(f'Linkcell.create() begins')
        if box.shape==(3,3):
            box=np.diagonal(box)
        self.cutoff=cutoff
        self.box=box
        self.origin=origin
        self.ncells=np.floor(self.box/self.cutoff).astype(int)
        self.celldim=box/self.ncells
        self.cells=np.zeros((*self.ncells,3)) # array of 3-float cell corners
        self.cellndx=np.array(list(product(*[np.arange(x) for x in self.ncells]))) # list of 3-index tuples
        for t in self.cellndx:
            i,j,k=t
            self.cells[i,j,k]=self.celldim*np.array([i,j,k])+self.origin
        self.make_neighborlists()
        logging.debug(f'Linkcell.create() ends; structure has {len(self.cellndx)} cells ({self.ncells}) dim {self.celldim}')

    def point_in_box(self,R):
        LL=self.origin
        UU=self.origin+self.box
        gt=all(R>=LL)
        lt=all(R<UU)
        return lt and gt
    
    def wrap_point(self,R):
        ''' wraps a point that lies outside the overall box back into the box '''
        if self.point_in_box(R):
            return R
        nR=R.copy()
        for d in range(0,3):
            if nR[d]<self.origin[d]:
                nR[d]+=self.box[d]
            elif nR[d]>(self.origin[d]+self.box[d]):
                nR[d]-=self.box[d]
        return nR

    def cellndx_of_point(self,R):
        # if any(np.isclose(R,self.box)):
        #     logging.debug(f'point {R} close to box {self.box}')
        #     np.printoptions(precision=10)
        #     d=R-self.box
        #     for x in d:
        #         logging.debug(f'{x:.7e}')
        C=np.floor(R*np.reciprocal(self.celldim)).astype(int)
        lowdim=(C<0).astype(int) # will never happen if R is wrapped
        hidim=(C>=self.ncells).astype(int) # could happen if exactly there
        C+=lowdim
        C-=hidim
        # if any(C>=self.ncells):
        #     logging.debug(f'cellndx of {R} is {C} out of bounds of {self.ncells}')
        return C

    def point_in_cellndx(self,R,C):
        LL,UU=self.corners_of_cellndx(C)
        gt=all(R>=LL)
        lt=all(R<UU)
        return lt and gt

    def corners_of_cellndx(self,C):
        LL=self.cells[C[0],C[1],C[2]]
        UU=LL+self.celldim
        return np.array([LL,UU])

    def cellndx_in_structure(self,C):
        return all(np.zeros(3)<=C) and all(C<self.ncells)

    def ldx_of_cellndx(self,C):
        nc=self.ncells
        xc=C[0]*nc[1]*nc[2]+C[1]*nc[1]+C[2]
        return xc

    def populate(self,Coordinates):
        ''' Set the linkcell-idx attribute of every atom and create
            each cell's list of members (each element is a globalIdx) '''
        N=Coordinates.A.shape[0]  # Coordinates.N should also work
        logging.debug(f'Linkcell.populate() assigning cell indices to {N} atoms in {self.box}...')
        Coordinates.set_atomset_attribute('linkcell-idx',-1*np.ones(N).astype(int))
        self.memberlists=[[] for _ in range(self.cellndx.shape[0])]
        for i in range(N):
            R=Coordinates.get_R(i+1)
            C=self.cellndx_of_point(self.wrap_point(R))
            idx=self.ldx_of_cellndx(C)
            Coordinates.set_atom_attribute('linkcell-idx',idx,{'globalIdx':i+1})
            try:
                self.memberlists[idx].append(i+1)
            except:
                logging.debug(f'{idx} out of range? cellndx.shape[0] is {self.cellndx.shape[0]}')
                logging.debug(f'C {C} R {R} wrapped(R) {self.wrap_point(R)}')
                logging.debug(f'pointinbox {self.point_in_box(R)}')
        amm=np.array([0.,1e9,-1e9])
        for i in range(len(self.memberlists)):
            c=len(self.memberlists[i])
            amm[0]+=c
            if c<amm[1]:
                amm[1]=c 
            elif c>amm[2]:
                amm[2]=c
        amm[0]/=len(self.memberlists)
        logging.debug(f'Linkcell.populate() ends. Max: {amm[2]:.0f} min {amm[1]:.0f} avg {amm[0]:0.3f}')

    def make_neighborlists(self):
        self.neighborlists=[[] for _ in range(self.cellndx.shape[0])]
        for C in self.cellndx:
            idx=self.ldx_of_cellndx(C)
            for D in self.neighbors_of_cellndx(C):
                self.neighborlists[idx].append(self.ldx_of_cellndx(D))

    def make_memberlists(self,cdf):
        self.memberlists=[[] for _ in range(self.cellndx.shape[0])]
        logging.debug(f'Generated {len(self.memberlists)} empty memberlists {len(self.memberlists[0])}.')
        for i,r in cdf.iterrows():
            cidx=r['linkcell-idx']
            idx=r['globalIdx']
            self.memberlists[cidx].append(idx)
            # logging.debug(f'Added {idx} as element {len(self.memberlists[cidx])} to cell {cidx}')
        amm=np.array([0.,1e9,-1e9])
        for i in range(len(self.memberlists)):
            c=len(self.memberlists[i])
            # logging.debug(f'Cell {i} has {c} members.')
            # for j in range(len(self.memberlists[i])):
                # logging.debug(f'->{j}: {self.memberlists[i][j]}')
            amm[0]+=c
            if c<amm[1]:
                amm[1]=c 
            elif c>amm[2]:
                amm[2]=c
        assert amm[0]==cdf.shape[0]
        amm[0]/=len(self.memberlists)
        logging.debug(f'Linkcell.make_memberlists() ends. Max: {amm[2]:.0f} min {amm[1]:.0f} avg {amm[0]:0.3f}')

    def neighbors_of_cellndx(self,Ci):
        assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
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

    def are_cellndx_neighbors(self,Ci,Cj):
        assert self.cellndx_in_structure(Ci),f'Error: cell {Ci} outside of cell structure {self.ncells}'
        assert self.cellndx_in_structure(Cj),f'Error: cell {Cj} outside of cell structure {self.ncells}'
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

    def are_ldx_neighbors(self,ildx,jldx):
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
        assert n==6
    print(f'check2 {all(check2)}')
