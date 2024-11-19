"""

.. module:: ring
   :synopsis: handles ring-piercing determinations
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
# Pierced rings
# Cameron F. Abrams cfa22@drexel.edu
#
# How to use (suggested):
#
# 1. Create a list of rings from a coordinate snapshot:
#    Suppose X is an Nx3 numpy array ordered such
#    that each consecutive group of six elements are
#    the positions of carbon atoms of a phenyl ring.
#    Then
#
#    Rings=[]
#    for i in range(0,len(X),6):
#        Rings.append(Ring(np.array([X[i+j] for j in range(6)]))
#
#    This can be done once per CURE iteration as long as R is
#    in scope.
#
# 2. Now, suppose B is a 2x3 numpy array containing coordinates
#    of two atoms that are a potential bond.  Cast this as a 
#    Segment, and then loop over Rings until a piercing is found
#   
#    S=Segment(B)
#    for R in Rings:
#        pierced,P=R.segint(S)
#        if pierced:
#            # print some message
#            # set some "not allowed to bond" flag
#            break
#    
# IMPORTANT NOTE: no MIC is used here, so all 
# coordinates must be unwrapped into the *same*
# periodic image!
#
import numpy as np
from collections import UserList
from functools import singledispatchmethod
import networkx as nx
from copy import deepcopy
import logging
logger=logging.getLogger(__name__)

def lawofcos(a,b):
    """lawofcos return the cosine of the angle defined by vectors a and b if they share a vertex (the LAW OF COSINES)

    :param a: a vector
    :type a: numpy.ndarray(3,float)
    :param b: another vector
    :type b: numpy.ndarray(3,float)
    :return: cosine of the angle formed by a and b
    :rtype: float
    """
    return np.dot(a,b)/np.sqrt(np.dot(a,a)*np.dot(b,b))

class Segment:
    """ a segment object owns a list of Points P with two elements representing segment endpoints, and a vector that points from the first point to the second, V
    """
    def __init__(self,P):
        """__init__ generates a new Segment object from the points in container P

        :param P: listlike container of two points, each of which is a 3-dimensional numpy array
        :type P: list
        """
        self.P=P.copy()
        # will need to recompute this when molecules are shifted
        self.V=self.P[1]-self.P[0] # p1=p0+t*(p1-p0)

class Ring:
    def __init__(self,idx):
        """__init__ generates a Ring object from the list of atom globalIdx

        A ring is a sequence of integers that is treated as cyclic and bidirectional.

        So, the list [1,2,3,4,5] is "equal" to the following lists:

        [2,3,4,5,1]
        [3,4,5,1,2]
        [4,5,1,2,3]
        [5,1,2,3,4]
        [5,4,3,2,1]
        [4,3,2,1,5]
        [3,2,1,5,4]
        [2,1,5,4,3]
        [1,5,4,3,2]

        The first four elements in the list above are "treadmilled" versions of the
        parent list.  The final five elements are the reverse of the parent list
        and all treadmilled version of that reversed list.

        :param P: list of ints
        :type P: list
        """
        self.idx=idx.copy()
        assert(all([type(x)==int for x in self.idx]))

    def copy(self):
        newring=deepcopy(self)
        return newring

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        self.P=np.array(A[A[idx_key].isin(self.idx)][pos_key].values)
        # logger.debug(f'P {self.P}')
        self.O=np.mean(self.P,axis=0)
        # logger.debug(f'O {self.O}')
        iR=Ring(list(range(len(self.idx))))
        a=iR.treadmill()
        self.B=[]
        for i,j in zip(iR.idx,next(a)):
            b=self.P[i]-self.P[j]
            # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {b}')
            self.B.append(b)
            # logger.debug(f'R-{self.idx[i]}-{self.idx[j]}: {self.B[-1]}')
        self.B=np.array(self.B)
        # logger.debug(f'B {self.B}')
        # Approximate the normal unit vector of the plane of the ring
        # by averaging cross-products of all adjacent origin-to-vertex
        # vectors.  This orients the unit vector such that it is 
        # in the +z direction in the local coordinate system defined
        # by the ring with vertices ordered in the counter-clockwise
        # direction in the plane.
        a=iR.treadmill()
        self.C=[]
        for i,j in zip(iR.idx,next(a)):
            c=np.cross(self.B[i],self.B[j])
            # logger.debug(f'C-{self.idx[i]}-{self.idx[j]}: {c}')
            self.C.append(c)
        self.C=np.array(self.C)
        # logger.debug(f'C {self.C}')
        n=np.sum(self.C,axis=0)
        self.n=n/np.linalg.norm(n)
        # logger.debug(f'n {self.n}')
        # compute planarity as average of all
        # cross-i-cross-i+1 dot products
        a=iR.treadmill()
        self.planarity=0
        for i,j in zip(iR.idx,next(a)):
            ci=self.C[i]
            cj=self.C[j]
            self.planarity+=lawofcos(ci,cj)
        self.planarity/=len(iR.idx)
        # get the d for n-dot-r + d = 0 equation of the plane
        # n[0]*(x-O[0])+n[1]*(y-O[1])+n[2]*(z-O[2])=0
        # n[0]*x + n[1]*y + n[2]*z - (n[0]*O[0]+n[1]*O[1]+n[2]*O[2]) = 0
        self.d=-np.dot(self.n,self.O)        
        # self.vP=[]
        # for v in self.P:
        #     r=v-self.O
        #     p=r-np.dot(r,self.n)*self.n
        #     newv=p+self.O
        #     self.vP.append(newv)
        # self.vP=np.array(self.vP)
        # logger.debug(f'vP {self.vP}')

    def treadmill(self):
        """ yield the treadmilled versions of the list """
        for i in range(1,len(self.idx)):
            yield self.idx[i:]+self.idx[:i]

    def __eq__(self,other):
        check1=any([self.idx==other.idx] + [self.idx==x for x in other.treadmill()])
        check2=any([self.idx==other.idx[::-1]]+[self.idx==x[::-1] for x in other.treadmill()])
        # logger.debug(f'checking ring eq: {str(self)}=={str(other)}? {check1} and {check2}')
        return check1 or check2
    
    def __str__(self):
        return '-'.join([str(x) for x in self.idx])

    def shift(self,shift):
        self.idx=[x+shift for x in self.idx]
        return self

    def remap(self,mapper):
        new_idx=[mapper[x] for x in self.idx]
        self.idx=new_idx

    def unwrap(self,P,unwrapf=None,pbc=[1,1,1]):
        r=self.copy()
        for i in range(len(r.P)):
            r.P[i]=unwrapf(r.P[i],P,pbc=pbc)
        return r

    def pierced_by(self,P):
        """determines if segment with endpoints P[0] and P[1] pierces ring; 
        uses ray projection method and fact that scaled length must be between 
        0 and 1 for a plane intersection

        :param P: a 2-element numpy array of 3-space points
        :type S: numpy.ndarray
        :return: True if P[0]-P[1] pierces self's ring, along with the intersection point
        :rtype: tuple (boolean, numpy.ndarray(3))
        """
        V=P[0]-P[1]
        # If these two points are on either side of the ring plane,
        # then the length of the normal vector from P[0] to the plane
        # will be less than the length of the vector from P[1] to P[0] projected
        # onto the plane normal. 
        denom=np.linalg.norm(np.dot(self.n,V)*self.n) # length of vector from P[1] to P[0] projected onto normal
        OP=np.dot(self.n,P[0]-self.O)*self.n # normal vector from plane to P[0]
        numer=np.linalg.norm(OP) # length of normal vector from plane to P[0]
        # logger.debug(f'numer {numer} denom {denom}')
        t=-1.0
        if denom>1.e-13:
            t=numer/denom # a fraction if P[0] and P[1] are on either side of the plane
        # logger.debug(f't {t}')
        if 0<t<1:
            # compute point in ring plane that marks intersection with this vector
            # V=P[0]-P[1]
            # P[1]=P[0]-V
            # intersection point = P[0]-t*V
            PP=P[0]-t*V
            # logger.debug(f'PP {PP}')
            # determine if PP is inside ring:
            # compute the series of unit-vector cross-products v(i)-PP-v((i+1)%N)
            # sum will have a large component along normal vector if yes, essentially 0 if no
            iR=Ring(list(range(len(self.idx))))
            a=iR.treadmill()
            sumC=np.zeros(3)
            for i,j in zip(iR.idx,next(a)):
                V1=PP-self.P[i]
                V2=PP-self.P[j]
                c=np.cross(V1/np.linalg.norm(V1),V2/np.linalg.norm(V2))
                sumC+=c/np.linalg.norm(c)
            tst=np.dot(self.n,sumC/len(iR.idx))
            is_inside=np.all(np.isclose([tst],[1.0]))
            # logger.debug(f'tst {tst} is_inside {is_inside}')
            return is_inside,PP
        return False,np.ones(3)*np.nan

class RingList(UserList):
    @singledispatchmethod
    def __init__(self,input_obj):
        self.data=input_obj
    @__init__.register(nx.Graph)
    def _from_graph(self,G):
        L=[]
        for ll in nx.chordless_cycles(G):
            L.append(Ring(ll))
        super().__init__(L)

    def shift(self,shift):
        for item in self:
            item.shift(shift)
        return self
    
    def all_atoms(self):
        retlist=[]
        for item in self:
            retlist.extend(item.idx)
        return retlist

    def injest_coordinates(self,A,idx_key='globalIdx',pos_key=['posX','posY','posZ']):
        for item in self:
            item.injest_coordinates(A,idx_key=idx_key,pos_key=pos_key)
    
    def filter(self,idxlist):
        retL=RingList([])
        for item in self:
            if any([x in idxlist for x in item.idx]):
                retL.append(item)
        return retL

    def remap(self,mapper):
        for item in self:
            item.remap(mapper)

    def __str__(self):
        return ';'.join([str(x) for x in self])