"""

.. module:: test_rings
   :synopsis: tests ringss
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger=logging.getLogger(__name__)
from HTPolyNet.ring import *
import networkx as nx
import pandas as pd
import numpy as np
from HTPolyNet.matrix4 import Matrix4

def regular_polygon(nsides,a):
    dangle=2*np.pi/nsides
    v=[]
    angle=0.0
    for s in range(nsides):
        x=a*np.cos(angle)
        y=a*np.sin(angle)
        v.append(np.array([x,y,0.0]))
        angle+=dangle
    return np.array(v)

def hexagon(a):
    return regular_polygon(6,a)

def pentagon(a):
    return regular_polygon(5,a)

class TestRings(unittest.TestCase):
    def test_ring_treadmill(self):
        r=Ring([2,4,6,8,10,12])
        l=[x for x in r.treadmill()]
        self.assertEqual(len(l),len(r.idx)-1)
        self.assertEqual(l[0],[ 4, 6, 8,10,12, 2])
        self.assertEqual(l[1],[ 6, 8,10,12, 2, 4])
        self.assertEqual(l[2],[ 8,10,12, 2, 4, 6])
        self.assertEqual(l[3],[10,12, 2, 4, 6, 8])
        self.assertEqual(l[4],[12, 2, 4, 6, 8,10])
    def test_ring_eq(self):
        r=Ring([101,201,301,401,501,601])
        rr=[Ring(x) for x in r.treadmill()]
        self.assertTrue(all([r==x for x in rr]))
        q=Ring(r.idx[::-1])
        self.assertEqual(r,q)
    def test_ring_list_graph(self):
        g=nx.Graph()
        g.add_edge(1,2)
        g.add_edge(2,3)
        g.add_edge(3,4)
        g.add_edge(4,5)
        g.add_edge(5,1)
        g.add_edge(5,6)
        g.add_edge(6,7)
        g.add_edge(7,8)
        g.add_edge(8,4)
        L=RingList(g)
        self.assertEqual(len(L),2)
        self.assertTrue(Ring([2,3,4,5,1]) in L)
        self.assertTrue(Ring([4,5,6,7,8]) in L)
        self.assertTrue(Ring([5,4,3,2,1]) in L)
        self.assertTrue(Ring([8,7,6,5,4]) in L)

    def test_ring_list_basic(self):
        ll=[101,201,301,401,501,601]
        L=RingList([])
        for i in range(10):
            L.append(Ring(ll))
            ll=[x+1 for x in ll]
        self.assertEqual(len(L),10)
        self.assertTrue(Ring([101,201,301,401,501,601]) in L)
        self.assertFalse(Ring([11,12,13,14]) in L)

    def test_ring_injest_coordinates(self):
        df=pd.DataFrame({
            'globalIdx':[1,2,3,4,5],
            'posX':np.array([1,0,-1,-10,0]),
            'posY':np.array([0,1,0,-10,-1]),
            'posZ':np.array([0,0,0,0,0])
        })
        r=Ring([1,2,3,5])
        r.injest_coordinates(df)
        self.assertTrue(np.all(r.P[0]==np.array([1.,0.,0.])))
        self.assertEqual(r.planarity,1.0)

    def test_ring_pierce(self):
        # some true tests for rings of various sizes
        for nsides in range(3,8):
            # ring in x-y plane
            b=Ring(list(range(1,nsides+1)))
            P=regular_polygon(nsides,1.5)
            B=np.array([[0.,0.,1.],[0.,0.,-1.]])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertTrue(didit)

            # ring and bond rotated/translated quasi-randomly
            M=Matrix4()
            M.rotate_axis(45.0,np.array([1,-2,1.5])).translate(np.array([10.0,10.0,11.0]))
            P=np.array([M.transform(p) for p in P])
            B=np.array([M.transform(p) for p in B])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertTrue(didit)

        # false tests: both ends on same side of plane
        for nsides in range(3,8):
            b=Ring(list(range(1,nsides+1)))
            P=regular_polygon(nsides,1.5)
            B=np.array([[-1.,0.,1.],[0.,1.,1.]])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertFalse(didit)

            # ring and bond rotated/translated quasi-randomly
            M=Matrix4()
            M.rotate_axis(45.0,np.array([1,-2,1.5])).translate(np.array([10.0,10.0,11.0]))
            P=np.array([M.transform(p) for p in P])
            B=np.array([M.transform(p) for p in B])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertFalse(didit)

        # false tests: bond intersection point on plane not inside ring
        for nsides in range(3,8):
            b=Ring(list(range(1,nsides+1)))
            P=regular_polygon(nsides,1.5)
            B=np.array([[1.7,1.7,1.],[1.7,1.7,-1.]])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertFalse(didit)

            # ring and bond rotated/translated quasi-randomly
            M=Matrix4()
            M.rotate_axis(45.0,np.array([1,-2,1.5])).translate(np.array([10.0,10.0,11.0]))
            P=np.array([M.transform(p) for p in P])
            B=np.array([M.transform(p) for p in B])
            b.injest_coordinates(pd.DataFrame({
                'globalIdx':b.idx,
                'posX':np.array([x[0] for x in P]),
                'posY':np.array([x[1] for x in P]),
                'posZ':np.array([x[2] for x in P])
            }))
            didit,point=b.pierced_by(B)
            self.assertFalse(didit)
