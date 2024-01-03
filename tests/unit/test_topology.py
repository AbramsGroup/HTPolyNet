"""

.. module:: test_topology
   :synopsis: tests topology
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger=logging.getLogger(__name__)
import HTPolyNet.topology as tp
import os
import pandas as pd
import numpy as np

class TestTopology(unittest.TestCase):
    def test_read_top(self):
        fn='test.top' # this top file has blank entries for all bond, angle, and dihedral parameters
        T=tp.Topology.read_top(fn) # blanks are padded
        self.assertTrue('bonds' in T.D)
        self.assertEqual(tp._PAD_,-99.99)
        self.assertTrue(all(T.D['bonds']['c1']==tp._PAD_))
        self.assertTrue(all(T.D['bonds']['c0']==tp._PAD_))
        self.assertTrue(all(T.D['angles']['c1']==tp._PAD_))
        self.assertTrue(all(T.D['angles']['c0']==tp._PAD_))
        for c in [f'c{i}' for i in range(6)]:
            self.assertTrue(all(T.D['dihedrals'][c]==tp._PAD_))
    def test_write_top(self):
        if os.path.exists('write_test.top'):
            os.remove('write_test.top')
        fn='test.top'
        T=tp.Topology.read_top(fn)
        T.write_top('write_test.top')
        self.assertTrue(os.path.exists('write_test.top'))
        W=tp.Topology.read_top('write_test.top',pad=pd.NA)
        os.remove('write_test.top')
        self.assertTrue(all(W.D['bonds']['c1'].isna()))
        self.assertTrue(all(W.D['bonds']['c0'].isna()))
        self.assertTrue(all(W.D['angles']['c1'].isna()))
        self.assertTrue(all(W.D['angles']['c0'].isna()))
        for c in [f'c{i}' for i in range(6)]:
            self.assertTrue(all(W.D['dihedrals'][c].isna()))
        self.assertTrue(W.D['atoms'].shape==T.D['atoms'].shape)

