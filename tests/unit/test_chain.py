"""

.. module:: test_chain
   :synopsis: tests chain
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import pytest
import unittest
import logging
logger=logging.getLogger(__name__)
from HTPolyNet.chain import *
import networkx as nx
import pandas as pd
import numpy as np

class TestChain(unittest.TestCase):
    def testchain_create(self):
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        self.assertEqual(len(cm.chains),2)

    def testchain_nocreate(self):
        cm=ChainManager()
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        self.assertEqual(len(cm.chains),0)

    def testchain_findatom(self):
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        c1=cm.chain_of(1)
        c2=cm.chain_of(2)
        c3=cm.chain_of(3)
        c4=cm.chain_of(4)
        self.assertEqual(c1,c2)
        self.assertEqual(c3,c4)
        self.assertNotEqual(c1,c3)
    
    def testchain_merge(self):
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        self.assertEqual(len(cm.chains),2)        
        cm.injest_bond(2,3)
        self.assertEqual(len(cm.chains),1)
        c=cm.chain_of(1)
        self.assertTrue(c.is_head(1))
        self.assertTrue(c.is_tail(4))

    def testchain_cycle(self):
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        cm.injest_bond(2,3)
        self.assertEqual(len(cm.chains),1)        
        cm.injest_bond(4,1)
        self.assertTrue(cm.chains[0].is_cyclic)
        c=cm.chain_of(1)
        self.assertFalse(c.is_head(1))
        self.assertFalse(c.is_tail(4))
    
    def testchain_cycle(self):
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        cm.injest_bond(2,3)
        self.assertEqual(len(cm.chains),1)        
        cm.injest_bond(4,1)
        self.assertTrue(cm.chains[0].is_cyclic)

    def testchain_error1(self):
        with pytest.raises(Exception) as e_info:
            cm=ChainManager(create_if_missing=True)
            cm.injest_bond(1,2)
            cm.injest_bond(2,3)
            # atom 3 is not in a chain
        self.assertEqual(e_info.value.args[0],'This is a bug - no j-chain!')
    
    def testchain_error2(self):
        with pytest.raises(Exception) as e_info:
            cm=ChainManager(create_if_missing=True)
            cm.injest_bond(1,2)
            cm.injest_bond(5,1)
            # atom 5 is not a chain
        self.assertEqual(e_info.value.args[0],'This is a bug - no i-chain!')

    def testchain_todataframe(self):
        df=pd.DataFrame({'globalIdx':[1,2,3,4,5,6],'bondchain':[-1]*6,'bondchain_idx':[-1]*6})
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        cm.injest_bond(5,6)
        cm.to_dataframe(df)
        df.set_index('globalIdx',inplace=True)
        self.assertEqual(df.loc[1,'bondchain'],0)
        self.assertEqual(df.loc[1,'bondchain_idx'],0)
        self.assertEqual(df.loc[2,'bondchain'],0)
        self.assertEqual(df.loc[2,'bondchain_idx'],1)
        self.assertEqual(df.loc[3,'bondchain'],1)
        self.assertEqual(df.loc[3,'bondchain_idx'],0)
        self.assertEqual(df.loc[4,'bondchain'],1)
        self.assertEqual(df.loc[4,'bondchain_idx'],1)
        self.assertEqual(df.loc[5,'bondchain'],2)
        self.assertEqual(df.loc[5,'bondchain_idx'],0)
        self.assertEqual(df.loc[6,'bondchain'],2)
        self.assertEqual(df.loc[6,'bondchain_idx'],1)
                
    def testchain_todataframe_merge(self):
        df=pd.DataFrame({'globalIdx':[1,2,3,4,5,6],'bondchain':[-1]*6,'bondchain_idx':[-1]*6})
        cm=ChainManager(create_if_missing=True)
        cm.injest_bond(1,2)
        cm.injest_bond(3,4)
        cm.injest_bond(5,6)
        cm.injest_bond(2,3)
        cm.injest_bond(4,5)
        cm.to_dataframe(df)
        df.set_index('globalIdx',inplace=True)
        self.assertEqual(df.loc[1,'bondchain'],0)
        self.assertEqual(df.loc[1,'bondchain_idx'],0)
        self.assertEqual(df.loc[2,'bondchain'],0)
        self.assertEqual(df.loc[2,'bondchain_idx'],1)
        self.assertEqual(df.loc[3,'bondchain'],0)
        self.assertEqual(df.loc[3,'bondchain_idx'],2)
        self.assertEqual(df.loc[4,'bondchain'],0)
        self.assertEqual(df.loc[4,'bondchain_idx'],3)
        self.assertEqual(df.loc[5,'bondchain'],0)
        self.assertEqual(df.loc[5,'bondchain_idx'],4)
        self.assertEqual(df.loc[6,'bondchain'],0)
        self.assertEqual(df.loc[6,'bondchain_idx'],5)

    def testchain_merge(self):
        cm1=ChainManager(create_if_missing=True)
        cm2=ChainManager(create_if_missing=True)
        cm1.injest_bond(1,2)
        cm1.injest_bond(3,4)
        cm1.injest_bond(2,3)
        cm2.injest_bond(11,12)
        cm2.injest_bond(13,14)
        cm2.injest_bond(12,13)
        cm1.injest_other(cm2)
        self.assertEqual(cm1.chain_of(1).idx,0)
        self.assertEqual(cm1.chain_of(11).idx,1)

    def testchain_from_dataframe(self):
        df=pd.DataFrame({'i':[1,2,3,4,5,6],'bondchain':[0,0,1,1,2,2],'bondchain_idx':[0,1,0,1,0,1]})
        cm=ChainManager()
        cm.from_dataframe(df)
        self.assertEqual(len(cm.chains),3)
        self.assertEqual(cm.chain_of(1).idx,0)
        self.assertEqual(cm.chain_of(3).idx,1)
        self.assertEqual(cm.chain_of(5).idx,2)

    def testchain_from_dataframe_reindex(self):
        df=pd.DataFrame({'i':[1,2,3,4,5,6],'bondchain':[10,10,101,101,2003,2003],'bondchain_idx':[0,1,0,1,0,1]})
        cm=ChainManager()
        cm.from_dataframe(df)
        self.assertEqual(len(cm.chains),3)
        self.assertEqual(cm.chain_of(1).idx,0)
        self.assertEqual(cm.chain_of(3).idx,1)
        self.assertEqual(cm.chain_of(5).idx,2)

    def testchain_from_dataframe_nochains(self):
        df=pd.DataFrame({'i':[1,2,3,4,5,6],'bondchain':[-1]*6,'bondchain_idx':[-1]*6})
        cm=ChainManager()
        cm.from_dataframe(df)
        self.assertEqual(len(cm.chains),0)
