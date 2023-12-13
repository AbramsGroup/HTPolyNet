"""

.. module:: test_dataframetools
   :synopsis: tests dataframetools
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import os
import logging
logger=logging.getLogger(__name__)
from HTPolyNet.dataframetools import *
import pandas as pd

class TestDataframeTools(unittest.TestCase):
    def test_getrow(self):
        df=pd.DataFrame({
            'a':[ 1, 2, 3, 4, 5],
            'b':[ 6, 7, 8, 9,10],
            'c':[11,12,13,14,15]
        }
        )
        qdict={'a':3,'b':8}
        row=get_row(df,qdict)
        ans=pd.Series({'a':3,'b':8,'c':13})
        res=row==ans
        self.assertTrue(res.all())
