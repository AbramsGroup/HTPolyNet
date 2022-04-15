from sympy import li
from HTPolyNet.projectfilesystem import Library
import logging

logging.basicConfig(filename='test.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)

l=Library.system(verbose=True)

l.checkout('FDE-p.mol2')
l.checkin('FDE-p.mol2')
l.checkout('DFA@N1-KENLAUSUCKS#C1-p.gro')
l.checkout('DFA@N1-KENLAUSUCKS#C1-p.gor')
