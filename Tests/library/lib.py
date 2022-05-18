#import importlib.resources
import HTPolyNet.projectfilesystem as pfs
import logging

logging.basicConfig(filename='test.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)

pfs.pfs_setup(mock=True)
pfs.checkout('molecules/parameterized/KENLAUSUCKS.mol2')
