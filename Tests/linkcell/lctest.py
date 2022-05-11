from HTPolyNet.linkcell import Linkcell
from HTPolyNet.coordinates import Coordinates
import logging
logging.basicConfig(filename='cfa.log',encoding='utf-8',filemode='w',format='%(asctime)s %(message)s',level=logging.DEBUG)
c=Coordinates.read_gro('init.gro')
c.wrap_coords()
l=Linkcell()
l.create(0.5,c.box)
l.populate(c,ncpu=8)

for i,I in enumerate(l.cellndx):
    logging.debug(f'{i} {I} {l.ldx_of_cellndx(I)}')
