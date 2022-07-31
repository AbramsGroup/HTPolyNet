from HTPolyNet.molecule import Molecule
from HTPolyNet.topocoord import TopoCoord

m=Molecule()
m.TopoCoord=TopoCoord(grofilename='FDE.gro',topfilename='FDE.top')

m.flip_stereocenter(12)