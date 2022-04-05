from HTPolyNet.coordinates import Coordinates
from HTPolyNet.configuration import CappingBond
from HTPolyNet.ambertools import GAFFParameterize

c=Coordinates.read_mol2('FDE.mol2')
cappers=[CappingBond({"pair":["O1","C1"], "bondorder":1}),
CappingBond({"pair":["O2","C2"], "bondorder":1})]
c.cap(cappers)
GAFFParameterize('FDE-capped-unminimized','FDE-capped-unminimized-parameterized',force=False,parmed_save_inline=False)
GAFFParameterize('FDE-capped-minimized','FDE-capped-parameterized',force=False,parmed_save_inline=False)