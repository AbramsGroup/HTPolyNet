from HTPolyNet.coordinates import Coordinates
from HTPolyNet.configuration import CappingBond
from HTPolyNet.ambertools import GAFFParameterize
from HTPolyNet.libraries import Library
from HTPolyNet.gromacs import grompp_and_mdrun

l=Library()
print('Library')
print(l)
c=Coordinates.read_mol2('FDE.mol2')
cappers=[CappingBond({"pair":["O1","C1"], "bondorder":1}),
CappingBond({"pair":["O2","C2"], "bondorder":1})]
c.cap(cappers)

GAFFParameterize('FDE-capped-unminimized','FDE-capped-unminimized-parameterized',force=False,parmed_save_inline=False)
pc=Coordinates.read_mol2('FDE-capped-unminimized-parameterized.mol2')
l.fetch('em-single-molecule.mdp')
grompp_and_mdrun(gro='FDE-capped-unminimized-parameterized',
                 top='FDE-capped-unminimized-parameterized',
                 mdp='em-single-molecule',
                 out='FDE-capped-minimized',boxSize=[15,15,15])
mc=Coordinates.read_gro('FDE-capped-minimized.gro')
pc.copy_coords(mc)
pc.write_mol2('FDE-capped-minimized.mol2')
GAFFParameterize('FDE-capped-minimized','FDE-capped-minimized-parameterized',force=False,parmed_save_inline=False)