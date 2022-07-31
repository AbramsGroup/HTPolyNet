from HTPolyNet.coordinates import Coordinates

c=Coordinates()
print('A empty?',c.A.empty)
print('mol2_bonds empty?',c.mol2_bonds.empty)

c=Coordinates.read_mol2('FDE.mol2')
print('A empty?',c.A.empty)
print('mol2_bonds empty?',c.mol2_bonds.empty)
c.write_gro('test.gro')

d=Coordinates.read_gro('FDE-p.gro')
print('A empty?',d.A.empty)
print('mol2_bonds empty?',d.mol2_bonds.empty)
d.write_mol2('test.mol2')

d.merge(c)
d.write_mol2('test2.mol2')


