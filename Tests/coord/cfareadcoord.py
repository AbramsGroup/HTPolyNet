from HTPolyNet.coordinates import Coordinates

c=Coordinates.read_mol2('test.mol2')

c.delete_atoms([4,5])

c.write_mol2('out.mol2')


