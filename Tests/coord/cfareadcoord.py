from HTPolyNet.coordinates import Coordinates

c=Coordinates.read_mol2('test.mol2')
d=Coordinates.read_mol2('VEA.mol2')
d.translate_coords([10.,0.,0.])
c.merge(d)
c.write_mol2('out.mol2')


