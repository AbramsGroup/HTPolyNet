from HTPolyNet.coordinates import Coordinates

c=Coordinates.fromMol2File('test.mol2')
print(c.DF)
if 'BDF' in c.__dict__:
    print(c.BDF)


