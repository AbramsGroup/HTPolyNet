from HTPolyNet.coordinates import Coordinates

c=Coordinates.fromGroFile('init.gro')

print(c.DF)

c.delete_atoms(list(range(2000,2100)))

with open('test.gro','w') as f:
    f.write(str(c))
