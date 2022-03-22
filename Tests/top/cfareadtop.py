import parmed

s=parmed.gromacs.GromacsTopologyFile('init.itp')

print(s.atomtypes)