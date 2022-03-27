from HTPolyNet.topology import Topology,typeorder,typedata


s1=Topology()

for _ in range(10):
    s2=Topology.from_topfile('STY.itp')
    s1.merge(s2)
    s2=Topology.from_topfile('VEA.itp')
    s1.merge(s2)

print(f'Total charge: {s1.total_charge():.6e}')
# print('Bondlist:')
# print(str(s1.bondlist))
print(s1.get_atom(8))
print(s1.get_atom(26))
s1.add_bonds([(8,26)])

with open('test.top','w') as f:
    f.write(str(s1))

# a=('cc','ca','ca','ca')
# print(typeorder(a))