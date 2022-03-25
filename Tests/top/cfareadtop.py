from HTPolyNet.topology import Topology,typeorder,typedata


s1=Topology()

for _ in range(10):
    s2=Topology.from_topfile('STY.itp')
    s1.merge(s2)
    s2=Topology.from_topfile('VEA.itp')
    s1.merge(s2)

with open('test.top','w') as f:
    f.write(str(s1))

# a=('cc','ca','ca','ca')
# print(typeorder(a))