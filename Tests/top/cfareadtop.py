from HTPolyNet.topology import Topology,typeorder

s=Topology.from_topfile('init.itp')

for k,v in s.D.items():
    print(f'Directive {k} section has {len(v)} entries')
    print('---')

print()

with open('test.top','w') as f:
    f.write(str(s))

# a=('cc','ca','ca','ca')
# print(typeorder(a))