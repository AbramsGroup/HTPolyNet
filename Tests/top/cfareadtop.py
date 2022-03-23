from HTPolyNet.topology import Topology

s=Topology.from_topfile('init.itp')

for k,v in s.D.items():
    print(f'Directive {k} has {len(v)} stanza(s):')
    for vv in v:
        print(vv)
    print('---')

print()

with open('test.top','w') as f:
    f.write(str(s))
