import HTPolyNet.processTop as pt

r=pt.cfaReadTop('FDE-un.top')

for k,d in r.items():
    print(k)
    print(d)