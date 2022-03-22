import HTPolyNet.processTop as pt

r=pt.GromacsTopToDataFrameDict('testDBEGA.itp')

for k,d in r.items():
    print(k)
    print(d)