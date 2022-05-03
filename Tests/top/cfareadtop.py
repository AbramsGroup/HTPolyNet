from HTPolyNet.topology import Topology,typeorder,typedata
 
s1=Topology.read_gro('DFAFDE.top')
atdf=s1.D['atoms']

t=s1.bondtree_as_list((1,43),4)
mappables=[]
for b in t:
    for i in range(2):
        if not b[i] in mappables:
            mappables.append(b[i])
madf=atdf[atdf['nr'].isin(mappables)]
print(f'Bond {b} mappable atoms:\n{madf.to_string()}')