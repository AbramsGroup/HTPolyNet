#from HTPolyNet.topology import Topology,typeorder,typedata
from HTPolyNet.topocoord import TopoCoord

s1=TopoCoord(topfilename='DFA.top',grofilename='DFA.gro')

a,b,c=s1.get_gro_attribute_by_attributes(['atomName','resNum','resName'],{'globalIdx':3})
print(a,b,c)
#print(s1.D['dihedrals'].set_index(['ai','aj','ak','al']).to_string())



# atdf=s1.D['atoms']

# t=s1.bondtree_as_list((1,43),4)
# mappables=[]
# for b in t:
#     for i in range(2):
#         if not b[i] in mappables:
#             mappables.append(b[i])
# madf=atdf[atdf['nr'].isin(mappables)]
# print(f'Bond {b} mappable atoms:\n{madf.to_string()}')