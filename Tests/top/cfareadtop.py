from HTPolyNet.topology import Topology,typeorder,typedata
 
s1=Topology.read_gro('DFA.top')
a=s1.D['atoms']
p=s1.D['pairs']
np=p.shape[0]
pai=p.ai.to_list()
paj=p.aj.to_list()
d=s1.D['dihedrals']
nd=d.shape[0]
print(f'{np} pairs {nd} dihedrals')
dd=[]
for pi,pj in zip(pai,paj):
    dwpi=d[((d.ai==pi)&(d.al==pj))|((d.ai==pj)&(d.al==pi))].index.to_list()
    dd.extend(dwpi)
print(f'{d.loc[dd].shape[0]} dihedrals with i-l in [ pairs ]')
d=d.loc[dd]
ai=d.ai
al=d.al
print(a[a.nr.isin(ai)|a.nr.isin(al)]['atom'])

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