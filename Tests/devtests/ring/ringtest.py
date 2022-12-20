from http.cookiejar import CookiePolicy
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.ring import Ring, Segment
import numpy as np
import pandas as pd

pbc=[1,1,1]

c=TopoCoord('DFA.top','DFA.gro')
d=TopoCoord('a_seg.top','a_seg.gro')
c.ring_detector()
c.label_ring_atoms()
c.make_ringlist()
R=Ring(c.Coordinates.ringlist[1][:,1:])
R.analyze()
#print(R)
#print(f'origin {R.O}')
#print(f'normal {R.n}')
S=Segment([R.O+0.25*R.n,R.O-0.25*R.n])
print(f'segint {R.segint(S)}')
d.Coordinates.A.loc[0,['posX','posY','posZ']]=S.P[0]
d.Coordinates.A.loc[1,['posX','posY','posZ']]=S.P[1]
d.Coordinates.A.loc[0,'cycle-idx']=0
d.Coordinates.A.loc[1,'cycle-idx']=0
c.merge(d)
#c.linkcell_initialize(0.5)

c.write_gro('pierced.gro')

B=c.Coordinates.box.diagonal().copy()
ntests=15
for i in range(ntests):
    L=np.random.random(3)*B
    c.translate(L)
    c.wrap_coords()
    c.linkcell_initialize(0.5,force_repopulate=True)
    c.make_ringdflist()

    # c.write_gro('pierced-wrapped.gro')

    i=c.Coordinates.A.shape[0]-1
    j=c.Coordinates.A.shape[0]

    Ri=c.get_R(i)
    #Rj=c.get_R(j)
    Rjp=c.get_R(j)
    print(f'R{i} {Ri}\nR{j}p {Rjp}\nBox {B}')
#    print(c.Coordinates.unwrap(Rjp,Ri,pbc))
    # Rij=c.Coordinates.mic(Ri-Rj,pbc)
    # rij=np.sqrt(Rij.dot(Rij))
    # generate the nearest periodic image of Rj to Ri
    # Rjp=Ri-Rij
    # return array of atom coordinates of ring pierced by this bond, if any
    C=c.Coordinates.ringpierce_testing(Ri,Rjp,pbc)
    if type(C)==pd.DataFrame:  # this is a ring
        print(C.to_string())
        print('PIERCE')
        # all this generate a special output file for inspection
        cidx=C['globalIdx'].to_list()#C[:,0].astype(int) # get globalIdx's
        idx=[i,j]
        idx.extend(cidx) # list of globalIdx's for this output
        sub=c.Coordinates.subcoords(c.Coordinates.A[c.Coordinates.A['globalIdx'].isin(idx)].copy())
        sub.write_gro(f'ring-{i}-{j}='+'-'.join([f'{x}' for x in cidx])+'.gro')
    c.translate(-1*L)