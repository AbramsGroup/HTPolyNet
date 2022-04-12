from pytrr import GroTrrReader
import numpy as np
# idea: from a long NVT MD, use mean interatomic distance matrix
# to determine atom symmetry classes
#
# seems to work!  need to implement as an htpolynet command to
# allow for standalone symmetry class detection
#
import matplotlib.pyplot as plt

def symm(d):
    d=np.sort(d,axis=0)
    c=np.array([d[:,i] for i in range(d.shape[1])])
    l=[]
    for i,ri in enumerate(c):
        for j,rj in enumerate(c):
            ij=ri-rj
            if i<j: 
                l.append((i,j,np.sqrt(ij.dot(ij))))
    L=sorted(l,key=lambda x: x[2])
    # TODO: must analyze this to derive clusters/classes and assign
    # class id to each atom
    return L

with GroTrrReader('nvt.trr') as trrfile:
    # print(trrfile.filename)
    d=np.array((0,))
    print(d.shape)
    for mobydick in trrfile:
        frame_natoms=mobydick['natoms']
        if d.shape==(1,):
            d=np.zeros((frame_natoms,frame_natoms))
            natoms=frame_natoms
        data=trrfile.get_data()
        for i,ri in enumerate(data['x']):
            for j,rj in enumerate(data['x']):
                rij=ri-rj
                dist=np.sqrt(rij.dot(rij))
                d[i][j]+=dist
    d/=natoms
    L=symm(d)
    x=np.arange(100)
    fig,ax=plt.subplots(1,1)
    ax.scatter(x,[j[2] for j in L[:100]])
    plt.savefig('plt.png')
        