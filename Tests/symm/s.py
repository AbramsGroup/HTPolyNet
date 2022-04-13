from pytrr import GroTrrReader
import argparse as ap
import numpy as np
# idea: from a long NVT MD, use mean interatomic distance matrix
# to determine atom symmetry classes
#
# seems to work!  need to implement as an htpolynet command to
# allow for standalone symmetry class detection
#
import matplotlib.pyplot as plt

def encluster(i,j,c):
    if c[i]==c[j]:
        return True
    if c[j]>c[i]:
        for k in range(len(c)):
            if c[k]==c[j]:
                c[k]=c[i]
    elif c[j]<c[i]:
        for k in range(len(c)):
            if c[k]==c[i]:
                c[k]=c[j]
    return False        

def symm(d,thresh=0.10):
    na=d.shape[0]
    # print(na)
    cluster_ids=np.arange(na)
    # print(cluster_ids)
    d=np.sort(d,axis=0)
    c=np.array([d[:,i] for i in range(d.shape[1])])
#    print(c[15],c[16],c[15]-c[16],(c[15]-c[16]).dot(c[15]-c[16]))
    l=[]
    for i,ri in enumerate(c):
        for j,rj in enumerate(c):
            if i<j: 
                ij=ri-rj
                l.append((i,j,np.sqrt(ij.dot(ij))))
    L=sorted(l,key=lambda x: x[2])
    # sixteen=[l for l in L if l[0]==16 or l[1]==16]
    # print(sixteen)
    lim=40
    x=np.arange(lim)
    fig,ax=plt.subplots(1,1,figsize=(5,15))
    ax.scatter(x,[j[2] for j in L[:lim]])
    for X,Y,i,j in zip(x,[j[2] for j in L[:lim]],[j[0] for j in L[:lim]],[j[1] for j in L[:lim]]):
        ax.text(X,Y,f'{i}-{j}')
    plt.savefig('plt.png')
        
    all_done=False
    cpass=0
    while not all_done:
        all_done=True
        for i,j,r in L:
            if r<thresh:
                all_done=encluster(i,j,cluster_ids)
        cpass+=1
    cpop={}
    for i,c in enumerate(cluster_ids):
        if not c in cpop:
            cpop[c]=[]
        cpop[c].append(i)
    rc=set(cpop.keys())
    dpop={}
    for i,r in enumerate(rc):
        dpop[i]=cpop[r]
        for j in dpop[i]:
            cluster_ids[j]=i
    return dpop,cluster_ids

parser=ap.ArgumentParser()
parser.add_argument('-t',type=float,default=0.2)
args=parser.parse_args()

with GroTrrReader('nvt1000.trr') as trrfile:
    # print(trrfile.filename)
    d=np.array((0,))
    nframes=0
    for mobydick in trrfile:
        frame_natoms=mobydick['natoms']
        if d.shape==(1,):
            print(frame_natoms)
            d=np.zeros((frame_natoms,frame_natoms))
            natoms=frame_natoms
        data=trrfile.get_data()
        for i,ri in enumerate(data['x']):
            for j,rj in enumerate(data['x']):
                rij=ri-rj
                dist=np.sqrt(rij.dot(rij))
                d[i][j]+=dist
        nframes+=1
    print(nframes)
    d/=nframes
    c,cids=symm(d,thresh=args.t)
    for i,v in c.items():
        print(i,v)
