"""

.. module:: unused_symmetry_stuff
   :synopsis: tools for automated detection of molecular symmetry; as the name implies, HTPolyNet does not use these routines currently.  In pre-release versions, we attempted to use these routines along with very high-T MD simulations of flexible molecules to identify symmetry-equivalent atoms using a distance-matrix approach.  Now since we force the user to explicitly declare symmetry-sets, these routines are defunct.  They didn't work reliably anyway, but they're kept here to stimulate future efforts to make this feature work.
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import numpy as np
import os
from pytrr import GroTrrReader
import logging

logger=logging.getLogger(__name__)

def encluster(i,j,c):
    """encluster enforce objects i and j to have the same cluster index c[i] and c[j]

    :param i: an index
    :type i: int
    :param j: another index
    :type j: int
    :param c: array of cluster indices
    :type c: list
    :return: True if i and j already have same cluster index, False otherwise
    :rtype: bool
    """
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

def symm(d,thresh=0.1,outfile=None):
    """symm Builds and returns an atom-idx-ordered list of sea-cluster indexes.
        Any two atoms with the same sea-cluster-index are considered
        symmetry equivalent.

    :param d: interatomic distance matrix
    :type d: np.ndarray((N,N),float)
    :param thresh: threshold value of euclidean norm of difference between rank-ordered distance matrix columns, defaults to 0.1
    :type thresh: float, optional
    :param outfile: name of outpufile, defaults to None
    :type outfile: str, optional
    :return: symmetry set indices for each atom
    :rtype: list
    """
    na=d.shape[0]
    # initialize per-atom cluster id's to atom indexes (0-)
    cluster_ids=np.arange(na)
    # sort every column of distance matrix
    d=np.sort(d,axis=0)
    # make a list of columns
    c=np.array([d[:,i] for i in range(d.shape[1])])
    l=[]
    # compute magnitude of distance between columns for 
    # all unique pairs of columns, store as a list
    # of 3-tuples (i,j,distance)
    for i,ri in enumerate(c):
        for j,rj in enumerate(c):
            if i<j: 
                ij=ri-rj
                l.append((i,j,np.sqrt(ij.dot(ij))))

    # sort the list by distances (not strictly necessary)
    L=sorted(l,key=lambda x: x[2])
    if outfile:
        with open(outfile,'w') as f:
            for i,j,r in L:
                f.write(f'{i} {j} {r:0.8f}\n')

    # clusterize the i,j entries of each element of L
    all_done=False
    cpass=0
    while not all_done:
        all_done=True
        for i,j,r in L:
            if r<thresh:
                all_done=encluster(i,j,cluster_ids)
        cpass+=1

    # reindex cluster id 
    cpop={}
    for i,c in enumerate(cluster_ids):
        if not c in cpop:
            cpop[c]=[]
        cpop[c].append(i)
    rc=set(cpop.keys())
    for i,r in enumerate(rc):
        for j in cpop[r]:
            cluster_ids[j]=i
    return cluster_ids

def analyze_sea(deffnm,thresh=0.1):
    """analyze_sea Builds and returns an atom-idx-ordered list of sea-cluster indexes.
        Any two atoms with the same sea-cluster-index are considered
        symmetry equivalent.  The main job of this method is to compute the time-averaged interatomic
        distance matrix.  This matrix, if computed from a "hot" md simulation,
        should reveal atoms that are topologically symmetric, since the set of 
        average interatomic distances from atom A to all other atoms and the set
        of average interatomic distances from atom B to all other atoms are the 
        same if A and B are symmetry-equivalent.

    :param deffnm: gromacs mdrun deffnm
    :type deffnm: str
    :param thresh:  threshold value of euclidean norm of difference between rank-ordered distance matrix columns, defaults to 0.1
    :type thresh: float, optional
    """
    if not os.path.exists(f'{deffnm}.trr'):
        logger.error(f'{deffnm}.trr not found.')
        return []
    logger.debug(f'SEA analysis from {deffnm}.trr')
    with GroTrrReader(f'{deffnm}.trr') as trrfile:
        d=np.array((0,))
        nframes=0
        for mobydick in trrfile:
            frame_natoms=mobydick['natoms']
            if d.shape==(1,):
                d=np.zeros((frame_natoms,frame_natoms))
            data=trrfile.get_data()
            # tally all interatomic distances
            for i,ri in enumerate(data['x']):
                for j,rj in enumerate(data['x']):
                    rij=ri-rj
                    dist=np.sqrt(rij.dot(rij))
                    d[i][j]+=dist
            nframes+=1
        # averages over frames
        d/=nframes
        logger.debug(f'{deffnm}.trr: {nframes} frames')
        # send the distance matrix to be processed, return
        # the atom-ordered list of sea-cluster-idx's
        return symm(d,thresh=thresh,outfile=f'{deffnm}-symmanalysis.dat')
