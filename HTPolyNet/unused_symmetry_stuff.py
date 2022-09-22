"""

.. module:: unused_symmetry_stuff
   :synopsis: tools for automated detection of molecular symmetry; as the name implies, HTPolyNet does not use these routines currently.  In pre-release versions, we attempted to use these routines along with very high-T MD simulations of flexible molecules to identify symmetry-equivalent atoms using a distance-matrix approach.  Now since we force the user to explicitly declare symmetry-sets, these routines are defunct.  They didn't work reliably anyway, but they're kept here to stimulate future efforts to make this feature work.  I've also included a couple of methods used for enumeration of angles and dihedrals that have not proven necessary.
   
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

# defunct old funtions from topology module that I can't just quit
from HTPolyNet.topology import typeorder,repeat_check,idxorder,_GromacsTopologyDirectiveHeaders_
import pandas as pd
def add_enumerated_angles(self,newbonds,ignores=[],quiet=True):       
    at=self.D['atoms']
    ijk=self.D['angletypes'].set_index(['i','j','k']).sort_index()
    newangles=[]
    for b in newbonds:
        ''' new angles due to other neighbors of b[0] '''
        for ai in [i for i in self.bondlist.partners_of(b[0]) if (i!=b[1] and not i in ignores)]:
            aj=b[0]
            ak=b[1]
            it=at.iloc[ai-1].type
            jt=at.iloc[aj-1].type
            kt=at.iloc[ak-1].type
            idx=typeorder((it,jt,kt))
            i,j,k=idxorder((ai,aj,ak))
            repeat_check((i,j,k))
            if idx in ijk.index:
                angletype=ijk.loc[idx,'func']  # why no .values[0]
            else:
                if not quiet:
                    logger.warning(f'Angle type {idx} ({ai}-{aj}-{ak}) not found.')
                angletype=1
            h=_GromacsTopologyDirectiveHeaders_['angles']
            data=[i,j,k,angletype,pd.NA,pd.NA]
            assert len(h)==len(data), 'Error: not enough data for new angle?'
            angledict={k:[v] for k,v in zip(h,data)}
            self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
            newangles.append([i,j,k])
        ''' new angles due to other neighbors of b[1] '''
        for ak in [k for k in self.bondlist.partners_of(b[1]) if (k!=b[0] and not k in ignores)]:
            ai=b[0]
            aj=b[1]
            it=at.iloc[ai-1].type
            jt=at.iloc[aj-1].type
            kt=at.iloc[ak-1].type
            idx=typeorder((it,jt,kt))
            i,j,k=idxorder((ai,aj,ak))
            repeat_check((i,j,k))
            if idx in ijk.index:
                angletype=ijk.loc[idx,'func'] # why no .values[0]
            else:
                if not quiet:
                    logger.warning(f'Angle type {idx} ({ai}-{aj}-{ak}) not found.')
                angletype=1
            h=_GromacsTopologyDirectiveHeaders_['angles']
            data=[i,j,k,angletype,pd.NA,pd.NA]
            assert len(h)==len(data), 'Error: not enough data for new angle?'
            angledict={k:[v] for k,v in zip(h,data)}
            self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
            newangles.append([i,j,k])
    return newangles
    
def add_enumerated_dihedrals(self,newbonds,ignores=[],quiet=True):
    newdihedrals=[]
    newpairs=[]
    at=self.D['atoms']
    ijkl=self.D['dihedraltypes'].set_index(['i','j','k','l']).sort_index()

    ''' new proper dihedrals for which the new bond is the central j-k bond '''
    for b in newbonds:
        aj,ak=idxorder(b)
        for ai in [i for i in self.bondlist.partners_of(aj) if (i!=ak and not i in ignores)]:
            for al in [l for l in self.bondlist.partners_of(ak) if (l!=aj and not l in ignores)]:
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                kt=at.iloc[ak-1].type
                lt=at.iloc[al-1].type
                idx=typeorder((it,jt,kt,lt))
                i,j,k,l=idxorder((ai,aj,ak,al))
                repeat_check((i,j,k,l),msg=f'central {j}-{k}')
                if idx in ijkl.index:
                    dihedtype=ijkl.loc[idx,'func'].values[0] # why values[0]
                else:
                    if not quiet:
                        logger.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                    dihedtype=9
                h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                diheddict={k:[v] for k,v in zip(h,data)}
                self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                newdihedrals.append([i,j,k,l])
                ''' i-l is a new 1-4 pair '''
                h=_GromacsTopologyDirectiveHeaders_['pairs']
                data=[i,l,1]
                pairdict={k:[v] for k,v in zip(h,data)}
                self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                newpairs.append((i,l))
        ''' new proper dihedrals for which the new bond is the i-j or j-i bond '''
        for ai,aj in zip(b,reversed(b)):
            for ak in [k for k in self.bondlist.partners_of(aj) if (k!=ai and not k in ignores)]:
                for al in [l for l in self.bondlist.partners_of(ak) if (l!=aj and not l in ignores)]:
                    it=at.iloc[ai-1].type
                    jt=at.iloc[aj-1].type
                    kt=at.iloc[ak-1].type
                    lt=at.iloc[al-1].type
                    idx=typeorder((it,jt,kt,lt))
                    i,j,k,l=idxorder((ai,aj,ak,al))
                    repeat_check((i,j,k,l),msg=f'i-j neighbor of j-k {j}-{k}')
                    if idx in ijkl.index:
                        dihedtype=ijkl.loc[idx,'func'].values[0]
                        # dihedtype=ijkl.loc[idx,'func']
                    else:
                        if not quiet:
                            logger.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                        dihedtype=9
                    h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                    data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                    assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                    diheddict={k:[v] for k,v in zip(h,data)}
                    self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                    newdihedrals.append([i,j,k,l])
                    ''' i-l is a new 1-4 pair '''
                    h=_GromacsTopologyDirectiveHeaders_['pairs']
                    data=[i,l,1]
                    pairdict={k:[v] for k,v in zip(h,data)}
                    self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                    newpairs.append((i,l))

        ''' new proper dihedrals for which the new bond is the k-l or l-k bond '''
        for ak,al in zip(b,reversed(b)):
            for aj in [j for j in self.bondlist.partners_of(ak) if (j!=al and not j in ignores)]:
                for ai in [i for i in self.bondlist.partners_of(aj) if (i!=ak and not i in ignores)]:
                    it=at.iloc[ai-1].type
                    jt=at.iloc[aj-1].type
                    kt=at.iloc[ak-1].type
                    lt=at.iloc[al-1].type
                    idx=typeorder((it,jt,kt,lt))
                    i,j,k,l=idxorder((ai,aj,ak,al))
                    repeat_check((i,j,k,l),msg=f'k-l neighbor of j-k {j}-{k}')
                    if idx in ijkl.index:
                        dihedtype=ijkl.loc[idx,'func'].values[0]
                    else:
                        if not quiet:
                            logger.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                        dihedtype=9
                    h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                    data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                    assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                    diheddict={k:[v] for k,v in zip(h,data)}
                    self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                    newdihedrals.append([i,j,k,l])
                    h=_GromacsTopologyDirectiveHeaders_['pairs']
                    ''' i-l is a new 1-4 pair '''
                    data=[i,l,1]
                    pairdict={k:[v] for k,v in zip(h,data)}
                    self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                    newpairs.append((i,l))
    return newdihedrals,newpairs