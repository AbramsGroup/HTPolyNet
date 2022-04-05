from HTPolyNet.configuration import Monomer
from itertools import combinations_with_replacement, product

class OligoNode:
    def __init__(self):
        pass

class OligoTree:
    def __init__(self,baseconn):
        self.baseconn=baseconn
    
def get_conn(mol):
    conn=[]
    for mrname,mr in mol.reactive_atoms.items():
        ''' check to see if this atom's symmetry partners are already on the list '''
        donotadd=False
        for s in mr.sym:
            if s in conn:
                donotadd=True
        if not donotadd:
            conn.append(mrname)
    return conn

def nextconn(partners,conns):
    hic=conns[-1]
    for i in range(len(partners)):
        if partners[i]!=hic:
            break
    else:
        return False
    partners[i]=conns[0] if partners[i]==[] else conns[conns.index(partners[i])+1]
    print('inside',partners)
    return True

def react_mol2(m,n,minimize=True):
    ''' m and n are Monomer instances (defined in configuration.py) '''
    if not isinstance(m,Monomer) or not isinstance(n,Monomer):
        raise Exception(f'react_mol2 needs Monomers not {type(m)} and {type(n)}')
    if not 'active' in m.Topology or not 'active' in n.Topology:
        raise Exception('react_mol2 needs Monomers with active topologies')
    if minimize and (not 'active' in m.Coords or not 'active' in n.Coords):
        raise Exception('react_mol2 needs Monomers with active coordinates')

    print(f'react_mol2: {m.name} and {n.name}')
    for basemol,othermol in zip([m,n],[n,m]):
        basera=get_conn(basemol)
        baseconn=[basemol.reactive_atoms[i].z for i in basera]
        otherra=get_conn(othermol)
        otherconn=[othermol.reactive_atoms[i].z for i in otherra]
        otherra=['']+otherra
        basearr=[]
        for b,z in zip(basera,baseconn):
            basearr.append(list(combinations_with_replacement(otherra,z)))
        print(basearr)
        print(list(product(*basearr)))

    return []