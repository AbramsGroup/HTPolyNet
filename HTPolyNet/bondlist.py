''' a simple class for constructing a two-way bondlist dictionary '''
import numpy as np
import networkx as nx

class Bondlist:
    def __init__(self):
        self.B={}
    
    @classmethod
    def fromDataFrame(cls,df):
        inst=cls()
        inst.update(df)
        return inst

    def update(self,df):
        if not 'ai' in df.columns and not 'aj' in df.columns:
            raise Exception('Bondlist expects a dataframe with columns "ai" and "aj".')
        aiset=set(df.ai)
        ajset=set(df.aj)
        keyset=aiset.union(ajset)
        keys=sorted(list(keyset))
        self.B.update({k:[] for k in keys})
        for i,row in df.iterrows():
            self.B[row.loc['ai']].append(row.loc['aj'])
            self.B[row.loc['aj']].append(row.loc['ai'])

    def __str__(self):
        retstr=''
        for k,v in self.B.items():
            retstr+=f'{k}: '+' '.join(str(vv) for vv in v)+'\n'
        return retstr
    
    def partners_of(self,idx):
        if idx in self.B:
            return self.B[idx]
        return []

    def are_bonded(self,idx,jdx):
        if idx in self.B and jdx in self.B:
            return jdx in self.B[idx]
        return False

    def append(self,pair=[]):
        if len(pair)==2:
            ai,aj=min(pair),max(pair)
            if not ai in self.B:
                self.B[ai]=[]
            self.B[ai].append(aj)
            if not aj in self.B:
                self.B[aj]=[]
            self.B[aj].append(ai)

    def delete_atoms(self,idx=[]):
        ''' delete entries '''
        for i in idx:
            if i in self.B:
                del self.B[i]
        ''' delete members in other entries '''
        for k,v in self.B.items():
            for i in idx:
                if i in v:
                    self.B[k].remove(i)
    
    def adjacency_matrix(self):
        N=len(self.B)
        A=np.zeros((N,N)).astype(int)
        for i,n in self.B.items():
            for j in n:
                A[i-1,j-1]=1
                A[j-1,i-1]=1
        return A

    def graph(self):
        g=nx.DiGraph()
        N=len(self.B)
        g.add_nodes_from(list(range(N)))
        for i,n in self.B.items():
            for j in n:
                g.add_edge(i-1,j-1)
        return g
