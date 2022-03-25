''' a simple class for constructing a two-way bondlist dictionary '''

class Bondlist:
    def __init__(self):
        self.B={}
    
    @classmethod
    def fromDataFrame(cls,df):
        if not 'ai' in df.columns and not 'aj' in df.columns:
            raise Exception('Bondlist expects a dataframe with columns "ai" and "aj".')
        inst=cls()
        aiset=set(df.ai)
        ajset=set(df.aj)
        keyset=aiset.union(ajset)
        keys=sorted(list(keyset))
        inst.B={k:[] for k in keys}
        for i,row in df.iterrows():
            inst.B[row.ai].append(row.aj)
            inst.B[row.aj].append(row.ai)
        return inst

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
    
