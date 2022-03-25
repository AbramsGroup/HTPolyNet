''' a simple class for constructing a bondlist dictionary '''

class Bondlist:
    def __init__(self):
        self.B={}
    
    @classmethod
    def fromDataFrame(cls,df):
        if not 'ai' in df.columns and not 'aj' in df.columns:
            raise Exception('Cannot build bondlist from non-bonds DF.')
        inst=cls()
        keys=list(set(df.ai)).sort()
        inst.B={k:[] for k in keys}
        for i,row in df.iterrows():
            inst.B[row.ai].append(row.aj)
        return inst

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
    
