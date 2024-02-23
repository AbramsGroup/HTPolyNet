"""

.. module:: chain
   :synopsis: Manages data structure that keeps track of polymerized chains
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
logger=logging.getLogger(__name__)

class Chain:
    def __init__(self,head=-1,tail=-1,idx=-1):
        self.idx_list=[]
        if head!=-1:
            self.idx_list.append(head)
        if tail!=-1:
            self.idx_list.append(tail)
        self.is_cyclic=False
        self.idx=idx
    def is_head(self,a):
        return not self.is_cyclic and a==self.idx_list[0]
    def is_tail(self,a):
        return not self.is_cyclic and a==self.idx_list[-1]
    def merge(self,other):
        # tail of self bonds to head of other
        self.idx_list.extend(other.idx_list)

    def remap_idx(self,idx_mapper):
        self.idx_list=[idx_mapper[x] for x in self.idx_list]
    def shift(self,shift):
        self.idx_list=[x+shift for x in self.idx_list]

class ChainManager:
    def __init__(self,**kwargs):
        self.chains=[]
        self.create_if_missing=kwargs.get('create_if_missing',False)

    def chain_of(self,a):
        for c in self.chains:
            if a in c.idx_list:
                return c
        return None
    
    def new_chain(self,head,tail):
        self.chains.append(Chain(head=head,tail=tail,idx=len(self.chains)))

    def mergechains(self,dst,src):
        dst.merge(src)
        self.chains.remove(src)
        self.reindex()

    def injest_other(self,other):
        self.chains.extend(other.chains)
        self.reindex()

    def reindex(self):
        for i,c in enumerate(self.chains):
            c.idx=i

    def shift(self,shift):
        for c in self.chains:
            c.shift(shift)

    def remap(self,idx_mapper):
        for c in self.chains:
            c.remap_idx(idx_mapper)

    def injest_bond(self,ai,aj):
        ic=self.chain_of(ai)
        jc=self.chain_of(aj)
        if not ic and not jc:
            if self.create_if_missing:
                logger.debug(f'New chain with head {ai} and tail {aj}')
                self.new_chain(ai,aj)
            else:
                pass
        elif not ic:
            logger.debug(f'j-Atom {aj} (chain {jc.idx}) bonds to atom {ai}, but {ai} is not yet in a chain.')
            logger.debug(f'This should never happen in a C=C polymerizing system; all reactive Cs are in chains from the very beginning!')
            raise Exception('This is a bug - no i-chain!')
        elif not jc:
            logger.debug(f'i-Atom {ai} (chain {ic.idx}) bonds to atom {aj}, but {aj} is not yet in a chain.')
            logger.debug(f'This should never happen in a C=C polymerizing system; all reactive Cs are in chains from the very beginning!')
            raise Exception('This is a bug - no j-chain!')
        elif ic==jc:
            assert (ic.is_head(ai) and ic.is_tail(aj)) or (ic.is_head(aj) and ic.is_tail(ai))
            ic.is_cyclic=True
            logger.debug(f'New bond between {ai} and {aj} creates a cyclic chain of length ({len(ic.idx_list)})')
        elif ic.is_head(ai) and jc.is_tail(aj):
            assert not ic.is_cyclic and not jc.is_cyclic
            self.mergechains(jc,ic)
        elif jc.is_head(aj) and ic.is_tail(ai):
            assert not ic.is_cyclic and not jc.is_cyclic
            self.mergechains(ic,jc)

    def injest_bonds(self,bondlist):
        for b in bondlist:
            self.injest_bond(b[0],b[1])

    def to_dataframe(self,D,lidx_col='bondchain_idx',gidx_col='bondchain'):
        assert lidx_col in D and gidx_col in D,f'Dataframe missing column {lidx_col} or {gidx_col}'
        # wipe out df
        D[lidx_col]=[-1]*D.shape[0]
        D[gidx_col]=[-1]*D.shape[0]
        for gidx,c in enumerate(self.chains):
            for lidx,atidx in enumerate(c.idx_list):
                df_idx=atidx-1 # should be true always
                D.loc[df_idx,lidx_col]=lidx
                D.loc[df_idx,gidx_col]=gidx

    def from_dataframe(self,D,lidx_col='bondchain_idx',gidx_col='bondchain'):
        # overwrites whole thing
        self.chains=[]
        heads=D[D[lidx_col]==0].index+1 # always!
        headchains=D[D[lidx_col]==0][gidx_col]
        for h,hc in zip(heads,headchains):
            self.chains.append(Chain(head=h,idx=hc))
        self.chains.sort(key=lambda x: x.idx)
        old_idx=[x.idx for x in self.chains]
        self.reindex()
        assert all([(x==y) for x,y in zip([x.idx for x in self.chains],list(range(len(self.chains))))])
        next_idx=1
        while next_idx<D.shape[0]:
            if D[D[lidx_col]==next_idx].shape[0]==0:
                break
            next_ats=D[D[lidx_col]==next_idx].index+1 # always
            next_atchains=D[D[lidx_col]==next_idx][gidx_col]
            for a,ac in zip(next_ats,next_atchains):
                acidx=old_idx.index(ac)
                self.chains[acidx].idx_list.append(a)
            next_idx+=1


    