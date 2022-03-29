import pandas as pd
from HTPolyNet.bondlist import Bondlist

def typeorder(a):
    ''' correctly order the tuple of atom types for particular 
        interaction types to maintain sorted type dataframes '''
    assert type(a)==tuple, 'error: typeorder() requires a tuple argument'
    if len(a)==2: # bond
        return a if a[0]<a[1] else a[::-1]
    elif len(a)==3: # angle
        return a if a[0]<a[2] else a[::-1]
    elif len(a)==4: # dihedral
        return a if a[1]<a[2] else a[::-1]
idxorder=typeorder  # same syntax to order global atom indices in an interaction index


_GromacsIntegers_=('nr','atnum','resnr','ai','aj','ak','al','#nmols','nrexcl','funct','func','nbfunc','comb-rule')
_GromacsFloats_=('charge','mass','chargeB','massB',*tuple([f'c{i}' for i in range(5)]),
                 'b0','kb','th0','cth','rub','kub','phase','kd','pn','fudgeLJ','fudgeQQ')
def typedata(h,s):
    if h in _GromacsIntegers_:
        return int(s)
    if h in _GromacsFloats_:
        return float(s)
    return s

# 'impropers' is not a real gromacs directive, but we use it for our purposes here
_GromacsExtensiveDirectives_=('atoms','pairs','bonds','angles','dihedrals','impropers')

_GromacsTopologyDirectiveHeaders_={
    'atoms':['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass','typeB', 'chargeB', 'massB'],
    'pairs':['ai', 'aj', 'funct'],
    'bonds':['ai', 'aj', 'funct', 'c0', 'c1'],
    'angles':['ai', 'aj', 'ak', 'funct', 'c0', 'c1'],
    'dihedrals':['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'],
    'atomtypes':['name', 'atnum', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'],
    'moleculetype':['name', 'nrexcl'],
    'bondtypes':['i','j','func','b0','kb'],
    'angletypes':['i','j','k','func','th0','cth','rub','kub'],
    'dihedraltypes':['i','j','k','l','func','phase','kd','pn'],
    'system':['Name'],
    'molecules':['Compound','#mols'],
    'defaults':['nbfunc','comb-rule','gen-pairs','fudgeLJ','fudgeQQ']
    }
# dihedral funct==(2,4) means improper

class Topology:

    def __init__(self,name=''):
        self.name=name
        ''' D: a dictionay keyed on Gromacs topology directives with values that are lists of
               one or more pandas dataframes corresponding to sections '''
        self.D={}
        ''' bondlist: a class that owns a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist=Bondlist()

    @classmethod
    def from_topfile(cls,filename):
        inst=cls()
        '''
        Reads a Gromacs-style topology file 'filename' and returns a dictionary keyed on directive names.
        Each value in the dictionary is a pandas dataframe.  Each
        dataframe represents an individual section found with its directive in the file, with columns corresponding to the fields in the section.  Note that the directive 'dihedrals' can appear twice, and we assume a *second* 'dihedrals'
        section enumerates improper dihedrals.
        '''
        inst.D={}
        with open(filename,'r') as f:
            data=f.read().split('[')
            stanzas=[('['+x).split('\n') for x in data][1:]
            for s in stanzas:
                directive=s[0].split()[1].strip()
                contentlines=[l for l in s[1:] if not (len(l.strip())==0 or l.strip().startswith(';'))]
                if not directive in _GromacsTopologyDirectiveHeaders_:
                    raise KeyError(f'unrecognized topology directive "{directive}"')
                header=_GromacsTopologyDirectiveHeaders_[directive]
                series={k:[] for k in header}
                for line in contentlines:
                    if directive!='system':  # no need to split line
                        tokens=[x.strip() for x in line.split()]
                    else:
                        tokens=[line]
                    padded=[typedata(k,_) for _,k in zip(tokens,header)]
                    # pad with NaN's so that it is same length as series
                    for _ in range(len(tokens),len(header)):
                        padded.append(pd.NA)
                    assert len(padded)==len(header), f'Error: Padding solution does not work! {directive} {len(tokens)}:{len(padded)}!={len(header)} {",".join(tokens)} {",".join(header)}'
                    for k,v in zip(header,padded):
                        series[k].append(v)
                tdf=pd.DataFrame(series)
                if directive=='dihedrals':
                    # if there is already a dihedrals section, assume
                    # this new one is for impropers
                    # TODO: Check this against the funct Series
                    if directive in inst.D:
                        inst.D['impropers']=tdf
                    else:
                        inst.D['dihedrals']=tdf
                else:
                    inst.D[directive]=tdf
            inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
            # we must assume the 'atoms' are sorted by global index; however, all other
            # sections need not be sorted.  For convenience, we will keep them sorted by
            # atom indices or atom type name, where appropriate.
            if 'atomtypes' in inst.D:
                inst.D['atomtypes'].sort_values(by='name',inplace=True)
            if 'bonds' in inst.D:
                inst.D['bonds'].sort_values(by=['ai','aj'],inplace=True)
            if 'bondtypes' in inst.D:
                inst.D['bondtypes'].sort_values(by=['i','j'],inplace=True)
            if 'pairs' in inst.D:
                inst.D['pairs'].sort_values(by=['ai','aj'],inplace=True)
            if 'pairtypes' in inst.D:
                inst.D['pairtypes'].sort_values(by=['i','j'],inplace=True)
            if 'angles' in inst.D:
                # central atom (aj) is the primary index
                inst.D['angles'].sort_values(by=['aj','ai','ak'],inplace=True)
            if 'angletypes' in inst.D:
                inst.D['angletypes'].sort_values(by=['j','i','k'],inplace=True)
            if 'dihedrals' in inst.D:
                # central atoms (aj,ak) are the primary index
                inst.D['dihedrals'].sort_values(by=['aj','ak','ai','al'],inplace=True)
            if 'dihedraltypes' in inst.D:
                inst.D['dihedraltypes'].sort_values(by=['j','k','i','l'],inplace=True)
            if 'impropers' in inst.D:
                inst.D['impropers'].sort_values(by=['aj','ak','ai','al'],inplace=True)
            return inst

    @classmethod
    def from_ex(cls,other):
        ''' make a new Topology instance by copying only the extensive dataframes 
            from an existing topology '''
        inst=cls()
        for t in _GromacsExtensiveDirectives_:
            inst.D[t]=other.D[t].copy()
        return inst

    def __str__(self):
        ''' Generates a string in the proper top format '''
        retstr=''
        for k,vv in self.D.items():
            if vv.empty:
                ''' an empty stanza will not be output '''
                continue
            ''' our internal directive 'impropers' must be
                reported as 'dihedrals' in the top file '''
            retstr+='[ '+(k if k!='impropers' else 'dihedrals')+' ]\n'
            retstr+='; '+'\t'.join(vv.columns)+'\n'
            for i,row in vv.iterrows():
                ''' assumes NaN's are only allowed in trailing columns '''
                retstr+='\t'.join(row[~row.isna()].apply(str))+'\n'
        return retstr

    def write(self,outfile):
        with open(outfile,'w') as f:
            f.write(str(self))

    def total_charge(self):
        if 'atoms' in self.D:
            return self.D['atoms']['charge'].sum()
        return 0.0

    def atomcount(self):
        if 'atoms' in self.D:
            return len(self.D['atoms'])
        return 0

    def replicate(self,n):
        spawned=Topology()
        for i in range(n):
            s=Topology.from_ex(self)
            spawned.merge(s)
        self.merge(spawned)

    def add_bonds(self,pairs=[]):
        ''' add bonds to a topology
            pairs:  list of 2-tuples of atom global indices '''
        at=self.D['atoms']
        ij=self.D['bondtypes'].set_index(['i','j'])
        bmi=self.D['bonds'].set_index(['ai','aj']).index
        pmi=self.D['pairs'].set_index(['ai','aj']).index
        newbonds=[]
        for b in pairs:
            ai,aj=idxorder(b)
            # if this bond is not in the topology
            if not (ai,aj) in bmi:
                newbonds.append([ai,aj])
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                idx=typeorder((it,jt))
                if idx in ij.index:
                    bt=ij.loc[idx,'func']
                else:
                    raise Exception(f'no bondtype {idx} found.')
                # add a new bond!
                h=_GromacsTopologyDirectiveHeaders_['bonds']
                data=[ai,aj,bt,pd.NA,pd.NA]
                assert len(h)==len(data), 'Error: not enough data for new bond?'
                bonddict={k:[v] for k,v in zip(h,data)}
                pd.concat((self.D['bonds'],pd.DataFrame(bonddict)),ignore_index=True)
                # update the bondlist
                self.bondlist.append([ai,aj])
                # remove this pair from pairs if it's in there (it won't be)
                if idx in pmi:
                    d=self.D['pairs']
                    indexes_to_drop=d[(d.ai.isin(idx))&(d.aj.isin(idx))].index
                    indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                    self.D['pairs']=d.take(list(indexes_to_keep)).reset_index(drop=True)
            else:
                raise Exception(f'attempt to add already existing bond {ai}-{aj}')
        ''' maintain sort of bonds df if new bonds were added '''
        if len(newbonds)>0:
            self.D['bonds'].sort_values(by=['ai','aj'],inplace=True)

        ijk=self.D['angletypes'].set_index(['i','j','k'])
        ijkl=self.D['dihedraltypes'].set_index(['i','j','k','l'])
        newangles=[]
        newdihedrals=[]
        for b in newbonds:
            ''' new angles due to other neighbors of b[0] '''
            for ai in [i for i in self.bondlist.partners_of(b[0]) if i!=b[1]]:
                aj=b[0]
                ak=b[1]
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                kt=at.iloc[ak-1].type
                idx=typeorder((it,jt,kt))
                if idx in ijk.index:
                    angletype=ijk.loc[idx,'func']
                    i,j,k=idxorder((ai,aj,ak))
                    h=_GromacsTopologyDirectiveHeaders_['angles']
                    data=[i,j,k,angletype,pd.NA,pd.NA]
                    assert len(h)==len(data), 'Error: not enough data for new angle?'
                    angledict={k:[v] for k,v in zip(h,data)}
                    pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                    newangles.append([i,j,k])
                else:
                    # no longer exception but warning
                    raise Exception(f'Angle type {idx} not found.')
            ''' new angles due to other neighbors of b[1] '''
            for ak in [k for k in self.bondlist.partners_of(b[1]) if k!=b[0]]:
                ai=b[0]
                aj=b[1]
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                kt=at.iloc[ak-1].type
                idx=typeorder((it,jt,kt))
                if idx in ijk.index:
                    angletype=ijk.loc[idx,'func']
                    i,j,k=idxorder((ai,aj,ak))
                    h=_GromacsTopologyDirectiveHeaders_['angles']
                    data=[i,j,k,angletype,pd.NA,pd.NA]
                    assert len(h)==len(data), 'Error: not enough data for new angle?'
                    angledict={k:[v] for k,v in zip(h,data)}
                    pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                    newangles.append([i,j,k])
                else:
                    raise Exception(f'Angle type {idx} not found.')

            ''' DIHEDRALS i-j-k-l 
                j-k is the torsion bond
                angle is measured between the i-j-k plane and the j-k-l plane 
                a new bond could be one of i-j, j-k, or k-l
            '''

            ''' new proper dihedrals for which the new bond is the central j-k bond '''
            aj,ak=idxorder(b)        
            for ai in [i for i in self.bondlist.partners_of(aj) if i!=ak]:
                for al in [l for l in self.bondlist.partners_of(ak) if l!=aj]:
                    it=at.iloc[ai-1].type
                    jt=at.iloc[aj-1].type
                    kt=at.iloc[ak-1].type
                    lt=at.iloc[al-1].type
                    idx=typeorder((it,jt,kt,lt))
                    if idx in ijkl.index:
                        dihedtype=ijkl.loc[idx,'func']
                        i,j,k,l=idxorder((ai,aj,ak,al))
                        h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                        data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                        assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                        diheddict={k:[v] for k,v in zip(h,data)}
                        pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                        newdihedrals.append([i,j,k,l])
                else:
                    raise Exception(f'Dihedral type {idx} not found.')   
        
            ''' new proper dihedrals for which the new bond is the i-j or j-i bond '''
            for ai,aj in zip(b,reversed(b)):
                for ak in [k for k in self.bondlist.partners_of(aj) if k!=ai]:
                    for al in [l for l in self.bondlist.partners_of(ak) if l!=ak]:
                        it=at.iloc[ai-1].type
                        jt=at.iloc[aj-1].type
                        kt=at.iloc[ak-1].type
                        lt=at.iloc[al-1].type
                        idx=typeorder((it,jt,kt,lt))
                        if idx in ijkl.index:
                            dihedtype=ijkl.loc[idx,'func']
                            i,j,k,l=idxorder((ai,aj,ak,al))
                            h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                            data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                            assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                            diheddict={k:[v] for k,v in zip(h,data)}
                            pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                            newdihedrals.append([i,j,k,l])
                        else:
                            raise Exception(f'Dihedral type {idx} not found.')   
        
            ''' new proper dihedrals for which the new bond is the k-l or l-k bond '''
            for ak,al in zip(b,reversed(b)):
                for aj in [j for j in self.bondlist.partners_of(ak) if j!=al]:
                    for ai in [i for i in self.bondlist.partners_of(aj) if i!=ak]:
                        it=at.iloc[ai-1].type
                        jt=at.iloc[aj-1].type
                        kt=at.iloc[ak-1].type
                        lt=at.iloc[al-1].type
                        idx=typeorder((it,jt,kt,lt))
                        if idx in ijkl.index:
                            dihedtype=ijkl.loc[idx,'func']
                            i,j,k,l=idxorder((ai,aj,ak,al))
                            h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                            data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                            assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                            diheddict={k:[v] for k,v in zip(h,data)}
                            pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                            newdihedrals.append([i,j,k,l])
                        else:
                            raise Exception(f'Dihedral type {idx} not found.')   

            ''' new improper dihedrals '''
            # TODO
            # I don't know if this is necessary

    def delete_atoms(self,idx=[],reindex=True):
        d=self.D['atoms']
        indexes_to_drop=d[d.nr.isin(idx)].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['atoms']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['atoms']
            oldGI=d['nr'].copy()
            d['old_nr']=oldGI
            d['nr']=d.index+1
            mapper={k:v for k,v in zip(d['old_nr'],d['nr'])}
            d['nr_shift']=d['nr']-oldGI  # probably not necessary
        d=self.D['bonds']
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['bonds']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['bonds']
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
            self.bondlist=Bondlist.fromDataFrame(d)
        d=self.D['pairs']
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['pairs']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['pairs']
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
        d=self.D['angles']
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))|(d.ak.isin(idx))].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['angles']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['angles']
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
            d.ak=d.ak.map(mapper)
        for four_body_type in ['dihedrals','impropers']:
            d=self.D[four_body_type]
            indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))|(d.ak.isin(idx))|(d.al.isin(idx))].index
            indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
            self.D[four_body_type]=d.take(list(indexes_to_keep)).reset_index(drop=True)
            if reindex:
                d=self.D[four_body_type]
                d.ai=d.ai.map(mapper)
                d.aj=d.aj.map(mapper)
                d.ak=d.ak.map(mapper)
                d.al=d.al.map(mapper)

    def _myconcat(self,other,directive='',idxlabel=[],idxshift=0,drop_duplicates=False):
        if not directive in other.D:
            return
        if directive in self.D:
            # shift atom indices
            for i in idxlabel:
                other.D[directive][i]+=idxshift
            if drop_duplicates:
                self.D[directive]=pd.concat((self.D[directive],other.D[directive]),ignore_index=True).drop_duplicates()
            else:
                self.D[directive]=pd.concat((self.D[directive],other.D[directive]),ignore_index=True)
        else:
            self.D[directive]=other.D[directive]

    def merge(self,other):
        self.merge_ex(other)
        ''' merge types but drop duplicates '''
        for t in ['atomtypes','bondtypes','angletypes','dihedraltypes']:
            if t in self.D:
                self._myconcat(other,directive=t,drop_duplicates=True)
            else:
                self.D[t]=other.D[t]

    def merge_ex(self,other):
        ''' merge EXTENSIVE quantities '''
        idxshift=0 if 'atoms' not in self.D else len(self.D['atoms'])
        self._myconcat(other,directive='atoms',idxlabel=['nr'],idxshift=idxshift)
        self._myconcat(other,directive='bonds',idxlabel=['ai','aj'],idxshift=idxshift)
        self.bondlist.update(other.D['bonds'])
        self._myconcat(other,directive='pairs',idxlabel=['ai','aj'],idxshift=idxshift)
        self._myconcat(other,directive='angles',idxlabel=['ai','aj','ak'],idxshift=idxshift)
        self._myconcat(other,directive='dihedrals',idxlabel=['ai','aj','ak','al'],idxshift=idxshift)
        self._myconcat(other,directive='impropers',idxlabel=['ai','aj','ak','al'],idxshift=idxshift)

    def get_atom(self,idx):
        return self.D['atoms'].iloc[idx-1]