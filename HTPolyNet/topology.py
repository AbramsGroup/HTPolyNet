import pandas as pd
from HTPolyNet.bondlist import Bondlist

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
# merge_top does not edit molecule_type, system, molecules, defaults
#
class Topology:

    def __init__(self,name=''):
        self.name=name
        ''' D: a dictionay keyed on Gromacs topology directives with values that are lists of
               one or more pandas dataframes corresponding to sections '''
        self.D={}
        ''' bondlist: a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist={}

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
                    padded=tokens[:]
                    # pad with NaN's so that it is same length as series
                    for _ in range(len(tokens),len(header)):
                        padded.append(pd.NA)
                    assert len(padded)==len(header), f'Error: Padding solution does not work! {directive} {len(tokens)}:{len(padded)}!={len(header)} {",".join(tokens)} {",".join(header)}'
                    for k,v in zip(header,padded):
                        series[k].append(v)
                tdf=pd.DataFrame(series)
                if directive=='dihedrals':
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
                retstr+='\t'.join(row[~row.isna()])+'\n'
        return retstr

    def add_bonds(self,pairs=[]):
        at=self.D['atoms']
        ij=self.D['bondtypes'].set_index(['i','j'])
        bmi=self.D['bonds'].set_index(['ai','aj']).index
        newbonds=[]
        for b in pairs:
            ai,aj=min(b),max(b)
            # if this bond is not in the topology
            if not (ai,aj) in bmi:
                newbonds.append([ai,aj])
                it=at[ai-1].type
                jt=at[aj-1].type
                if (it,jt) in ij.index:
                    bt=ij.loc[(it,jt),'func']
                elif (jt,it) in ij.index:
                    bt=ij.loc[(it,jt),'func']
                else:
                    raise Exception(f'no bondtype for {it} {jt} found.')
                # add a new bond!
                self.D['bonds'].append({'ai':ai,'aj':aj,'func':bt,'c0':pd.NA,'c1':pd.NA})
                # update the bondlist
                self.bondlist.append([ai,aj])
            else:
                raise Exception(f'attempt to add already existing bond {ai}-{aj}')
        if len(newbonds)>0:
            self.D['bonds'].sort_values(by=['ai','aj'],inplace=True)
            bmi=self.D['bonds'].set_index(['ai','aj']).index

        ijk=self.D['angletypes'].set_index(['i','j','k'])
        newangles=[]
        for b in newbonds:
            # new angles due to other neighbors of b[0]
            for ai in [i for i in self.bondlist[b[0]] if i!=b[1]]:
                aj=b[0]
                ak=b[1]
                order=[ai,ak]
                ai,ak=min(order),max(order)
                # indices are ordered ai,aj,ak where ai<=ak and aj is the angle vertex
                it=at[ai-1].type
                jt=at[aj-1].type
                kt=at[ak-1].type
                order=[it,kt]
                it,kt=min(order),max(order)
                # types are ordered it,jt,kt where it<kt and jt is the angle vertex
                if (it,jt,kt) in ijk.index:
                    angletype=ijk.loc[(it,jt,kt),'func']
                    self.D['angles'].append({'ai':ai,'aj':aj,'ak':ak,'funct':angletype,'c0':pd.NA,'c1':pd.NA})
                    newangles.append([ai,aj,ak])
                else:
                    raise Exception(f'Angle type {it}-{jt}-{kt} not found.')
            # new angles due to other neighbors of b[1]
            for ak in [k for k in self.bondlist[b[1]] if k!=b[0]]:
                ai=b[0]
                aj=b[1]
                order=[ai,ak]
                ai,ak=min(order),max(order)
                # indices are ordered ai,aj,ak where ai<=ak and aj is the angle vertex
                it=at[ai-1].type
                jt=at[aj-1].type
                kt=at[ak-1].type
                order=[it,kt]
                it,kt=min(order),max(order)
                # types are ordered it,jt,kt where it<kt and jt is the angle vertex
                if (it,jt,kt) in ijk.index:
                    angletype=ijk.loc[(it,jt,kt),'func']
                    self.D['angles'].append({'ai':ai,'aj':aj,'ak':ak,'funct':angletype,'c0':pd.NA,'c1':pd.NA})
                    newangles.append([ai,aj,ak])
                else:
                    raise Exception(f'Angle type {it}-{jt}-{kt} not found.')

        # TODO: Identify and type all new dihedrals        
            for ak in [k for k in self.bondlist[ai] if k!=aj]:
                for al in [l for l in self.bondlist[aj] if l!=ai]:
                    di=ak
                    dj=ai
                    dk=aj
                    dl=al
                    dit=at[di-1].type
                    djt=at[dj-1].type
                    dkt=at[dk-1].type
                    dlt=at[dl-1].type
                    
                    # di,dj,dk,dl are the i,j,k,l defined by Gromacs for a dihedral around the j-k bond.

                    pass # TODO: process ak,ai,aj,al quad
                    # NOTE: the order of the indices for a dihedral entry are i,j,k,l
                    # where the angle is measured between the ijk and jkl planes (that is,
                    # around the j--k bond).  HERE, I'm using python identifiers ai and aj
                    # to represent the diehdral's bond, and ak is the neighbor of ai
                    # while al is the neighbor of aj.  So the order these indices
                    # MUST appear in the topology are ak,ai,aj,al.  To facilitate this, in
                    # the loop body, I use *local* dihedral atom indices di,dj,dk,dl.
            
        nt=self.D['angletypes']
        dt=self.D['dihedraltypes']
        mt=self.D['impropertypes']
        

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

    def merge(self,other):

        # TODO: must update bondlist
        pass
