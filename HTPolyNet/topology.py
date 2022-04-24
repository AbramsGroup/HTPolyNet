import pandas as pd
import logging
from HTPolyNet.bondlist import Bondlist
import os
from copy import deepcopy
from scipy.constants import physical_constants

def typeorder(a):
    ''' correctly order the tuple of atom types for particular
        interaction types to maintain sorted type dataframes '''
    assert type(a)==tuple, 'error: typeorder() requires a tuple argument'
    if len(a)==2: # bond
        return a if a[0]<a[1] else a[::-1]
    elif len(a)==3: # angle
        return a if a[0]<a[2] else a[::-1]
    elif len(a)==4: # dihedral
        if a[1]==a[2]:
            return a if a[0]<a[3] else a[::-1]
        return a if a[1]<a[2] else a[::-1]
idxorder=typeorder  # same syntax to order global atom indices in an interaction index

def df_typeorder(df,typs):
    for i in df.index:
        df.loc[i,typs]=typeorder(tuple(df.loc[i,typs]))


_GromacsIntegers_=('nr','atnum','resnr','ai','aj','ak','al','#mols','nrexcl','funct','func','nbfunc','comb-rule')
_GromacsFloats_=('charge','mass','chargeB','massB',*tuple([f'c{i}' for i in range(5)]),
                 'b0','kb','th0','cth','rub','kub','phase','kd','pn','fudgeLJ','fudgeQQ')
def typedata(h,s):
    if h in _GromacsIntegers_:
        return int(s)
    if h in _GromacsFloats_:
        return float(s)
    return s

_GromacsExtensiveDirectives_=('atoms','pairs','bonds','angles','dihedrals')
_NonGromacsExtensiveDirectives_=('mol2_bonds')
_GromacsTopologyDirectiveOrder_=['defaults','atomtypes','bondtypes','angletypes','dihedraltypes','moleculetype','atoms','pairs','bonds','angles','dihedrals','system','molecules']
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
_GromacsTopologyHashables_={
    'atoms':['nr'],
    'pairs':['ai', 'aj'],
    'bonds':['ai', 'aj'],
    'angles':['ai', 'aj', 'ak'],
    'dihedrals':['ai', 'aj', 'ak', 'al',],
    'atomtypes':['name'],
    'bondtypes':['i','j'],
    'angletypes':['i','j','k'],
    'dihedraltypes':['i','j','k','l']
    }
_GromacsTopologyDataFields_={
    'atoms':['type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass','typeB', 'chargeB', 'massB'],
    'pairs':['funct'],
    'bonds':['funct', 'c0', 'c1'],
    'angles':['funct', 'c0', 'c1'],
    'dihedrals':['funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'],
    'atomtypes':['atnum', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'],
    'bondtypes':['func','b0','kb'],
    'angletypes':['func','th0','cth','rub','kub'],
    'dihedraltypes':['func','phase','kd','pn']
    }
_GromacsTopologyDirectiveDefaults_={
    'system':['A_generic_system'],
    'molecules':['None',1],
    'moleculetype':['None',3],
#    'defaults':[1,2,'no',0.5,0.83333333]
    'defaults':[1,2,'yes',0.5,0.83333333]
}
# dihedral funct==(2,4) means improper

class Topology:
    ''' Gromacs topology handler '''
    def __init__(self,system=''):
        ''' D: a dictionay keyed on Gromacs topology directives with values that are lists of
               one or more pandas dataframes corresponding to sections '''
        self.D={}
        for k,v in _GromacsTopologyDirectiveDefaults_.items():
            hdr=_GromacsTopologyDirectiveHeaders_[k]
            dfdict={k:[a] for k,a in zip(hdr,v)}
            self.D[k]=pd.DataFrame(dfdict)
        self.D['system']=pd.DataFrame({'name':[system]})
        ''' bondlist: a class that owns a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist=Bondlist()
        self.empty=True

    @classmethod
    def read_gro(cls,filename):
        assert os.path.exists(filename), f'Error: {filename} not found.'
        inst=cls()
        inst.filename=filename
        '''
        Reads a Gromacs-style topology file 'filename' and returns a dictionary keyed on directive names.
        Each value in the dictionary is a pandas dataframe.  Each
        dataframe represents an individual section found with its directive in the file, with columns corresponding to the fields in the section.  Note that the we allow for input topology/itp files to have two 'dihedrals' and 'dihedraltypes' sections; these
        are merged in the result.
        '''
        inst.D={}
        inst.includes=[]
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
                    if line.startswith('#include'):
                        tokens=[x.strip() for x in line.split()]
                        inst.includes.append(tokens[1].replace('"',''))
                        continue
                    # ignore anything after a ';'
                    line=line.split(';')[0].strip()
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
                if directive=='dihedraltypes':
                    if directive in inst.D:
                        logging.info(f'Found second set of {len(tdf)} [ dihedraltypes ] in {inst.filename}; merging into set of {len(inst.D["dihedraltypes"])} types already read in...')
                        # we have already read-in a dihedraltypes section
                        # so let's append this one
                        inst.D['dihedraltypes']=pd.concat([inst.D['dihedraltypes'],tdf],ignore_index=True).drop_duplicates()
                        logging.info(f'    -> now there are {len(inst.D["dihedraltypes"])} dihedral types.')
                    else:
                        inst.D[directive]=tdf
                elif directive=='dihedrals':
                    # if there is already a dihedrals section, assume
                    if directive in inst.D:
                        inst.D['dihedrals']=pd.concat([inst.D['dihedrals'],tdf],ignore_index=True)
                    else:
                        inst.D[directive]=tdf
                else:
                    inst.D[directive]=tdf
            # we must assume the 'atoms' are sorted by global index; however, all other
            # sections need not be sorted.  For convenience, we will keep them sorted by
            # atom indices or atom type name, where appropriate.
            # TODO: sort each row's hashables before sorting whole dataframes
            if 'atomtypes' in inst.D:
                inst.D['atomtypes'].sort_values(by='name',inplace=True)
            if 'bonds' in inst.D:
                inst.D['bonds'].sort_values(by=['ai','aj'],inplace=True)
                inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
            if 'bondtypes' in inst.D:
                df_typeorder(inst.D['bondtypes'],typs=['i','j'])
                inst.D['bondtypes'].sort_values(by=['i','j'],inplace=True)
            if 'pairs' in inst.D:
                inst.D['pairs'].sort_values(by=['ai','aj'],inplace=True)
            if 'pairtypes' in inst.D:
                df_typeorder(inst.D['pairtypes'],typs=['i','j'])
                inst.D['pairtypes'].sort_values(by=['i','j'],inplace=True)
            if 'angles' in inst.D:
                # central atom (aj) is the primary index
                inst.D['angles'].sort_values(by=['aj','ai','ak'],inplace=True)
            if 'angletypes' in inst.D:
                df_typeorder(inst.D['angletypes'],typs=['i','j','k'])
                inst.D['angletypes'].sort_values(by=['j','i','k'],inplace=True)
            if 'dihedrals' in inst.D:
                # central atoms (aj,ak) are the primary index
                inst.D['dihedrals'].sort_values(by=['aj','ak','ai','al'],inplace=True)
            if 'dihedraltypes' in inst.D:
                df_typeorder(inst.D['dihedraltypes'],typs=['i','j','k','l'])
                # print(f'    -> pre sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
                inst.D['dihedraltypes'].sort_values(by=['j','k','i','l'],inplace=True)
                # print(f'    -> post sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
            for f in inst.includes:
                # print(f'reading included topology {f}')
                inst.merge(Topology.read_gro(f))
            # if replicate>0:
            #     inst.rep_ex(replicate)
            # print(f'{filename}',inst.D.keys())
            #assert 'defaults' in inst.D, f'Error: no [ defaults ] in {filename}'
#            logging.debug(f'Checking for duplicates in just-read-in gro topology {filename}')
            inst.dup_check(die=False)
 #           logging.debug(f'read_gro ends')
            inst.empty=False
            return inst

    def bond_source_check(self):
        if 'bonds' in self.D and 'mol2_bonds' in self.D:
            logging.info(f'Bond data source check requested.')
            bmi=self.D['bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
            mbmi=self.D['mol2_bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
            check=all([x==y for x,y in zip(bmi,mbmi)])
            logging.info(f'Result: {check}')
            logging.info(f'GROMACS:')
            logging.info(self.D['bonds'].to_string())
            logging.info(f'MOL2:')
            logging.info(self.D['mol2_bonds'].to_string())

    def has_bond(self,pair):
        bmi=self.D['bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
        mbmi=self.D['mol2_bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
        print(bmi,mbmi)
        return pair in bmi and pair in mbmi

    def shiftatomsidx(self,idxshift,directive,rows=[],idxlabels=[]):
        ''' shift all global atom indices (referenced by labels in idxlables[]) '''
        if directive in self.D:
            cols=self.D[directive].columns.get_indexer(idxlabels)
            # print(f'directive {directive} idxlabels {idxlabels} idxshift {idxshift} rows {rows} cols {cols}')
            self.D[directive].iloc[rows[0]:rows[1],cols]+=idxshift

    def rep_ex(self,count=0):
        ''' replicate extensive components (atoms, pairs, bonds, angles, dihedrals) '''
        if count>0:
            counts={k:0 for k in _GromacsExtensiveDirectives_}
            for t in _GromacsExtensiveDirectives_:
                if t in self.D:
                    counts[t]=len(self.D[t])
                for t in _NonGromacsExtensiveDirectives_:
                    counts[t]=len(self.D[t])
            try:
                idxshift=counts['atoms']
            except:
                raise Exception(f'Error: expected an "atoms" dataframe')
            for t in _GromacsExtensiveDirectives_:
                if t in self.D:
                    # print(f'replicating {t} by {count}')
                    self.D[t]=pd.concat([self.D[t]]*count,ignore_index=True)
            # print(f'new raw atom count (pre-index-shifted) {len(self.D["atoms"])}')
            for c in range(1,count):
                # if c%100 == 0:
                # print(f'  -> shifting indices in replica {c}...')
                # print(f'     -> atoms')
                self.shiftatomsidx(idxshift*c,'atoms',rows=[(c*counts['atoms']),((c+1)*counts['atoms'])],idxlabels=['nr'])
                self.shiftatomsidx(c,'atoms',rows=[(c*counts['atoms']),((c+1)*counts['atoms'])],idxlabels=['resnr'])
                # print(f'     -> bonds')
                self.shiftatomsidx(idxshift*c,'bonds',rows=[(c*counts['bonds']),((c+1)*counts['bonds'])],idxlabels=['ai','aj'])
                self.shiftatomsidx(idxshift*c,'mol2_bonds',rows=[(c*counts['mol2_bonds']),((c+1)*counts['mol2_bonds'])],idxlabels=['ai','aj'])
                # print(f'     -> pairs')
                self.shiftatomsidx(idxshift*c,'pairs',rows=[(c*counts['pairs']),((c+1)*counts['pairs'])],idxlabels=['ai','aj'])
                self.shiftatomsidx(idxshift*c,'angles',rows=[(c*counts['angles']),((c+1)*counts['angles'])],idxlabels=['ai','aj','ak'])
                self.shiftatomsidx(idxshift*c,'dihedrals',rows=[(c*counts['dihedrals']),((c+1)*counts['dihedrals'])],idxlabels=['ai','aj','ak','al'])
            if 'bonds' in self.D:
                self.bondlist=Bondlist.fromDataFrame(self.D['bonds'])

    @classmethod
    def from_ex(cls,other):
        ''' make a new Topology instance by copying only the extensive dataframes
            from an existing topology '''
        inst=cls()
        for t in _GromacsExtensiveDirectives_:
            if t in other.D:
                inst.D[t]=other.D[t].copy()
        return inst

    def to_file(self,filename=''):
        if filename=='':
            return
        with open(filename,'w') as f:
            f.write('; Gromacs-format topology written by HTPolyNet\n')
        assert 'defaults' in self.D, 'Error: no [ defaults ] in topology?'
        for k in _GromacsTopologyDirectiveOrder_:
            if k in self.D:
                with open(filename,'a') as f:
                    f.write(f'[ {k} ]\n; ')
                if k in _GromacsTopologyHashables_:
                    odf=self.D[k].sort_values(by=_GromacsTopologyHashables_[k])
                    odf.to_csv(filename,sep=' ',mode='a',index=False,header=True,doublequote=False)
                else:
                    self.D[k].to_csv(filename,sep=' ',mode='a',index=False,header=True,doublequote=False)
                with open(filename,'a') as f:
                    f.write('\n')
        with open(filename,'a') as f:
            f.write('; end\n')

    def total_charge(self):
        if 'atoms' in self.D:
            return self.D['atoms']['charge'].sum()
        return 0.0

    def adjust_charges(self,desired_charge=0.0):
        apparent_charge=self.total_charge()
        overcharge=apparent_charge-desired_charge
        logging.info(f'Adjusting charges due to overcharge of {overcharge}')
        cpa=-overcharge/len(self.D['atoms'])
        self.D['atoms']['charge']+=cpa
        logging.info(f'New total charge after adjustment: {self.total_charge()}')
        return self
        
    def total_mass(self,units='gromacs'):
        fac=1.0
        if units=='SI':
            fac=physical_constants['atomic mass constant'][0]
        if 'atoms' in self.D:
            M_amu=self.D["atoms"]["mass"].sum()
            logging.info(f'total mass says {M_amu}')
            return M_amu*fac

    def atomcount(self):
        if 'atoms' in self.D:
            return len(self.D['atoms'])
        return 0

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
                newbonds.append((ai,aj))
                logging.debug(f'asking types of {ai} and {aj}; at.shape {at.shape}')
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
                self.D['bonds']=pd.concat((self.D['bonds'],pd.DataFrame(bonddict)),ignore_index=True)
                logging.info(f'just added {bonddict}')
                # update the bondlist
                self.bondlist.append([ai,aj])
                if 'mol2_bonds' in self.D:
                    data=[len(self.D['mol2_bonds']),ai,aj,1] # assume single bond
                    bonddict={k:[v] for k,v in zip(['bondIdx','ai','aj','type'],data)}
                    self.D['mol2_bonds']=pd.concat((self.D['mol2_bonds'],pd.DataFrame(bonddict)),ignore_index=True)
                # remove this pair from pairs if it's in there (it won't be)
                if idx in pmi:
                    d=self.D['pairs']
                    indexes_to_drop=d[(d.ai.isin(idx))&(d.aj.isin(idx))].index
                    indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                    self.D['pairs']=d.take(list(indexes_to_keep)).reset_index(drop=True)
            else:
                raise Exception(f'attempt to add already existing bond {ai}-{aj}')
        ''' maintain sort of bonds df if new bonds were added '''
        # if len(newbonds)>0:
        #     self.D['bonds'].sort_values(by=['ai','aj'],inplace=True)

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
                    self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                    newangles.append([i,j,k])
                else:
                    # no longer exception but warning
                    logging.warning(f'Angle type {idx} not found.  Hopefully you are about to parameterize!')
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
                    self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                    newangles.append([i,j,k])
                else:
                    logging.warning(f'Angle type {idx} not found.  Hopefully you are about to parameterize!')

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
                        self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                        newdihedrals.append([i,j,k,l])
                else:
                    logging.warning(f'Dihedral type {idx} not found. Hopefully you are about to parameterize!')

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
                            self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                            newdihedrals.append([i,j,k,l])
                        else:
                            logging.warning(f'Dihedral type {idx} not found. Hopefully you are about to parameterize!')

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
                            self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                            newdihedrals.append([i,j,k,l])
                        else:
                            logging.warning(f'Dihedral type {idx} not found. Hopefully you are about to parameterize!')


    def delete_atoms(self,idx=[],reindex=True):
        #logging.debug(f'Delete atoms: {idx}')
        d=self.D['atoms']
        indexes_to_drop=d[d.nr.isin(idx)].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['atoms']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['atoms']
            oldGI=d['nr'].copy()
            d['nr']=d.index+1
            mapper={k:v for k,v in zip(oldGI,d['nr'])}
            #logging.debug(f'delete_atoms: mapper {mapper}')
            #d['nr_shift']=d['nr']-oldGI  # probably not necessary
        ptt=['bonds','mol2_bonds','pairs']
        for pt in ptt:
            if pt in self.D:
                d=self.D[pt]
                #logging.debug(f'delete atom: {pt} df prior to deleting')
                #logging.debug(d.to_string())
                indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
                #logging.debug(f'dropping {pt} {indexes_to_drop}')
                indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                self.D[pt]=d.take(list(indexes_to_keep)).reset_index(drop=True)
                if reindex:
                    d=self.D[pt]
                    #logging.debug(f'delete atom: {pt} df prior to reindexing')
                    #logging.debug(d.to_string())
                    d.ai=d.ai.map(mapper)
                    d.aj=d.aj.map(mapper)
                    if pt=='bonds':
                        #logging.debug(f'delete atom: bondlist remake from')
                        #logging.debug(d.to_string())
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
        for four_body_type in ['dihedrals']:
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
        # TODO: duplicate dihedral type check
        logging.debug('Topology.merge begins')
        self.merge_ex(other)
        self.merge_types(other)
        logging.debug('Topology.merge ends')

    def force_overrides(self,other,types):
        # If there are any dihedraltypes in other with same atomtypes as any dihedraltypes
        # in self but with DIFFERENT parameters, IMMEDIATELY copy other's parameters to
        # the override fields of the actual dihedrals in other and then DROP the
        # dihedral types from other.
        logging.debug(f'force_overrides begins')
        for typ,ext in types:
            if typ not in self.D or typ not in other.D:
                continue
            typfields=_GromacsTopologyHashables_[typ]
            datafields=_GromacsTopologyDataFields_[typ]
            idxfields=_GromacsTopologyHashables_[ext]
            stdf=self.D[typ]
            otdf=other.D[typ]
            oedf=other.D[ext]
            mdf=pd.concat((stdf,otdf),ignore_index=True).drop_duplicates()
            dups=mdf.duplicated(subset=typfields,keep='first')
            if any(dups):
                logging.debug(f'Duplicate {typ}(s) with differing parameters:')
                logging.debug('\n'+mdf[dups].to_string())
                for i,p in mdf[dups].iterrows():
                    # get the atom types
                    if typ=='dihedraltypes' and p.func==9:
                        continue
                    typvalues=p[typfields]
                    logging.debug(f'Searching {ext} for type(s) {list(typvalues)}')
                    for j,q in oedf.iterrows():
#                        logging.debug(f'\ni {i}\nq {q}\nq[idxfields] {q[idxfields]}')
                        # find any extensive entry in other whose atom indices refer to atoms of the stipulated types
                        typs=all([self.get_atomtype(x)==y for x,y in zip(list(q[idxfields]),typvalues)])
                        if typs:
                            logging.debug(f'Found {list(q[idxfields])}')
                            # overwrite datafields in the extensive dataframe entry for other
                            logging.debug(f'Existing extensive entry:\n{oedf.iloc[j].to_string()}')
                            #oedf.iloc[j][datafields]=p[datafields]
                            # logging.debug(f'After overwrite:\n{oedf.iloc[j]}')
                            # # drop the type entry in the other dataframe
                            # otdf.set_index(typfields)
                            # idx=otdf[typ].index
                            # otdf.drop(idx,inplace=True)
                            # otdf.reset_index().reindex()
            logging.debug('\n'+otdf.to_string())
        logging.debug('force_overrides ends')

    def dup_check(self,die=True):
        logging.debug('dup_check begins')
        L=['atomtypes','bondtypes','angletypes','dihedraltypes']
        for t in L:
            i=_GromacsTopologyHashables_[t]
            ''' checking for types with duplicate atom-type indices '''
            dups=self.D[t].duplicated(subset=i,keep=False)
            if any(dups):
                logging.error(f'Duplicate {t} with different parameters detected\n'+self.D[t][dups].to_string())
                # TODO: for duplicates, find the actual dihedrals which contributed them and override parameters
                # how the hell am i going to go that??
                if die:
                    raise Exception('duplicate topology types with different parameters detected')
        logging.debug('dup_check ends')

    def merge_types(self,other):
        ''' merge types but drop duplicates '''
        L=['atomtypes','bondtypes','angletypes','dihedraltypes']
#        self.force_overrides(other,[('dihedraltypes','dihedrals')])
        for t in L:
            self._myconcat(other,directive=t,drop_duplicates=True)


    def merge_ex(self,other):
        # print('   extensive merging...')
        ''' merge EXTENSIVE quantities '''
        idxshift=0 if 'atoms' not in self.D else len(self.D['atoms'])
        rdxshift=0 if 'atoms' not in self.D else self.D['atoms'].iloc[-1]['resnr']
        if 'atoms' in other.D:
            other.D['atoms']['resnr']+=rdxshift
        self._myconcat(other,directive='atoms',idxlabel=['nr'],idxshift=idxshift)
        self._myconcat(other,directive='bonds',idxlabel=['ai','aj'],idxshift=idxshift)
        if 'bonds' in other.D:
            self.bondlist.update(other.D['bonds'])
        if 'mol2_bonds' in self.D and 'mol2_bonds' in other.D:
            other.D['mol2_bonds'].bondIdx+=self.D['mol2_bonds'].shape[0]
        self._myconcat(other,directive='mol2_bonds',idxlabel=['ai','aj'],idxshift=idxshift)
        self._myconcat(other,directive='pairs',idxlabel=['ai','aj'],idxshift=idxshift)
        self._myconcat(other,directive='angles',idxlabel=['ai','aj','ak'],idxshift=idxshift)
        self._myconcat(other,directive='dihedrals',idxlabel=['ai','aj','ak','al'],idxshift=idxshift)

    def get_atom(self,idx):
        return self.D['atoms'].iloc[idx-1]

    def get_atomtype(self,idx):
#        logging.debug(f'Asking get_atomtype for type of atom with index {idx}')
        return self.D['atoms'].iloc[idx-1].type


