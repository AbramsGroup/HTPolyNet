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

    def add_bonds(self,bondlist=[]):
        #
        pass

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
