import pandas as pd

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
        # allowable keys include all keys in _GromacsTopologyDirectiveHeaders_ plus "impropers"
        self.bondlist={}

    @classmethod
    def from_topfile(cls,filename):
        inst=cls()
        '''
        Reads a Gromacs-style topology file 'filename' and returns a dictionary keyed on directive names.
        Each value in the dictionary is a list containing one or more pandas dataframes.  Each
        dataframe represents an individual section found with its directive in the file.  Gromacs
        directives can be repeated in a file (well, *some* of them can), so this is the most
        general solution.
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
                # This is going to be inconvenient; most dictionary entries will have only one dataframe
                # 'dihedrals' can have two
                if not directive in inst.D:
                    inst.D[directive]=[]
                inst.D[directive].append(tdf)
            # TODO: create bondlist
            return inst

    def df(self,directive):
        if directive in self.D:
            if directive=='impropers':
                rdf=self.D['dihedrals'][1]
            elif directive=='dihedrals':
                rdf=self.D['dihedrals'][0]
            else:
                rdf=self.D[directive][0]
            return rdf
        return None

    def __str__(self):
        ''' Generates a string in the proper top format '''
        retstr=''
        for k,v in self.D.items():
            for vv in v:
                if vv.empty:
                    ''' an empty stanza will not be output '''
                    continue
                retstr+='[ '+k+' ]\n'
                retstr+='; '+'\t'.join(vv.columns)+'\n'
                for i,row in vv.iterrows():
                    ''' assumes NaN's are only allowed in trailing columns '''
                    retstr+='\t'.join(row[~row.isna()])+'\n'
        return retstr

    def add_bonds(self,bondlist=[]):
        # 
        pass

    def delete_atoms(self,idx=[],reindex=True):
        d=self.D['atoms'][0]
        indexes_to_drop=d[d.nr.isin(idx)].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['atoms'][0]=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['atoms'][0]
            oldGI=d['nr'].copy()
            d['old_nr']=oldGI
            d['nr']=d.index+1
            mapper={k:v for k,v in zip(d['old_nr'],d['nr'])}
            d['nr_shift']=d['nr']-oldGI  # probably not necessary
        d=self.D['bonds'][0]
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['bonds'][0]=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['bonds'][0]
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
        d=self.D['angles'][0]
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))|(d.ak.isin(idx))].index
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['angles'][0]=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['angles'][0]
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)      
            d.ak=d.ak.map(mapper)
        # TODO: dihedrals
        # TODO: update bondlist using idx

    def merge_top(self,other):

        # TODO: must update bondlist
        pass

    def make_bondlist(self):
        self.bondlist={}

        pass

    def delete_from_bondlist(self,idx=[]):
        # TODO: delete (ai in idx)'s entry in bondlist
        # and remove all instances of (ai in idx) in other atoms' bondlists
        pass