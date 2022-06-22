"""

.. module:: topology
   :synopsis: Class for managing gromacs .top file data
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""

from hashlib import new
import pandas as pd
import logging
from HTPolyNet.bondlist import Bondlist
import os
from copy import deepcopy
from scipy.constants import physical_constants
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import matplotlib.pyplot as plt
import json
from HTPolyNet.plot import network_graph

def typeorder(a):
    """typeorder correctly order the tuple of atom types for particular
        interaction types to maintain sorted type dataframes

    :param a: tuple of atom indicies/types from a [ bond ], [ pair ], [ angle ], or [ dihedral ] record
    :type a: tuple
    :return: same atom indices/types correctly ordered to allow for easy searching/sorting
    :rtype: tuple
    """
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

def repeat_check(t,msg=''):
    """repeat_check Check for repeated index tuples

    :param t: list of index tuples
    :type t: list
    :param msg: optional message, defaults to ''
    :type msg: str, optional
    """
    for i in range(len(t)):
        for j in range(i+1,len(t)):
            assert t[i]!=t[j],f'Error: repeated index in {len(t)}-tuple {t}: t({i})={t[i]}\n{msg}'

def df_typeorder(df,typs):
    """df_typeorder type-orders the atom type attributes in each row of dataframe df

    :param df: a Topology type-directive dataframe; [ atomtypes ], [ bondtypes ], etc.
    :type df: pandas.DataFrame
    :param typs: list of type-attribute names; typically ['i','j',...]
    :type typs: list
    """
    for i in df.index:
        df.loc[i,typs]=typeorder(tuple(df.loc[i,typs]))

def treadmill(L):
    """treadmill Move first element of list L to end and return new list

    :param L: a list
    :type L: list
    :return: a new list
    :rtype: list
    """
    nL=L[1:]
    nL.append(L[0])
    return nL

def treadmills(L):
    """treadmills perform one complete cycle of treadmilling increments and store each increment as its own list and return list of such lists

    :param L: a list
    :type L: list
    :return: list of new lists, each a treadmill increment of passed-in list
    :rtype: list of lists
    """
    N=len(L)
    nL=L
    r=[]
    for i in range(1,N):
        nnL=treadmill(nL)
        r.append(nnL)
        nL=nnL
    return r

def _get_unique_cycles_dict(G,min_length=-1):
    ucycles={}
    counts_by_length={}
    for u in nx.simple_cycles(G):
        sl=len(u)
        if min_length<=sl:
            logging.debug(f'a cycle {u}')
            if not sl in counts_by_length:
                counts_by_length[sl]=0
            counts_by_length[sl]+=1
            if not sl in ucycles:
                ucycles[sl]=[]
            #u=[x+1 for x in u] # why am i doing this?
            utl=treadmills(u)
            ur=list(reversed(u))
            urtl=treadmills(ur)
            eqv=utl+[ur]+urtl
            ll=min_length if min_length!=-1 else min(counts_by_length.keys())
            uu=max(counts_by_length.keys())
            for l in range(ll,uu+1):
                if len(u)==l:
                    found=False
                    for e in eqv:
                        if e in ucycles[l]:
                            found=True
                            break
                    if not found:
                        ucycles[l].append(u)
    return ucycles

def _present_and_contiguous(subL,L):
    """_present_and_contiguous returns True is elements in subL appear as a contiguous sub-block in L

    :param subL: a sublist
    :type subL: list
    :param L: a list
    :type L: list
    :param forward_and_reverse: _description_, defaults to True
    :type forward_and_reverse: bool, optional
    :param periodic: _description_, defaults to True
    :type periodic: bool, optional
    :return: True or false
    :rtype: boolean
    """
    pretest=all([x in L for x in subL])
    if not pretest:
        return False
    for A in [L,L[::-1]]:
        T=treadmills(A)
        for t in T:
            testL=t[:len(subL)]
            logging.debug(f'___ {subL} {testL}')
            if all([x==y for x,y in zip(subL,testL)]):
                return True
    return False

_GromacsIntegers_=('nr','atnum','resnr','ai','aj','ak','al','#mols','nrexcl','funct','func','nbfunc','comb-rule')
_GromacsFloats_=('charge','mass','chargeB','massB',*tuple([f'c{i}' for i in range(5)]),
                 'b0','kb','th0','cth','rub','kub','phase','kd','pn','fudgeLJ','fudgeQQ')
def typedata(h,s):
    if h in _GromacsIntegers_:
        return int(s)
    if h in _GromacsFloats_:
        return float(s)
    return s

_GromacsExtensiveDirectives_=['atoms','pairs','bonds','angles','dihedrals']
_NonGromacsExtensiveDirectives_=['mol2_bonds']
_GromacsTopologyDirectiveOrder_=['defaults','atomtypes','bondtypes','angletypes','dihedraltypes','moleculetype','atoms','pairs','bonds','angles','dihedrals','system','molecules']
_GromacsTopologyDirectiveHeaders_={
    'atoms':['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass','typeB', 'chargeB', 'massB'],
    'pairs':['ai', 'aj', 'funct', 'c0', 'c1'],
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
# _GromacsTopologyDirective_ExtraAttributes_={
#     'bonds':['needs_relaxing'],
#     'angles':['needs_relaxing'],
#     'dihedrals':['needs_relaxing']
# }
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
# _GromacsTopologyDataFields_={
#     'atoms':['type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass','typeB', 'chargeB', 'massB'],
#     'pairs':['funct', 'c0', 'c1'],
#     'bonds':['funct', 'c0', 'c1'],
#     'angles':['funct', 'c0', 'c1'],
#     'dihedrals':['funct', 'c0', 'c1', 'c2', 'c3', 'c4', 'c5'],
#     'atomtypes':['atnum', 'mass', 'charge', 'ptype', 'sigma', 'epsilon'],
#     'bondtypes':['func','b0','kb'],
#     'angletypes':['func','th0','cth','rub','kub'],
#     'dihedraltypes':['func','phase','kd','pn']
#     }
_GromacsTopologyDirectiveDefaults_={
    'system':['A_generic_system'],
    'molecules':['None',1],
    'moleculetype':['None',3],
#    'defaults':[1,2,'no',0.5,0.83333333]
    'defaults':[1,2,'yes',0.5,0.83333333]
}

class Topology:
    """ Class for handling gromacs top data
    """
    def __init__(self,system_name=''):
        """__init__ Constructor for Topology class

        :param system_name: optional name of system, defaults to ''
        :type system_name: str, optional
        """
        ''' D: a dictionay keyed on Gromacs topology directives with values that are lists of
               one or more pandas dataframes corresponding to sections '''
        self.D={}
        for k,v in _GromacsTopologyDirectiveDefaults_.items():
            hdr=_GromacsTopologyDirectiveHeaders_[k]
            dfdict={kk:[a] for kk,a in zip(hdr,v)}
            self.D[k]=pd.DataFrame(dfdict)
        self.D['system']=pd.DataFrame({'name':[system_name]})
        ''' bondlist: a class that owns a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist=Bondlist()
        self.residue_network=nx.DiGraph()
        self.polyethylenes=nx.DiGraph()
        self.empty=True

    @classmethod
    def read_gro(cls,filename):
        """read_gro Reads a Gromacs-style topology file 'filename' and returns a dictionary keyed on directive names.
        Each value in the dictionary is a pandas dataframe.  Each
        dataframe represents an individual section found with its directive in the file, with columns corresponding to the fields in the section.  Note that the we allow for input topology/itp files to have two 'dihedrals' and 'dihedraltypes' sections; these
        are merged in the result.

        :param filename: name of gromacs top file to read
        :type filename: str
        :raises KeyError: If an unrecognized topology directive is encountered, program exits on error
        :return: a Topology instance
        :rtype: Topology
        """
        assert os.path.exists(filename), f'Error: {filename} not found.'
        inst=cls()
        inst.filename=filename
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
            inst.null_check(msg=f'read from {filename}')
            if 'atomtypes' in inst.D:
                inst.D['atomtypes'].sort_values(by='name',inplace=True)
            if 'bonds' in inst.D:
                inst.D['bonds']=inst.D['bonds'].sort_values(by=['ai','aj']).reset_index(drop=True)
                inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
            if 'bondtypes' in inst.D:
                df_typeorder(inst.D['bondtypes'],typs=['i','j'])
                inst.D['bondtypes'].sort_values(by=['i','j'],inplace=True)
            if 'pairs' in inst.D:
                inst.D['pairs']=inst.D['pairs'].sort_values(by=['ai','aj']).reset_index(drop=True)
            if 'pairtypes' in inst.D:
                df_typeorder(inst.D['pairtypes'],typs=['i','j'])
                inst.D['pairtypes'].sort_values(by=['i','j'],inplace=True)
            if 'angles' in inst.D:
                # central atom (aj) is the primary index
                inst.D['angles']=inst.D['angles'].sort_values(by=['aj','ai','ak']).reset_index(drop=True)
            if 'angletypes' in inst.D:
                df_typeorder(inst.D['angletypes'],typs=['i','j','k'])
                inst.D['angletypes'].sort_values(by=['j','i','k'],inplace=True)
            if 'dihedrals' in inst.D:
                # central atoms (aj,ak) are the primary index
                inst.D['dihedrals']=inst.D['dihedrals'].sort_values(by=['aj','ak','ai','al']).reset_index(drop=True)
            if 'dihedraltypes' in inst.D:
                df_typeorder(inst.D['dihedraltypes'],typs=['i','j','k','l'])
                # print(f'    -> pre sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
                inst.D['dihedraltypes'].sort_values(by=['j','k','i','l'],inplace=True)
                # print(f'    -> post sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
            for f in inst.includes:
                # print(f'reading included topology {f}')
                inst.merge(Topology.read_gro(f))
            inst.empty=False
            return inst

    def bond_source_check(self):
        """bond_source_check Checks to ensure the 'bonds' dataframe and 'mol2_bonds' dataframe contain the same bonds.  A mol2 dataframe is only created when a mol2 file is read by the Coordinates module.
        """
        if 'bonds' in self.D and 'mol2_bonds' in self.D:
            logging.info(f'Consistency check between gromacs-top bonds and mol2-bonds requested.')
            grobonds=self.D['bonds'].sort_values(by=['ai','aj'])
            bmi=grobonds.set_index(['ai','aj']).index
            mol2bonds=self.D['mol2_bonds'].sort_values(by=['ai','aj'])
            mbmi=mol2bonds.set_index(['ai','aj']).index
            check=all([x==y for x,y in zip(bmi,mbmi)])
            logging.info(f'Result: {check}')
            if not check:
                logging.info(f'GROMACS:')
                logging.info(grobonds[['ai','aj']].head().to_string())
                logging.info(f'MOL2:')
                logging.info(mol2bonds[['ai','aj']].head().to_string())
                for x,y in zip(bmi,mbmi):
                    logging.info(f'{x} {y} {x==y}')

    # def has_bond(self,pair):
    #     """has_bond Determines whether or not the bond represented by the two atom indices in 
    #     pair exists in the Topology

    #     :param pair: two atom indexes
    #     :type pair: tuple
    #     :return: True if bond exists in both the 'bonds' and 'mol2_bonds' dataframes
    #     :rtype: boolean
    #     """
    #     bmi=self.D['bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
    #     mbmi=self.D['mol2_bonds'].sort_values(by=['ai','aj']).set_index(['ai','aj']).index
    #     return pair in bmi and pair in mbmi

    def shiftatomsidx(self,idxshift,directive,rows=[],idxlabels=[]):
        """shiftatomsidx shifts all atoms indexes in topology directive dataframe

        :param idxshift: integer index shift
        :type idxshift: int
        :param directive: name of gromacs topology directive ('atoms','bonds','pairs','angles','dihedrals')
        :type directive: string
        :param rows: row boundaries, defaults to []
        :type rows: list, optional
        :param idxlabels: names of columns that contain atom indexes, defaults to []
        :type idxlabels: list, optional
        """
        if directive in self.D:
            cols=self.D[directive].columns.get_indexer(idxlabels)
            self.D[directive].iloc[rows[0]:rows[1],cols]+=idxshift

    def ring_detector(self):
        """ring_detector Detects all 3 to 7-members rings in a topology for one residue using the bondlist and networkx
        sets the Cycles attribute of self to ring-lists of atom names

        :return: A dictionary of ring-lists keyed on ring size; each ring-list is a list of atom indexes
        :rtype: dict
        """
        g=self.bondlist.graph()
        cycles=_get_unique_cycles_dict(g,min_length=3)
        cycles_by_name={}
        for a in cycles:
            cycles_by_name[a]=[]
            for c in cycles[a]:
                atoms=[]
                for atom in c:
                    atoms.append(self.get_atom_attribute(atom,'atom'))
                cycles_by_name[a].append(atoms)
            if len(cycles_by_name[a])>0:
                logging.debug(f'{a}-rings: {cycles_by_name[a]}')
        self.Cycles=cycles_by_name
        return cycles

    def rep_ex(self,count=0):
        """Replicate extensive topology components (atoms, pairs, bonds, angles, dihedrals)

        :param count: number of replicas to generate, defaults to 0
        :type count: int, optional
        :raises Exception: Dies if self is missing an atoms dataframe
        """
        if count>0:
            counts={k:0 for k in _GromacsExtensiveDirectives_}
            counts.update({k:0 for k in _NonGromacsExtensiveDirectives_})
            for t in _GromacsExtensiveDirectives_:
                if t in self.D:
                    counts[t]=len(self.D[t])
                for t in _NonGromacsExtensiveDirectives_:
                    if t in self.D:
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
        """from_ex make a new Topology instance by copying only the extensive dataframes
            from an existing topology 

        :param other: the other topology
        :type other: Topology
        :return: a new Topology generated by the extensive dataframes of other
        :rtype: Topology
        """
        ''' '''
        inst=cls()
        for t in _GromacsExtensiveDirectives_:
            if t in other.D:
                inst.D[t]=other.D[t].copy()
        return inst

    def to_file(self,filename):
        """to_file Write topology to a gromacs-format file

        :param filename: name of top file to write
        :type filename: str
        """
        # prevent buggy writing of NaNs
        self.null_check(msg=f'writing {filename}')
        with open(filename,'w') as f:
            f.write('; Gromacs-format topology written by HTPolyNet\n')
        assert 'defaults' in self.D, 'Error: no [ defaults ] in topology?'
        for k in _GromacsTopologyDirectiveOrder_:
            if k in self.D:
                # columns=_GromacsTopologyDirectiveHeaders_[k]
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

    def null_check(self,msg=''):
        """Paranoid checking for NaNs in dataframe locations that SHOULD NEVER HAVE NANS

        :param msg: a nice message, defaults to ''
        :type msg: str, optional
        :raises Exception: exits if a NaN is found
        """
        check = True
        for k in _GromacsTopologyDirectiveOrder_:
            if k in self.D:
                if k in _GromacsTopologyHashables_:
                    for a in _GromacsTopologyHashables_[k]:
                        check=any(self.D[k][a].isnull())
                        if check:
                            logging.debug(f'{msg} null in {k} {a}\n{self.D[k].to_string()}')
                            raise Exception('NaN error')

    def total_charge(self):
        """Compute and return total system charge

        :return: charge
        :rtype: float
        """
        if 'atoms' in self.D:
            return self.D['atoms']['charge'].sum()
        return 0.0

    def adjust_charges(self,desired_charge=0.0,msg=''):
        """Adjust atom partial charges a tiny bit so that total system charge is zero

        :param desired_charge: target system charge, defaults to 0.0
        :type desired_charge: float, optional
        :param msg: A message to write if pre-adjusted system charge is too high, defaults to ''
        :type msg: str, optional
        :return: self topology
        :rtype: Topology
        """
        apparent_charge=self.total_charge()
        overcharge=apparent_charge-desired_charge
        logging.info(f'Adjusting charges due to overcharge of {overcharge}')
        if np.abs(overcharge)>0.1:
            logging.info(f'{msg}')
        cpa=-overcharge/len(self.D['atoms'])
        self.D['atoms']['charge']+=cpa
        logging.info(f'New total charge after adjustment: {self.total_charge()}')
        return self
        
    def total_mass(self,units='gromacs'):
        """Returns total mass of all atoms in the Topology.

        :param units: unit system designation; if 'SI' returns kg, defaults to 'gromacs'
        :type units: str, optional
        :return: mass (in amu if units is 'gromacs' or kg if units is 'SI')
        :rtype: float
        """
        fac=1.0
        if units=='SI':
            fac=physical_constants['atomic mass constant'][0]
        if 'atoms' in self.D:
            M_amu=self.D['atoms']['mass'].sum()
            return M_amu*fac
        return 0.0

    def atomcount(self):
        """atomcount Returns the total number of atoms

        :return: number of atoms
        :rtype: int
        """
        if 'atoms' in self.D:
            return len(self.D['atoms'])
        return 0

    # def add_pairs(self,pairsdf,kb=280160.0):
    #     pmi=self.D['pairs'].set_index(['ai','aj']).sort_index().index
    #     for i,p in pairsdf.iterrows():
    #         ai,aj=idxorder((p['ai'],p['aj']))
    #         l0=p['initial-distance']
    #         if not (ai,aj) in pmi: #this pair not already here, good!
    #             data=[ai,aj,1,l0,kb]
    #             h=_GromacsTopologyDirectiveHeaders_['pairs']
    #             pairdict={k:[v] for k,v in zip(h,data)}
    #             pairtoadd=pd.DataFrame(pairdict)
    #             self.D['pairs']=pd.concat((self.D['pairs'],pairtoadd),ignore_index=True)
    #         else:
    #             logging.debug(f'Warning: pair {ai}-{aj} already in [ pairs ].  This is bug.')
    
    def add_restraints(self,pairdf,typ=6,kb=300000.):
        """Add type-6 (non-topoogical) bonds to help drag atoms destined to be bonded
        closer together in a series of dragging simulations

        :param pairdf: dataframe of pairs ['ai','aj']
        :type pairdf: pandas DataFrame
        :param typ: bond type, defaults to 6
        :type typ: int, optional
        :param kb: bond spring constant (kJ/mol/nm^2), defaults to 300000
        :type kb: float, optional
        """
        bmi=self.D['bonds'].set_index(['ai','aj']).sort_index().index
        for i,b in pairdf.iterrows():
            ai,aj=idxorder((b['ai'],b['aj']))
            b0=b['initial-distance']
            if not (ai,aj) in bmi:
                h=_GromacsTopologyDirectiveHeaders_['bonds']
                data=[ai,aj,typ,b0,kb]  # this new bond will have override parameters
                bonddict={k:[v] for k,v in zip(h,data)}
                bdtoadd=pd.DataFrame(bonddict)
                self.D['bonds']=pd.concat((self.D['bonds'],bdtoadd),ignore_index=True)

    def remove_restraints(self,pairdf):
        """Remove all bonds represented in in pairdf.
        These are interpreted as non-topological
        restraints, so deleting these 'bonds' does 
        not influence angles or dihedrals

        :param pairdf: dataframe of pairs ['ai','aj']
        :type pairdf: pandas DataFrame
        """
        d=self.D['bonds']
        to_drop=[]
        for i,b in pairdf.iterrows():
            ai,aj=idxorder((b['ai'],b['aj']))
            to_drop.append(d[(d.ai==ai)&(d.aj==aj)].index.values[0])
        self.D['bonds']=self.D['bonds'].drop(to_drop)


    def update_polyethylenes(self,bdf,idx_mapper):
        self.polyethylenes=nx.relabel_nodes(self.polyethylenes,idx_mapper)
        for i,r in bdf.iterrows():
            a=r['ai']
            b=r['aj']
            self.polyethylenes.add_edge(a,b)
    
    def polyethylene_cycle(self,i,j):
        """polyethylene_cycle return True if adding edge i<->j generates a cycle in
           self.polythelenes

        :param i: one atom index
        :type i: int
        :param j: another atom index
        :type j: int
        """
        precycles=list(nx.simple_cycles(self.polyethylenes))
        makes_a_cycle=False
        clens={}
        for c in precycles:
            l=len(c)
            if not l in clens:
                clens[l]=0
            clens[l]+=1
        makes_a_cycle=any([l>2 for l in clens.keys()])
        if makes_a_cycle:
            logging.debug(f'there is a problem with the polyethylene cycles!')
            network_graph(self.polyethylenes,'pe_net_bad.png')
        assert not makes_a_cycle

        self.polyethylenes.add_edge(i,j)
        cycles=list(nx.simple_cycles(self.polyethylenes))
        makes_a_cycle=False
        clens={}
        for c in cycles:
            l=len(c)
            if not l in clens:
                clens[l]=0
            clens[l]+=1
        makes_a_cycle=any([l>2 for l in clens])
        self.polyethylenes.remove_edge(i,j)
        return makes_a_cycle

    def polyethylene_cycles_collective(self,B):
        """polyethylene_cycles_collective look for cycles that appear when two or more proposed bonds are created

        :param B: list of proposed bonds as "passbond" instances
        :type B: list of bondrecordtuples (ai,aj,rij)
        """
        '''
                self.bond=bondtuple
                self.reactantname=reactantname
                self.distance=distance
                self.probability=probability
                self.order=order
        
        '''

        for b in B:
            ai,aj=b.bond
            self.polyethylenes.add_edge(ai,aj)
            assert ai in self.polyethylenes.nodes
            assert aj in self.polyethylenes.nodes
        bad_cycles=[]
        bad_cycle_dict=_get_unique_cycles_dict(self.polyethylenes,min_length=4)
        logging.debug(f'bad_cycle_dict: {bad_cycle_dict}')
        for k,v in bad_cycle_dict.items():
            bad_cycles.extend(v)

        if len(bad_cycles)>0:
            logging.debug(f'Proposed bondset forms {len(bad_cycles)} disallowed cycles')
            for b in B:
                ai,aj=b.bond
                logging.debug(f' -- b {ai} - {aj}')
            for c in bad_cycles:
                logging.debug(f' || c {c}')
        
        bad_bonds=[]
        for b in B:
            ai,aj=b.bond
            # logging.debug(f'collective pe cycle testing bond {ai} {aj}')
            self.polyethylenes.remove_edge(ai,aj)

            # is this pair in any of the cycles?
            target_cycles=[]
            for bc in bad_cycles:
                if _present_and_contiguous([ai,aj],bc):
                    logging.debug(f'Proposed bond {ai}-{aj} lives in cycle {bc}')
                    target_cycles.append(bc)
                    bad_bonds.append(b)
                    bad_cycles.remove(bc)
                    break
        if len(bad_bonds)>1:
            logging.debug(f'Proposed bonds to remove: {bad_bonds}')
            for b in bad_bonds:
                B.remove(b)
        return B
        
    def add_bonds(self,pairs=[]):
        """add_bonds Adds bonds indicated in list pairs to the topology

        :param pairs: list of pairs of atom indexes, defaults to []
        :type pairs: list, optional
        :raises Exception: dies if an existing bond is in the list of pairs
        """
        # logging.debug('add_bonds begins')
        at=self.D['atoms']
        ij=self.D['bondtypes'].set_index(['i','j'])
        #mb=self.D['mol2_bonds']
        bmi=self.D['bonds'].set_index(['ai','aj']).sort_index().index
        pmi=self.D['pairs'].set_index(['ai','aj']).sort_index().index
        newbonds=[]
        for b in pairs:
            # logging.debug(f'{b}')
            bondtuple=(b[0],b[1])
            order=b[2]
            ai,aj=idxorder(bondtuple)
            # if this bond is not in the topology
            if not (ai,aj) in bmi:
                newbonds.append((ai,aj))
                # logging.debug(f'asking types of {ai} and {aj}; at.shape {at.shape}')
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                idx=typeorder((it,jt))
                if idx in ij.index:
                    bt=ij.loc[idx,'func']  # why don't i need need values[0]
                    kb=ij.loc[idx,'kb']
                    b0=ij.loc[idx,'b0']
                else:
                    logging.debug(f'no bondtype {idx} found\nI assume you are just making a mol2 file for template parameterization.')
                    bt=1
                    b0=1.5
                    kb=999999
                    # raise Exception(f'no bondtype {idx} found.')
                # add a new bond!
                h=_GromacsTopologyDirectiveHeaders_['bonds']
                data=[ai,aj,bt,b0,kb]  # this new bond will have override parameters
                assert len(h)==len(data), 'Error: not enough data for new bond?'
                bonddict={k:[v] for k,v in zip(h,data)}
                bdtoadd=pd.DataFrame(bonddict)
                self.D['bonds']=pd.concat((self.D['bonds'],bdtoadd),ignore_index=True)
                # logging.info(f'add_bond:\n{bdtoadd.to_string()}')
                # logging.info(f'just added {bonddict}')
                if 'mol2_bonds' in self.D:
                    data=[len(self.D['mol2_bonds']),ai,aj,1] # assume single bond
                    bonddict={k:[v] for k,v in zip(['bondIdx','ai','aj','type'],data)}
                    self.D['mol2_bonds']=pd.concat((self.D['mol2_bonds'],pd.DataFrame(bonddict)),ignore_index=True)
                # remove this pair from pairs if it's in there (it won't be)
                if idx in pmi:
                    logging.debug(f'Warning: new bond {ai}-{aj} was evidently in the [ pairs ]!')
                    d=self.D['pairs']
                    indexes_to_drop=d[(d.ai.isin(idx))&(d.aj.isin(idx))].index
                    logging.debug(f'Dropping [ pair ]:\n{d[indexes_to_drop].to_string()}')
                    indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                    self.D['pairs']=d.take(list(indexes_to_keep)).reset_index(drop=True)
            else:
                # TODO: need to allow for possibility of converting an existing single bond to a double bond
                def dec_order_atom_name_gaff(nm):
                    il=list(nm)
                    logging.debug(f'expect this? {il}')
                    elnm=''
                    nsfx=''
                    for i in range(len(il)):
                        logging.debug(f'query {il[i]} {il[i].isdigit()}')
                        if not il[i].isdigit():
                            elnm+=il[i]
                        else:
                            nsfx+=il[i]
                    logging.debug(f'and this? {elnm} {nsfx}')
                    n=int(nsfx)
                    n-=1
                    nn=elnm+str(n)
                    return nn

                ityp=at.loc[ai-1]['type']
                jtyp=at.loc[aj-1]['type']
                logging.debug(f'Need to set order of {ai}({ityp})-{aj}({jtyp}) to {order}')
                # nityp=dec_order_atom_name_gaff(ityp)
                # njtyp=dec_order_atom_name_gaff(jtyp)
                # logging.debug(f'Trying {ai}:{ityp}->{nityp} and {aj}{jtyp}->{njtyp}')
                # at.loc[ai-1,'type']=nityp
                # at.loc[aj-1,'type']=njtyp
                # logging.debug(f'Updated topology [ atoms ]\n:{at.to_string()}')
                if 'mol2_bonds' in self.D:
                    mb=self.D['mol2_bonds']
                    bi=(mb['ai']==ai)&(mb['aj']==aj)
                    mb.loc[bi,'type']=order
                else:
                    logging.warning(f'No way to update this bond since there is not a mol2_bonds attribute in the Topology.')
                    # logging.debug(f'Updated mol2_bonds:\n{mb.to_string()}')
                # raise Exception(f'attempt to add already existing bond {ai}-{aj}')
        # update the bondlist
        for b in newbonds:
            self.bondlist.append(b)
        # logging.debug(f'Added {len(newbonds)} new bonds')

    def add_enumerated_angles(self,newbonds,ignores=[],quiet=True):       
        at=self.D['atoms']
        ijk=self.D['angletypes'].set_index(['i','j','k']).sort_index()
        newangles=[]
        for b in newbonds:
            ''' new angles due to other neighbors of b[0] '''
            for ai in [i for i in self.bondlist.partners_of(b[0]) if (i!=b[1] and not i in ignores)]:
                aj=b[0]
                ak=b[1]
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                kt=at.iloc[ak-1].type
                idx=typeorder((it,jt,kt))
                i,j,k=idxorder((ai,aj,ak))
                repeat_check((i,j,k))
                if idx in ijk.index:
                    angletype=ijk.loc[idx,'func']  # why no .values[0]
                else:
                    if not quiet:
                        logging.warning(f'Angle type {idx} ({ai}-{aj}-{ak}) not found.')
                    angletype=1
                h=_GromacsTopologyDirectiveHeaders_['angles']
                data=[i,j,k,angletype,pd.NA,pd.NA]
                assert len(h)==len(data), 'Error: not enough data for new angle?'
                angledict={k:[v] for k,v in zip(h,data)}
                self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                newangles.append([i,j,k])
            ''' new angles due to other neighbors of b[1] '''
            for ak in [k for k in self.bondlist.partners_of(b[1]) if (k!=b[0] and not k in ignores)]:
                ai=b[0]
                aj=b[1]
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                kt=at.iloc[ak-1].type
                idx=typeorder((it,jt,kt))
                i,j,k=idxorder((ai,aj,ak))
                repeat_check((i,j,k))
                if idx in ijk.index:
                    angletype=ijk.loc[idx,'func'] # why no .values[0]
                else:
                    if not quiet:
                        logging.warning(f'Angle type {idx} ({ai}-{aj}-{ak}) not found.')
                    angletype=1
                h=_GromacsTopologyDirectiveHeaders_['angles']
                data=[i,j,k,angletype,pd.NA,pd.NA]
                assert len(h)==len(data), 'Error: not enough data for new angle?'
                angledict={k:[v] for k,v in zip(h,data)}
                self.D['angles']=pd.concat((self.D['angles'],pd.DataFrame(angledict)),ignore_index=True)
                newangles.append([i,j,k])
        return newangles
        
    def add_enumerated_dihedrals(self,newbonds,ignores=[],quiet=True):
        newdihedrals=[]
        newpairs=[]
        at=self.D['atoms']
        ijkl=self.D['dihedraltypes'].set_index(['i','j','k','l']).sort_index()

        ''' new proper dihedrals for which the new bond is the central j-k bond '''
        for b in newbonds:
            aj,ak=idxorder(b)
            for ai in [i for i in self.bondlist.partners_of(aj) if (i!=ak and not i in ignores)]:
                for al in [l for l in self.bondlist.partners_of(ak) if (l!=aj and not l in ignores)]:
                    it=at.iloc[ai-1].type
                    jt=at.iloc[aj-1].type
                    kt=at.iloc[ak-1].type
                    lt=at.iloc[al-1].type
                    idx=typeorder((it,jt,kt,lt))
                    i,j,k,l=idxorder((ai,aj,ak,al))
                    repeat_check((i,j,k,l),msg=f'central {j}-{k}')
                    if idx in ijkl.index:
                        dihedtype=ijkl.loc[idx,'func'].values[0] # why values[0]
                    else:
                        if not quiet:
                            logging.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                        dihedtype=9
                    h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                    data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                    assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                    diheddict={k:[v] for k,v in zip(h,data)}
                    self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                    newdihedrals.append([i,j,k,l])
                    ''' i-l is a new 1-4 pair '''
                    h=_GromacsTopologyDirectiveHeaders_['pairs']
                    data=[i,l,1]
                    pairdict={k:[v] for k,v in zip(h,data)}
                    self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                    newpairs.append((i,l))
            ''' new proper dihedrals for which the new bond is the i-j or j-i bond '''
            for ai,aj in zip(b,reversed(b)):
                for ak in [k for k in self.bondlist.partners_of(aj) if (k!=ai and not k in ignores)]:
                    for al in [l for l in self.bondlist.partners_of(ak) if (l!=aj and not l in ignores)]:
                        it=at.iloc[ai-1].type
                        jt=at.iloc[aj-1].type
                        kt=at.iloc[ak-1].type
                        lt=at.iloc[al-1].type
                        idx=typeorder((it,jt,kt,lt))
                        i,j,k,l=idxorder((ai,aj,ak,al))
                        repeat_check((i,j,k,l),msg=f'i-j neighbor of j-k {j}-{k}')
                        if idx in ijkl.index:
                            dihedtype=ijkl.loc[idx,'func'].values[0]
                            # dihedtype=ijkl.loc[idx,'func']
                        else:
                            if not quiet:
                                logging.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                            dihedtype=9
                        h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                        data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                        assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                        diheddict={k:[v] for k,v in zip(h,data)}
                        self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                        newdihedrals.append([i,j,k,l])
                        ''' i-l is a new 1-4 pair '''
                        h=_GromacsTopologyDirectiveHeaders_['pairs']
                        data=[i,l,1]
                        pairdict={k:[v] for k,v in zip(h,data)}
                        self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                        newpairs.append((i,l))

            ''' new proper dihedrals for which the new bond is the k-l or l-k bond '''
            for ak,al in zip(b,reversed(b)):
                for aj in [j for j in self.bondlist.partners_of(ak) if (j!=al and not j in ignores)]:
                    for ai in [i for i in self.bondlist.partners_of(aj) if (i!=ak and not i in ignores)]:
                        it=at.iloc[ai-1].type
                        jt=at.iloc[aj-1].type
                        kt=at.iloc[ak-1].type
                        lt=at.iloc[al-1].type
                        idx=typeorder((it,jt,kt,lt))
                        i,j,k,l=idxorder((ai,aj,ak,al))
                        repeat_check((i,j,k,l),msg=f'k-l neighbor of j-k {j}-{k}')
                        if idx in ijkl.index:
                            dihedtype=ijkl.loc[idx,'func'].values[0]
                        else:
                            if not quiet:
                                logging.warning(f'Dihedral type {idx} {ai}-{aj}-{ak}-{al} not found.')
                            dihedtype=9
                        h=_GromacsTopologyDirectiveHeaders_['dihedrals']
                        data=[i,j,k,l,dihedtype,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA,pd.NA]
                        assert len(h)==len(data), 'Error: not enough data for new  dihedral?'
                        diheddict={k:[v] for k,v in zip(h,data)}
                        self.D['dihedrals']=pd.concat((self.D['dihedrals'],pd.DataFrame(diheddict)),ignore_index=True)
                        newdihedrals.append([i,j,k,l])
                        h=_GromacsTopologyDirectiveHeaders_['pairs']
                        ''' i-l is a new 1-4 pair '''
                        data=[i,l,1]
                        pairdict={k:[v] for k,v in zip(h,data)}
                        self.D['pairs']=pd.concat((self.D['pairs'],pd.DataFrame(pairdict)),ignore_index=True)
                        newpairs.append((i,l))
        return newdihedrals,newpairs

    def delete_atoms(self,idx=[],reindex=True,return_idx_of=[]):
        """Delete atoms from topology

        :param idx: list of atom indexes to delete, defaults to []
        :type idx: list, optional
        :param reindex: reindex atoms after deleting, defaults to True
        :type reindex: bool, optional
        :param return_idx_of: list of old indices to report new indices of, defaults to []
        :type return_idx_of: list, optional
        :return: old-index-to-new-index mapper
        :rtype: dict
        """
        #logging.debug(f'Delete atoms: {idx}')
        self.null_check(msg='beginning of delete atoms')
        # logging.debug(f'idx {idx}')
        d=self.D['atoms']
        new_idx=[]
        indexes_to_drop=d[d.nr.isin(idx)].index
        logging.info(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ atoms ]')#:\n{d.loc[indexes_to_drop].to_string()}')
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['atoms']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        mapper={}
        if reindex:
            d=self.D['atoms']
            oldGI=d['nr'].copy()
            d['nr']=d.index+1
            mapper={k:v for k,v in zip(oldGI,d['nr'])}
            # logging.debug(f'mapper {mapper}')
            assert not any([x in mapper for x in idx]),f'Error: Some deleted atoms in mapper.'
            k=np.array(list(mapper.keys()))
            v=np.array(list(mapper.values()))
            if any(np.isnan(k)):
                logging.error('null in mapper keys')
            if any(np.isnan(v)):
                logging.error('null in mapper values')
            # logging.debug(f'delete_atoms: mapper {mapper}')
            if len(return_idx_of)>0:
                # logging.info(f'Asking for updated global indexes of {return_idx_of}')
                new_idx=[mapper[o] for o in return_idx_of]
            #d['nr_shift']=d['nr']-oldGI  # probably not necessary
        ptt=['bonds','mol2_bonds','pairs']
        for pt in ptt:
            if pt in self.D:
                d=self.D[pt]
                # logging.debug(f'delete atom: {pt} df prior to deleting')
                # logging.debug(d.to_string())
                indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
                logging.info(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ {pt} ]')#:\n{d.loc[indexes_to_drop].to_string()}')

                # logging.debug(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ {pt} ]')#:\n{d.loc[indexes_to_drop].to_string()}')
                # logging.debug(f'dropping {pt} {indexes_to_drop}')
                indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                self.D[pt]=d.take(list(indexes_to_keep)).reset_index(drop=True)
                if reindex:
                    d=self.D[pt]
                    assert not any([x in idx for x in d.ai]),f'Error: deleted atom survived in {pt} ai'
                    assert not any([x in idx for x in d.aj]),f'Error: deleted atom survived in {pt} aj'
                    assert all([x in mapper for x in d.ai]),f'Error: surviving {pt} atom ai old idx not in mapper'
                    assert all([x in mapper for x in d.aj]),f'Error: surviving {pt} atom aj old idx not in mapper'
                    #logging.debug(f'delete atom: {pt} df prior to reindexing')
                    #logging.debug(d.to_string())
                    # pairs deleted here were deleted because either ai or aj was among 
                    # the atoms to delete.  We assert than any dihedral for which
                    # the i atom is ai and l atom is aj will necessarily be deleted 
                    # below.  The only other pairs that should be deleted would
                    # be ones in which the dihedral j or k atom is among those to be 
                    # deleted.  
                    if pt!='pairs': # don't remap these yet; need to delete pairs that
                        # might arise from dihedrals that are deleted.
                        d.ai=d.ai.map(mapper)
                        d.aj=d.aj.map(mapper)
                    if pt=='bonds':
                        # logging.debug(f'Updating bondlist using\n{d.to_string()}')
                        self.bondlist=Bondlist.fromDataFrame(d)
                    if pt=='mol2_bonds':
                        nBonds=self.D[pt].shape[0]
                        self.D[pt]['bondIdx']=list(range(1,nBonds+1))
        d=self.D['angles']
        # assert d.ai.dtype==int,f'pre-delete lost angle ai dtype {d.ai.dtype}'
        # assert d.aj.dtype==int,f'pre-delete lost angle aj dtype {d.aj.dtype}'
        # assert d.ak.dtype==int,f'pre-delete lost angle ak dtype {d.ak.dtype}'
        # logging.debug(f'ai {d.ai.isin(idx).to_string()}')
        # logging.debug(f'aj {d.aj.isin(idx).to_string()}')
        # logging.debug(f'ak {d.ak.isin(idx).to_string()}')
        indexes_to_drop=d[(d['ai'].isin(idx))|(d['aj'].isin(idx))|(d['ak'].isin(idx))].index
        # extras=d[d['ak'].isin(idx)].index
        logging.debug(f'Deleting {len(indexes_to_drop)} [ angles ]')#:\n{d.loc[indexes_to_drop].to_string()}')
        # logging.debug(f'any? {list(extras)}')
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        # logging.debug(f'drop {list(sorted(list(set(indexes_to_drop))))}')
        # logging.debug(f'keep {indexes_to_keep}')
        self.D['angles']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        assert d.ai.dtype==int,f'post-delete lost angle ai dtype {d.ai.dtype}'
        assert d.aj.dtype==int,f'post-delete lost angle aj dtype {d.aj.dtype}'
        assert d.ak.dtype==int,f'post-delete lost angle ak dtype {d.ak.dtype}'
        self.null_check(msg='inside delete atoms before angles reindex')
        if reindex:
            d=self.D['angles']
            assert not any([x in idx for x in d.ai]),'Error: deleted atom survived in angle ai'
            assert not any([x in idx for x in d.aj]),'Error: deleted atom survived in angle aj'
            zombie_tags=[x in idx for x in d.ak]
            if any(zombie_tags):
                logging.debug(f'Zombie ak angles:\n{self.D["angles"][zombie_tags].to_string()}')
            assert not any(zombie_tags),'Error: deleted atom survived in angle ak'
            assert all([x in mapper for x in d.ai]),'Error: surviving angle atom ai old idx not in mapper'
            assert all([x in mapper for x in d.aj]),'Error: surviving angle atom aj old idx not in mapper'
            assert all([x in mapper for x in d.ak]),'Error: surviving angle atom ak old idx not in mapper'
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
            d.ak=d.ak.map(mapper)
        self.null_check(msg='inside delete atoms after angles reindex')
        d=self.D['dihedrals']
        indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))|(d.ak.isin(idx))|(d.al.isin(idx))].index
        logging.debug(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ dihedrals ]')#:\n{d.loc[indexes_to_drop].to_string()}')
        # TODO: make sure you delete pairs!!
        dp=self.D['pairs']
        ddp=d.loc[indexes_to_drop]  # these are dihedrals marked for deletion
        # determine pairs deriving from these dihedrals and delete them!
        pai=ddp.ai.to_list()
        pal=ddp.al.to_list()
        dd=[]
        for pi,pl in zip(pai,pal):
            dwpi=dp[((dp.ai==pi)&(dp.aj==pl))|((dp.ai==pl)&(dp.aj==pi))].index.to_list()
            dd.extend(dwpi)
        ptk=set(range(dp.shape[0]))-set(dwpi)
        logging.debug(f'Deleting {dp.loc[dwpi].shape[0]} [ pairs ] due to dihedral deletions')
        # Note that we expect this to be zero if we are only deleting H's, since
        # an H can never be a 'j' or 'k' in a dihedral!
        self.D['pairs']=dp.take(list(ptk)).reset_index(drop=True)
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['dihedrals']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            d=self.D['dihedrals']
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
            d.ak=d.ak.map(mapper)
            d.al=d.al.map(mapper)
            d=self.D['pairs']
            d.ai=d.ai.map(mapper)
            d.aj=d.aj.map(mapper)
        if len(return_idx_of)>0:
            return new_idx
        self.null_check(msg='end of delete atoms')
        # self.to_file('tmp.top')
        return mapper
        
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
        """Merge topologies

        :param other: a topology
        :type other: Topology
        """
        # logging.debug('Topology.merge begins')
        self.merge_ex(other)
        self.merge_types(other)
        # logging.debug('Topology.merge ends')

    def dup_check(self,die=True):
        """Check for duplicate type-like topology records

        :param die: flag telling HTPolyNet to exit if duplicate found, defaults to True
        :type die: bool, optional
        :raises Exception: Exception raised if duplicate found and die is True
        """
        L=['atomtypes','bondtypes','angletypes']
        Not=' not' if not die else ''
        logging.debug(f'Checking for duplicate {L}; will{Not} die if found.')
        for t in L:
            i=_GromacsTopologyHashables_[t]
            ''' checking for types with duplicate atom-type indices '''
            dups=self.D[t].duplicated(subset=i,keep=False)
            if any(dups):
                logging.error(f'Duplicate {t} with different parameters detected\n'+self.D[t][dups].to_string())
                if die:
                    raise Exception('duplicate topology types with different parameters detected')

    def merge_types(self,other):
        """Merge type-like topology dataframes from other to self

        :param other: topology containing attribute D, a dictionary of dataframes
        :type other: Topology
        """
        L=['atomtypes','bondtypes','angletypes','dihedraltypes']
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

    def get_atom_attribute(self,idx,attribute):
        """Return value of attribute of atom idx

        :param idx: global atom index
        :type idx: int
        :param attribute: atom attribute name
        :type attribute: str
        :return: atom attribute value
        :rtype: varies
        """
        return self.D['atoms'].iloc[idx-1][attribute]

    def get_atomtype(self,idx):
        """Get atom type of atom with global index idx

        :param idx: atom global index
        :type idx: int
        :return: atom typ
        :rtype: str
        """
#        logging.debug(f'Asking get_atomtype for type of atom with index {idx}')
        return self.D['atoms'].iloc[idx-1].type

    def make_resid_graph(self,json_file=None,draw=None):
        adf=self.D['atoms']
        N=adf.shape[0]
        self.residue_network=nx.DiGraph()
        residues=adf['residue'].unique()
        for rn in residues:
            rs=adf[adf['residue']==rn]['resnr'].unique()
            self.residue_network.add_nodes_from(rs,resName=rn)
        
        resnrs=adf['resnr'].unique()
        for i in resnrs:
            # identify any inter-residue bonds emanating from
            # this residue, and classify
            #  - cross
            #  - polyethylene -- a bond from a carbon that is attached to another carbon that itself pariticpates in
            # an inter-residue crosslink

            # atom global indexes in this resnr
            ats=adf[adf['resnr']==i]['nr'].to_list()
            connectors=[]
            for a in ats:
                an=self.bondlist.partners_of(a)
                natsrn=adf[adf['nr'].isin(an)]['resnr'].unique()
                if len(natsrn)>1:
                    for n in natsrn:
                        if n!=i:  # this is an inter-residue connection
                            connectors.append((a,n))
            for c in connectors:
                a,n=c
                bondtype="cross"
                if not self.residue_network.has_edge(i,n):
                    self.residue_network.add_edge(i,n,bondtype=bondtype)
        if json_file:
            the_data=json_graph.node_link_data(self.residue_network)
            assert type(the_data)==dict,f'Error: node_link_data returns a {type(the_data)} but should return a dict'
            logging.debug(f'writing graph node_link_data to {json_file}')
            the_data_str=str(the_data)
            the_data_str=the_data_str.replace('True','"True"')
            the_data_str=the_data_str.replace('False','"False"')
            if "'" in the_data_str:
                logging.debug(f'json_graph.node_link_data produces single-quoted dict keys -- this is not JSON standard')
                json_compatible_string=the_data_str.replace("'",'"')
                try:
                    the_data=json.loads(json_compatible_string)
                except Exception as msg:
                    logging.debug(str(msg))
                    logging.debug(f'json.loads fails to encode string:\n{json_compatible_string}')
            with open (json_file,'w') as f:
                try:
                    json.dump(the_data,f)
                except Exception as msg:
                    logging.debug(str(msg))
                    logging.debug(f'writing resid graph to JSON not currently supported')
                    f.write(str(the_data)+'\n')
        if draw:
            network_graph(self.residue_network,draw)

    def copy_bond_parameters(self,bonds):
        """Generate and return a copy of a bonds dataframe that contains all bonds
           listed in bonds

        :param bonds: dataframe of bonds managed by runtime, 'ai','aj','reactantName'
        :type bonds: pandas.DataFrame
        :return: [ bonds ] dataframe extracted from system with all parameters
        :rtype: pandas.DataFrame
        """
        bdf=self.D['bonds']
        saveme=pd.DataFrame(columns=bdf.columns)
        for i,b in bonds.iterrows():
            ai,aj=idxorder((b['ai'],b['aj']))
            # logging.debug(f'copy parameters for ai {ai} aj {aj}')
            saveme=pd.concat((saveme,bdf[(bdf['ai']==ai)&(bdf['aj']==aj)].copy()),ignore_index=True)
        # logging.info(f'saved bond override params\n{saveme.to_string()}')
        return saveme

    def attenuate_bond_parameters(self,bondsdf,stage,max_stages,minimum_distance=0.0,init_colname='initial-distance'):
        """Alter the kb and b0 parameters for new crosslink bonds according to the values prior to 
            relaxation (stored in lengths), their equilibrium values, and the ratio stage/max_stages.
            Let stage/max_stages be x, and 1/max_stages <= x <= 1.  The spring constant for each
            bond is multiplied by x and the distance is 1 xth of the way from its maximum value 
            to its equilibrium value.

        :param bonds: dataframe of bonds managed by runtime, 'ai','aj','reactantName'
        :type bonds: pandas.DataFrame
        :param stage: index of stage in the series of post-bond-formation relaxation ("R" of SCUR)
        :type stage: int
        :param max_stages: total number of relaxation stages for this iteration
        :type max_stages: int
        :param minimum_distance: minimum bondlegth allowed, overriding type-specific b0 (if greater than 0)
        :type lengths: float
        """
        bdf=self.D['bonds']
        factor=(stage+1)/max_stages
        logging.debug(f'Attenuating {bondsdf.shape[0]} bond{"s" if bondsdf.shape[0]>1 else ""} in stage {stage+1}/{max_stages}')
        for i,b in bondsdf.iterrows():
            ai,aj=idxorder((b['ai'],b['aj']))
            rij=b[init_colname]
            b0,kb=self.get_bond_parameters(ai,aj)
            if minimum_distance>0.0:
                b0=minimum_distance
            new_b0=rij-factor*(rij-b0)
            new_kb=kb*factor
            # logging.debug(f'bond attenuation target for {ai}-{aj}:\nb0 {b0:.5f} kb {kb:.2f}; using b0 {new_b0:.5f} kb {new_kb:.2f}')
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c0']=new_b0
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c1']=new_kb

    def get_bond_parameters(self,ai,aj):
        """Gets b0 and kb for bond between atoms with global indexes ai and aj

        :param ai: global atom index
        :type ai: int
        :param aj: global atom index
        :type aj: int
        :return: b0, kb -- equilibrium bond length and spring constant
        :rtype: 2-tuple
        """
        ai,aj=idxorder((ai,aj))
        adf=self.D['atoms']
        bdf=self.D['bonds']
        tdf=self.D['bondtypes']
        b0=bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c0'].values[0]
        kb=bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c1'].values[0]
        if pd.isna(b0) or pd.isna(kb):
            ''' no overrides for this bond, so take from types '''
            it=adf.loc[adf['nr']==ai,'type'].values[0]
            jt=adf.loc[adf['nr']==aj,'type'].values[0]
            it,jt=typeorder((it,jt))
            b0=tdf.loc[(tdf['i']==it)&(tdf['j']==jt),'b0'].values[0]
            kb=tdf.loc[(tdf['i']==it)&(tdf['j']==jt),'kb'].values[0]
        return b0,kb

    def restore_bond_parameters(self,df):
        """Copy data from all bonds in dataframe df to global dataframe

        :param df: dataframe of bonds ['ai','aj','c0','c1']
        :type df: pandas DataFrame
        """
        bdf=self.D['bonds']
        for i,r in df.iterrows():
            ai,aj=r['ai'],r['aj']
            c0,c1=r['c0'],r['c1']
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c0']=c0
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c1']=c1

    def attenuate_pair_parameters(self,pairsdf,stage,max_stages,draglimit_nm=0.3):
        """Alter the kb and b0 parameters for new pre-crosslink pairs according 
            to the values prior to dragging (stored in pairdf['initial-distances']), 
            the desired lower limit of interatomic distance 'draglimit_nm', 
            and the ratio stage/max_stages.
            
        :param pairdf: pairs dataframe (['ai'],['aj'],['initial-distance'])
        :type pairdf: pandas.DataFrame
        :param stage: index of stage in the series of pre-bond-formation dragging
        :type stage: int
        :param max_stages: total number of drag stages for this iteration
        :type max_stages: int
        :param draglimit_nm: lower limit of interatomic distance requested from drag
        :type draglimit_nm: float
        """
        pdf=self.D['pairs']
        ess='s' if pairsdf.shape[0]>1 else ''
        factor=(stage+1)/max_stages
        logging.debug(f'Attenuating {pairsdf.shape[0]} pair{ess} in stage {stage+1}/{max_stages}')
        for i,b in pairsdf.iterrows():
            ai,aj=idxorder((b['ai'],b['aj']))
            b0=b['initial-distance']
            kb=pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c1'].values[0]
            pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c0']=draglimit_nm-factor*(b0-draglimit_nm)
            pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c1']=kb*factor
