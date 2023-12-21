"""

.. module:: topology
   :synopsis: Class for managing gromacs .top file data
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import json
import os
import pandas as pd
import numpy as np
import networkx as nx
from scipy.constants import physical_constants
from networkx.readwrite import json_graph
from itertools import product
from HTPolyNet.bondlist import Bondlist

logger=logging.getLogger(__name__)

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

def _get_unique_rings_dict(G,min_length=-1):
    """_get_unique_rings_dict generates dictionary identifying all unique covalent rings (chordless cycles) 
    from the graph G

    :param G: a graph representing atomic connectivity
    :type G: networkx.Graph
    :param min_length: minimum ring length, defaults to -1 (no limit; 3 would be a better choice)
    :type min_length: int, optional
    :return: dictionary of ring keyed on ring length with values being lists of lists of atom globalIdx's
    :rtype: dict
    """
    urings={}
    counts_by_length={}
    for u in nx.chordless_cycles(G):
        sl=len(u)
        if min_length<=sl:
            # logger.debug(f'a ring {u}')
            if not sl in counts_by_length:
                counts_by_length[sl]=0
            counts_by_length[sl]+=1
            if not sl in urings:
                urings[sl]=[]
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
                        if e in urings[l]:
                            found=True
                            break
                    if not found:
                        urings[l].append(u)
    return urings

def _present_and_contiguous(subL,L):
    """_present_and_contiguous returns True is elements in subL appear as a contiguous sub-block in 
       periodic and bidirectional list L

    :param subL: a sublist
    :type subL: list
    :param L: a list, treated as periodic and bidirectional
    :type L: list
    :return: True if sublist appears as part L
    :rtype: boolean
    """
    pretest=all([x in L for x in subL])
    if not pretest:
        return False
    for A in [L,L[::-1]]:
        T=treadmills(A)
        for t in T:
            testL=t[:len(subL)]
            # logger.debug(f'___ {subL} {testL}')
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

_GromacsTopologyDirectiveDefaults_={
    'system':['A_generic_system'],
    'molecules':['None',1],
    'moleculetype':['None',3],
    'defaults':[1,2,'yes',0.5,0.83333333]
}

def select_topology_type_option(options,typename='dihedraltypes',rule='stiffest'):
    """select_topology_type_option select from a list of topological interaction options of type typename using the provided rule

    :param options: list of parameterization options for a particular interaction
    :type options: list
    :param typename: string designation of interaction type, defaults to 'dihedraltypes'
    :type typename: str, optional
    :param rule: string describing the selection rule, defaults to 'stiffest'
    :type rule: str, optional
    :return: the selection parameterization option
    :rtype: element of options (dict)
    """
    hashables=_GromacsTopologyHashables_[typename]
    headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
    for i in hashables:
        headers.remove(i)
    if rule=='stiffest' or rule=='softest':
        if typename=='bondtypes':
            parindex=headers.index('kb')
        elif typename=='angletypes':
            parindex=headers.index('cth')
        elif typename=='dihedraltypes':
            parindex=headers.index('kd')
    sorted_options=list(sorted(options,key=lambda x: x[parindex]))
    if rule=='stiffest':
        return sorted_options[-1]
    elif rule=='softest':
        return sorted_options[0]
    return []

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
        self.D:dict[str:pd.DataFrame]={}
        for k,v in _GromacsTopologyDirectiveDefaults_.items():
            hdr=_GromacsTopologyDirectiveHeaders_[k]
            dfdict={kk:[a] for kk,a in zip(hdr,v)}
            self.D[k]=pd.DataFrame(dfdict)
        self.D['system']=pd.DataFrame({'name':[system_name]})
        ''' bondlist: a class that owns a dictionary keyed on atom global index with values that are lists of global atom indices bound to the key '''
        self.bondlist=Bondlist()
        self.residue_network=nx.DiGraph()
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
        dirname=os.path.dirname(filename)
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
                        # logger.info(f'Found second set of {len(tdf)} [ dihedraltypes ] in {inst.filename}; merging into set of {len(inst.D["dihedraltypes"])} types already read in...')
                        # we have already read-in a dihedraltypes section
                        # so let's append this one
                        inst.D['dihedraltypes']=pd.concat([inst.D['dihedraltypes'],tdf],ignore_index=True).drop_duplicates()
                        # logger.info(f'    -> now there are {len(inst.D["dihedraltypes"])} dihedral types.')
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
                inst.D['atomtypes']=inst.D['atomtypes'].sort_values(by='name').reset_index(drop=True)
#                inst.D['atomtypes'].sort_values(by='name',inplace=True)
            if 'bonds' in inst.D:
                inst.D['bonds']=inst.D['bonds'].sort_values(by=['ai','aj']).reset_index(drop=True)
                inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
            if 'bondtypes' in inst.D:
                df_typeorder(inst.D['bondtypes'],typs=['i','j'])
#                inst.D['bondtypes'].sort_values(by=['i','j'],inplace=True)
                inst.D['bondtypes']=inst.D['bondtypes'].sort_values(by=['i','j']).reset_index(drop=True)
            if 'pairs' in inst.D:
                inst.D['pairs']=inst.D['pairs'].sort_values(by=['ai','aj']).reset_index(drop=True)
            if 'pairtypes' in inst.D:
                df_typeorder(inst.D['pairtypes'],typs=['i','j'])
#                inst.D['pairtypes'].sort_values(by=['i','j'],inplace=True)
                inst.D['pairtypes']=inst.D['pairtypes'].sort_values(by=['i','j']).reset_index(drop=True)
            if 'angles' in inst.D:
                # central atom (aj) is the primary index
                inst.D['angles']=inst.D['angles'].sort_values(by=['aj','ai','ak']).reset_index(drop=True)
            if 'angletypes' in inst.D:
                df_typeorder(inst.D['angletypes'],typs=['i','j','k'])
#                inst.D['angletypes'].sort_values(by=['j','i','k'],inplace=True)
                inst.D['angletypes']=inst.D['angletypes'].sort_values(by=['j','i','k']).reset_index(drop=True)
            if 'dihedrals' in inst.D:
                # central atoms (aj,ak) are the primary index
                inst.D['dihedrals']=inst.D['dihedrals'].sort_values(by=['aj','ak','ai','al']).reset_index(drop=True)
            if 'dihedraltypes' in inst.D:
                df_typeorder(inst.D['dihedraltypes'],typs=['i','j','k','l'])
                # print(f'    -> pre sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
#                inst.D['dihedraltypes'].sort_values(by=['j','k','i','l'],inplace=True)
                inst.D['dihedraltypes']=inst.D['dihedraltypes'].sort_values(by=['j','k','i','l']).reset_index(drop=True)
                # print(f'    -> post sort: now there are {len(inst.D["dihedraltypes"])} dihedral types.')
            for f in inst.includes:
                # print(f'reading included topology {f}')
                inst.merge(Topology.read_gro(os.path.join(dirname,f)))
            inst.empty=False
            return inst

    def bond_source_check(self):
        """bond_source_check Checks to ensure the 'bonds' dataframe and 'mol2_bonds' dataframe contain the same bonds.  A mol2 dataframe is only created when a mol2 file is read by the Coordinates module.
        """
        if 'bonds' in self.D and 'mol2_bonds' in self.D:
            # logger.debug(f'Consistency check between gromacs-top bonds and mol2-bonds requested.')
            grobonds=self.D['bonds'].sort_values(by=['ai','aj'])
            bmi=grobonds.set_index(['ai','aj']).index
            mol2bonds=self.D['mol2_bonds'].sort_values(by=['ai','aj'])
            mbmi=mol2bonds.set_index(['ai','aj']).index
            check=all([x==y for x,y in zip(bmi,mbmi)])
            # logger.debug(f'Result: {check}')
            if not check:
                logger.error(f'Gromacs/Mol2 bond inconsistency detected')
                logger.debug(f'GROMACS:')
                logger.debug(grobonds[['ai','aj']].head().to_string())
                logger.debug(f'MOL2:')
                logger.debug(mol2bonds[['ai','aj']].head().to_string())
                for x,y in zip(bmi,mbmi):
                    logger.debug(f'{x} {y} {x==y}')

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

    def detect_rings(self):
        """detect_rings detect unique rings in the topology

        :return: ring dict (length:list-of-rings-by-atom-indices)
        :rtype: dict
        """
        g=self.bondlist.graph()
        cycles=_get_unique_rings_dict(g,min_length=3)
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
                            logger.debug(f'{msg} null in {k} {a}\n{self.D[k].to_string()}')
                            raise Exception('NaN error')

    def total_charge(self):
        """Compute and return total system charge

        :return: charge
        :rtype: float
        """
        if 'atoms' in self.D:
            return self.D['atoms']['charge'].sum()
        return 0.0

    def adjust_charges(self,atoms=[],desired_charge=0.0,overcharge_threshhold=0.1,msg=''):
        """Adjust atom partial charges a tiny bit so that total system charge is zero

        :param desired_charge: target system charge, defaults to 0.0
        :type desired_charge: float, optional
        :param overcharge_threshhold: threshold overcharge that triggers a message, defaults to 0.1
        :type overcharge_threshhold: float, optional
        :param msg: A message to write if pre-adjusted system charge is too high, defaults to ''
        :type msg: str, optional
        :return: self topology
        :rtype: Topology
        """
        apparent_charge=self.total_charge()
        overcharge=apparent_charge-desired_charge
        logger.debug(f'Adjusting charges of {len(atoms)} atoms due to overcharge of {overcharge:.6f}')
        if np.abs(overcharge)>overcharge_threshhold:
            logger.debug(f'{msg}')
        cpa=-overcharge/len(atoms)
        logger.debug(f'Adjustment is {cpa:.4e} per atom')
        for i in atoms:
            self.D['atoms'].loc[i-1,'charge']+=cpa
        logger.debug(f'New total charge after adjustment: {self.total_charge():.6f}')
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
            logger.debug(f'mass {M_amu} fac {fac}')
            return M_amu*fac
        return 0.0

    def atomcount(self):
        """atomcount Returns the total number of atoms

        :return: number of atoms
        :rtype: int
        """
        if 'atoms' in self.D:
            return self.D['atoms'].shape[0]
        return 0
    
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
            b0=b['initial_distance']
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

    def add_bonds(self,pairs=[]):
        """add_bonds Adds bonds indicated in list pairs to the topology

        :param pairs: list of pairs of atom indexes, defaults to []
        :type pairs: list, optional
        :raises Exception: dies if an existing bond is in the list of pairs
        """
        # logger.debug('begins')
        at=self.D['atoms']
        ij=self.D['bondtypes'].set_index(['i','j'])
        #mb=self.D['mol2_bonds']
        bmi=self.D['bonds'].set_index(['ai','aj']).sort_index().index
        pmi=self.D['pairs'].set_index(['ai','aj']).sort_index().index
        newbonds=[]
        for b in pairs:
            # logger.debug(f'{b}')
            # assert type(b[0])==int
            bondtuple=(int(b[0]),int(b[1]))
            order=int(b[2])
            ai,aj=idxorder(bondtuple)
            assert type(ai)==int
            assert type(aj)==int
            '''
            if this bond is not in the topology, then add it
            '''
            if not (ai,aj) in bmi:
                newbonds.append((ai,aj))
                # logger.debug(f'asking types of {ai} and {aj}; at.shape {at.shape}')
                it=at.iloc[ai-1].type
                jt=at.iloc[aj-1].type
                idx=typeorder((it,jt))
                if idx in ij.index:
                    bt=ij.loc[idx,'func']  # why don't i need need values[0]
                    kb=ij.loc[idx,'kb']
                    b0=ij.loc[idx,'b0']
                else:
                    logger.debug(f'no bondtype {idx} found; using placeholder parameters')
                    bt=1
                    b0=0.15
                    kb=999999
                    # raise Exception(f'no bondtype {idx} found.')
                '''
                add a new bond!
                '''
                h=_GromacsTopologyDirectiveHeaders_['bonds']
                data=[ai,aj,bt,b0,kb]  # this new bond will have override parameters
                assert len(h)==len(data), 'Error: not enough data for new bond?'
                bonddict={k:[v] for k,v in zip(h,data)}
                bdtoadd=pd.DataFrame(bonddict)
                self.D['bonds']=pd.concat((self.D['bonds'],bdtoadd),ignore_index=True)
                # logger.info(f'add_bond:\n{bdtoadd.to_string()}')
                logger.debug(f'just added {bonddict}')
                if 'mol2_bonds' in self.D:
                    data=[len(self.D['mol2_bonds']),ai,aj,1] # assume single bond
                    bonddict={k:[v] for k,v in zip(['bondIdx','ai','aj','order'],data)}
                    self.D['mol2_bonds']=pd.concat((self.D['mol2_bonds'],pd.DataFrame(bonddict)),ignore_index=True)
                # remove this pair from pairs if it's in there (it won't be)
                if idx in pmi:
                    logger.debug(f'Warning: new bond {ai}-{aj} was evidently in the [ pairs ]!')
                    d=self.D['pairs']
                    indexes_to_drop=d[(d.ai.isin(idx))&(d.aj.isin(idx))].index
                    logger.debug(f'Dropping [ pair ]:\n{d[indexes_to_drop].to_string()}')
                    indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                    self.D['pairs']=d.take(list(indexes_to_keep)).reset_index(drop=True)
            else:
                ''' 
                if it is, do nothing; it will be templated; if mol2_bonds are present (usually
                because a Topology is part of a molecule being parameterized), update the order
                of the bond.
                '''
                if 'mol2_bonds' in self.D:
                    mb=self.D['mol2_bonds']
                    bi=(mb['ai']==ai)&(mb['aj']==aj)
                    mb.loc[bi,'order']=order
        '''
        update the bondlist
        '''
        for b in newbonds:
            self.bondlist.append(b)
        logger.debug(f'Added {len(newbonds)} new bonds')


    def delete_atoms(self,idx=[],reindex=True,return_idx_of=[],**kwargs):
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
        #logger.debug(f'Delete atoms: {idx}')
        paranoid_about_pairs=kwargs.get('paranoid_about_pairs',False)
        self.null_check(msg='beginning of delete atoms')
        # logger.debug(f'idx {idx}')
        d=self.D['atoms']
        new_idx=[]
        indexes_to_drop=d[d.nr.isin(idx)].index
        total_missing_charge=d.loc[indexes_to_drop]['charge'].sum()
        logger.debug(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ atoms ]; charge to make up: {total_missing_charge:.4f}')#:\n{d.loc[indexes_to_drop].to_string()}')
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        self.D['atoms']=d.take(list(indexes_to_keep)).reset_index(drop=True)
        mapper={}
        if reindex:
            d=self.D['atoms']
            oldGI=d['nr'].copy()
            d['nr']=d.index+1
            mapper={k:v for k,v in zip(oldGI,d['nr'])}
            # logger.debug(f'mapper {mapper}')
            assert not any([x in mapper for x in idx]),f'Error: Some deleted atoms in mapper.'
            k=np.array(list(mapper.keys()))
            v=np.array(list(mapper.values()))
            if any(np.isnan(k)):
                logger.error('null in mapper keys')
            if any(np.isnan(v)):
                logger.error('null in mapper values')
            # logger.debug(f'delete_atoms: mapper {mapper}')
            if len(return_idx_of)>0:
                # logger.info(f'Asking for updated global indexes of {return_idx_of}')
                new_idx=[mapper[o] for o in return_idx_of]
            #d['nr_shift']=d['nr']-oldGI  # probably not necessary
        ptt=['bonds','mol2_bonds','pairs']
        for pt in ptt:
            if pt in self.D:
                d=self.D[pt]
                # logger.debug(f'delete atom: {pt} df prior to deleting')
                # logger.debug(d.to_string())
                indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
                logger.debug(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ {pt} ]')
                indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
                self.D[pt]=d.take(list(indexes_to_keep)).reset_index(drop=True)
                if reindex:
                    d=self.D[pt]
                    assert not any([x in idx for x in d.ai]),f'Error: deleted atom survived in {pt} ai'
                    assert not any([x in idx for x in d.aj]),f'Error: deleted atom survived in {pt} aj'
                    assert all([x in mapper for x in d.ai]),f'Error: surviving {pt} atom ai old idx not in mapper'
                    assert all([x in mapper for x in d.aj]),f'Error: surviving {pt} atom aj old idx not in mapper'
                    #logger.debug(f'delete atom: {pt} df prior to reindexing')
                    #logger.debug(d.to_string())
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
                        # logger.debug(f'Updating bondlist using (first 10 shown)\n{d.head(10).to_string()}')
                        self.bondlist=Bondlist.fromDataFrame(d)
                    if pt=='mol2_bonds':
                        nBonds=self.D[pt].shape[0]
                        self.D[pt]['bondIdx']=list(range(1,nBonds+1))
        d=self.D['angles']
        # assert d.ai.dtype==int,f'pre-delete lost angle ai dtype {d.ai.dtype}'
        # assert d.aj.dtype==int,f'pre-delete lost angle aj dtype {d.aj.dtype}'
        # assert d.ak.dtype==int,f'pre-delete lost angle ak dtype {d.ak.dtype}'
        # logger.debug(f'ai {d.ai.isin(idx).to_string()}')
        # logger.debug(f'aj {d.aj.isin(idx).to_string()}')
        # logger.debug(f'ak {d.ak.isin(idx).to_string()}')
        indexes_to_drop=d[(d['ai'].isin(idx))|(d['aj'].isin(idx))|(d['ak'].isin(idx))].index
        # extras=d[d['ak'].isin(idx)].index
        logger.debug(f'Deleting {len(indexes_to_drop)} [ angles ]')
        indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
        # logger.debug(f'drop {list(sorted(list(set(indexes_to_drop))))}')
        # logger.debug(f'keep {indexes_to_keep}')
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
                logger.debug(f'Zombie ak angles:\n{self.D["angles"][zombie_tags].to_string()}')
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
        logger.debug(f'Deleting {d.loc[indexes_to_drop].shape[0]} [ dihedrals ]')
        # if the atoms we have deleted are truly just H's, then there will be no other
        # spurious pairs after all dihedrals containing deleted atoms are deleted.
        # However, we may want to still search for such pairs, so let's leave this
        # as an option:
        if paranoid_about_pairs:
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
            if dp.loc[dwpi].shape[0]>0:
                logger.debug(f'  -> and deleting {dp.loc[dwpi].shape[0]} [ pairs ] from those dihedrals')
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
            tp=[]
            pdrops=[]
            for i,r in d.iterrows():
                ai,aj=min([r.ai,r.aj]),max([r.ai,r.aj])
                if not [ai,aj] in tp:
                    tp.append([ai,aj])
                else:
                    pdrops.append(i)
            logger.debug(f'Deleting {len(pdrops)} duplicate 1-4 pair descriptors -- this is likely due to a bug somewhere')
            self.D['pairs']=d.drop(pdrops).reset_index(drop=True)
        self.null_check(msg='end of delete atoms')
        logger.debug('finished.')
        if len(return_idx_of)>0:
            return new_idx
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
        # logger.debug('Topology.merge begins')
        # look for duplicated types between self and other.  If any are found, delete those types from other and copy their parameters into the explicit interactions they correspond to.
        self.merge_ex(other)
        self.merge_types(other)
        # logger.debug('Topology.merge ends')

    # def handle_duplicate_types(self,other,copy_directive='other_to_self',typename='',funcidx=4):
#         if not typename in self.D or not typename in other.D:
#             return
#         stdf=self.D[typename]
#         otdf=other.D[typename]
#         hashables=_GromacsTopologyHashables_[typename]
#         headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
#         logger.debug(f'{hashables} {headers}')
#         for i in hashables:
#             headers.remove(i)
#         common=[]
#         fi=headers.index('func')
#         copy_idx_pairs=[]
#         for idx,r in otdf.iterrows():
#             typidx=typeorder(tuple([r[i] for i in hashables]))
#             idata=[r[i] for i in headers]
#             for jdx,q in stdf.iterrows():
#                 typjdx=typeorder(tuple([q[i] for i in hashables]))
#                 jdata=[q[i] for i in headers]
#                 if typidx==typjdx and idata==jdata:
#                     common.append(typidx)
#             for jdx,q in stdf.iterrows():
#                 typjdx=typeorder(tuple([q[i] for i in hashables]))
#                 jdata=[q[i] for i in headers]
#                 if typidx==typjdx and idata!=jdata and not typidx in common and idata[fi]==funcidx and jdata[fi]==funcidx:
#                     logger.debug(f'duplicate {typename} {typidx}')
#                     logger.debug(f'o {idx} {idata}')
#                     logger.debug(f's {jdx} {jdata}')
#                     if not (idx,jdx) in copy_idx_pairs:
#                         copy_idx_pairs.append((idx,jdx))
#         if copy_directive=='other_to_self' and len(copy_idx_pairs)>0:
#             for o,s in copy_idx_pairs:
#                 self.D[typename].iloc[s]=other.D[typename].iloc[o]
#         elif copy_directive=='self_to_other' and len(copy_idx_pairs)>0:
#             for o,s in copy_idx_pairs:
#                 other.D[typename].iloc[o]=self.D[typename].iloc[s]
#           # logger.debug(f'{typidx}')
# #        pass
    def report_type(self,typidx_q,typename='',funcidx=-1):
        if not typename in self.D:
            return
        stdf=self.D[typename]
        hashables=_GromacsTopologyHashables_[typename]
        headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
        for i in hashables:
            headers.remove(i)
        fi=headers.index('func')
        for idx,r in stdf.iterrows():
            typidx=typeorder(tuple([r[i] for i in hashables]))
            idata=[r[i] for i in headers]
            if typidx==typidx_q and idata[fi]==funcidx:
                return [r[i] for i in headers]
        return []
    '''
    'bonds':        ['ai', 'aj', 'funct', 'c0', 'c1'],
    'bondtypes':    [ 'i',  'j', 'func',  'b0', 'kb'],
    'angles':       ['ai', 'aj', 'ak', 'funct', 'c0',  'c1'],
    'angletypes':   [ 'i',  'j',  'k', 'func',  'th0', 'cth', 'rub', 'kub'],
    'dihedrals':    ['ai', 'aj', 'ak', 'al', 'funct', 'c0',    'c1', 'c2', 'c3', 'c4', 'c5'],
    'dihedraltypes':[ 'i',  'j',  'k',  'l', 'func',  'phase', 'kd', 'pn'],
    '''
    def reset_override_from_type(self,interactionname,typename,inst_idx):
        if not typename in self.D:
            return
        typ_hashables=_GromacsTopologyHashables_[typename]
        stdf=self.D[typename].set_index(typ_hashables)
        sidf=self.D[interactionname]
        typ_headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
        ins_hashables=_GromacsTopologyHashables_[interactionname]
        ins_headers=_GromacsTopologyDirectiveHeaders_[interactionname].copy()
        typidx=[self.D['atoms'].iloc[x-1]['type'] for x in inst_idx]
        typidx=typeorder(tuple(typidx))
        for i in typ_hashables:
            typ_headers.remove(i)
        for i in ins_hashables:
            ins_headers.remove(i)
        iidx=idxorder(tuple(inst_idx))
        idx=-1
        for i,r in sidf.iterrows():
            jdx=tuple([r[a] for a in ins_hashables])
            if iidx==jdx:
                idx=i
                break
        assert idx!=-1
        num_data=min([len(typ_headers),len(ins_headers)])
        typ_headers=typ_headers[:num_data]
        ins_headers=ins_headers[:num_data]
        typrec=stdf.loc[typidx][typ_headers]
        cols=self.D[interactionname].columns.get_indexer(ins_headers)
        logger.debug(f'Resetting override in {interactionname} for {inst_idx} from')
        for ln in self.D[interactionname].iloc[idx,cols].to_string().split('\n'):
            logger.debug(ln)
        logger.debug('to')
        self.D[interactionname].iloc[idx,cols]=typrec
        for ln in self.D[interactionname].iloc[idx,cols].to_string().split('\n'):
            logger.debug(ln)
        

    def reset_type(self,typename,typidx_t,values):
        if not typename in self.D:
            return
        stdf=self.D[typename]
        hashables=_GromacsTopologyHashables_[typename]
        headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
        for i in hashables:
            headers.remove(i)
        cols=self.D[typename].columns.get_indexer(headers)
        idxs=[]
        for idx,r in stdf.iterrows():
            typidx=typeorder(tuple([r[i] for i in hashables]))
            if typidx==typidx_t:
                idxs.append(idx)
        assert len(headers)==len(values)
        logger.debug(f'Resetting {len(idxs)} entries {headers} to {values}')
        for idx in idxs:
            self.D[typename].iloc[idx,cols]=values
            for ln in self.D[typename].iloc[idx][headers].to_string().split('\n'):
                logger.debug(ln)

    def report_duplicate_types(self,other,typename='',funcidx=4):
        if not typename in self.D or not typename in other.D:
            return
        stdf=self.D[typename]
        otdf=other.D[typename]
        hashables=_GromacsTopologyHashables_[typename]
        headers=_GromacsTopologyDirectiveHeaders_[typename].copy()
        # logger.debug(f'{hashables} {headers}')
        for i in hashables:
            headers.remove(i)
        common=[]
        fi=headers.index('func')
        true_duplicate_types=[]
        for idx,r in otdf.iterrows():
            typidx=typeorder(tuple([r[i] for i in hashables]))
            idata=[r[i] for i in headers]
            for jdx,q in stdf.iterrows():
                typjdx=typeorder(tuple([q[i] for i in hashables]))
                jdata=[q[i] for i in headers]
                if typidx==typjdx and idata==jdata:
                    common.append(typidx)
        for idx,r in otdf.iterrows():
            typidx=typeorder(tuple([r[i] for i in hashables]))
            idata=[r[i] for i in headers]
            for jdx,q in stdf.iterrows():
                typjdx=typeorder(tuple([q[i] for i in hashables]))
                jdata=[q[i] for i in headers]
                if typidx==typjdx and idata!=jdata and not typidx in common and idata[fi]==funcidx and jdata[fi]==funcidx:
                    # logger.debug(f'duplicate {typename} {typidx}')
                    true_duplicate_types.append(typidx)
        return true_duplicate_types

    def dup_check(self,die=True):
        """Check for duplicate type-like topology records

        :param die: flag telling HTPolyNet to exit if duplicate found, defaults to True
        :type die: bool, optional
        :raises Exception: Exception raised if duplicate found and die is True
        """
        L=['atomtypes','bondtypes','angletypes']
        Not=' not' if not die else ''
        logger.debug(f'Checking for duplicate {L}; will{Not} die if found.')
        for t in L:
            i=_GromacsTopologyHashables_[t]
            ''' checking for types with duplicate atom-type indices '''
            dups=self.D[t].duplicated(subset=i,keep=False)
            if any(dups):
                logger.error(f'Duplicate {t} with different parameters detected\n'+self.D[t][dups].to_string())
                if die:
                    raise Exception('duplicate topology types with different parameters detected')

    def merge_types(self,other):
        """Merge type-like topology dataframes from other to self

        :param other: topology containing attribute D, a dictionary of dataframes
        :type other: Topology
        """
        # self.handle_duplicate_types(other,typename='dihedraltypes',funcidx=4,drop_directive='drop_from_self')
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
#        logger.debug(f'Asking get_atomtype for type of atom with index {idx}')
        return self.D['atoms'].iloc[idx-1].type

    def build_interresidue_graph(self,G,ri):
        if ri in G:
            return
        adf=self.D['atoms']
        ri_at_idx=adf[adf['resnr']==ri]['nr'].to_list()
        residues=adf['resnr'].unique()
        ri_partners=[]
        for rn in residues:
            skip=False
            rj_at_idx=adf[adf['resnr']==rn]['nr'].to_list()
            for i,j in product(ri_at_idx,rj_at_idx):
                # logger.debug(f'### {i} in {self.bondlist.partners_of(j)}?')
                if i in self.bondlist.partners_of(j):
                    if not rn in ri_partners:
                        ri_partners.append(rn)
                        skip=True
                        continue
                if skip: continue
        # logger.debug(f'rn partners {ri_partners}')
        G.add_node(ri)
        for rj in ri_partners:
            self.build_interresidue_graph(G,rj)

    def local_resid_cluster(self,ri):
        G=nx.DiGraph()
        # logger.debug(str(self.bondlist))
        self.build_interresidue_graph(G,ri)
        return [x for x in G]

    def make_resid_graph(self,json_file=None):
        adf=self.D['atoms']
        N=adf.shape[0]
        self.residue_network=nx.DiGraph()
        residues=adf['residue'].unique()
        for rn in residues:
            rs=adf[adf['residue']==rn]['resnr'].unique()
            self.residue_network.add_nodes_from(rs,resName=rn)
        
        resnrs=adf['resnr'].unique()
        for i in resnrs:
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
            logger.debug(f'writing graph node_link_data to {json_file}')
            the_data_str=str(the_data)
            the_data_str=the_data_str.replace('True','"True"')
            the_data_str=the_data_str.replace('False','"False"')
            if "'" in the_data_str:
                logger.debug(f'json_graph.node_link_data produces single-quoted dict keys -- this is not JSON standard')
                json_compatible_string=the_data_str.replace("'",'"')
                try:
                    the_data=json.loads(json_compatible_string)
                except Exception as msg:
                    logger.debug(str(msg))
                    logger.debug(f'json.loads fails to encode string:\n{json_compatible_string}')
            with open (json_file,'w') as f:
                try:
                    json.dump(the_data,f)
                except Exception as msg:
                    logger.debug(str(msg))
                    logger.debug(f'writing resid graph to JSON not currently supported')
                    f.write(str(the_data)+'\n')
        # if draw:
        #     network_graph(self.residue_network,draw)

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
        for b in bonds.itertuples():
            ai,aj=idxorder((b.ai,b.aj))
            # logger.debug(f'copy parameters for ai {ai} aj {aj}')
            saveme=pd.concat((saveme,bdf[(bdf['ai']==ai)&(bdf['aj']==aj)].copy()),ignore_index=True)
        # logger.info(f'saved bond override params\n{saveme.to_string()}')
        return saveme

    def attenuate_bond_parameters(self,bondsdf,stage,max_stages,minimum_distance=0.0,init_colname='initial_distance'):
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
        logger.debug(f'Attenuating {bondsdf.shape[0]} bond{"s" if bondsdf.shape[0]>1 else ""} in stage {stage+1}/{max_stages}')
        jdx=list(bondsdf.columns).index(init_colname)
        for b in bondsdf.itertuples(index=False):
            ai,aj=idxorder((b.ai,b.aj))
            rij=b[jdx]
            b0,kb=self.get_bond_parameters(ai,aj)
            if minimum_distance>0.0:
                b0=minimum_distance
            new_b0=rij-factor*(rij-b0)
            new_kb=kb*factor
            # logger.debug(f'bond attenuation target for {ai}-{aj}:\nb0 {b0:.5f} kb {kb:.2f}; using b0 {new_b0:.5f} kb {new_kb:.2f}')
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
        for r in df.itertuples():
            ai,aj=r.ai,r.aj
            c0,c1=r.c0,r.c1
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c0']=c0
            bdf.loc[(bdf['ai']==ai)&(bdf['aj']==aj),'c1']=c1

    def attenuate_pair_parameters(self,pairsdf,stage,max_stages,draglimit_nm=0.3):
        """Alter the kb and b0 parameters for new pre-crosslink pairs according 
            to the values prior to dragging (stored in pairdf['initial_distances']), 
            the desired lower limit of interatomic distance 'draglimit_nm', 
            and the ratio stage/max_stages.
            
        :param pairdf: pairs dataframe (['ai'],['aj'],['initial_distance'])
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
        logger.debug(f'Attenuating {pairsdf.shape[0]} pair{ess} in stage {stage+1}/{max_stages}')
        for b in pairsdf.itertuples():
            ai,aj=idxorder((b.ai,b.aj))
            b0=b.initial_distance
            kb=pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c1'].values[0]
            pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c0']=draglimit_nm-factor*(b0-draglimit_nm)
            pdf.loc[(pdf['ai']==ai)&(pdf['aj']==aj),'c1']=kb*factor
