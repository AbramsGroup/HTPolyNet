'''
coordinates.py -- simple class for handling atom coordinates from gro and mol2 files
'''

import pandas as pd
import numpy as np
from io import StringIO
from copy import deepcopy
import hashlib
import logging

from HTPolyNet.bondlist import Bondlist

_ANGSTROM_='Ångström'

class Coordinates:
    gro_colnames = ['resNum', 'resName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    gro_colunits = ['*','*','*','*','nm','nm','nm','nm/ps','nm/ps','nm/ps']
    mol2_atom_colnames = ['globalIdx','atomName','posX','posY','posZ','type','resNum','resName','charge']
    mol2_atom_colunits = ['*','*',_ANGSTROM_,_ANGSTROM_,_ANGSTROM_,'*','*']
    mol2_bond_colnames = ['globalIdx','ai','aj','type']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_colnames, [int, int, int, str])}
    coord_aux_attributes = ['sea-idx','reactive-role']

    def __init__(self,name=''):
        self.name=name
        self.format='gro'
        self.units={}
        self.units['length']='nm'
        self.units['velocity']='nm/ps'
        self.title=''
        self.metadat={}
        self.N=0
        ''' dataframes: 'atoms' and 'bonds' '''
        self.D={}
        self.box=np.zeros((3,3))
        
    @classmethod
    def read_gro(cls,filename=''):
        inst=cls(filename)
        if filename!='':
            with open(filename,'r') as f:
                data=f.read().split('\n')
                while '' in data:
                    data.remove('')
                inst.title=data[0]
                inst.N=int(data[1])
                series={k:[] for k in cls.gro_colnames}
                for x in data[2:-1]:
                    series['resNum'].append(int(x[0:5].strip()))
                    series['resName'].append(x[5:10].strip())
                    series['atomName'].append(x[10:15].strip())
                    ''' if formatted correctly, globalIdx is row index + 1 always! '''
                    series['globalIdx'].append(int(x[15:20].strip()))
                    numbers=list(map(float,[y.strip() for y in x[20:].split()]))
                    series['posX'].append(numbers[0])
                    series['posY'].append(numbers[1])
                    series['posZ'].append(numbers[2])
                    if len(numbers)==6:
                        series['velX'].append(numbers[3])
                        series['velY'].append(numbers[4])
                        series['velZ'].append(numbers[5])
                if len(series['velX'])==0:
                    del series['velX']
                    del series['velY']
                    del series['velZ']
                assert inst.N==len(series['globalIdx']), f'Atom count mismatch inside {filename}'
                # for k,v in series.items():
                #     logging.debug(f'in coordinates.read_gro: {k} has {len(v)} items.')
                inst.D['atoms']=pd.DataFrame(series)
                boxdataline=data[-1]
                n=10
                boxdata=list(map(float,[boxdataline[i:i+n].strip() for i in range(0,len(boxdataline),n)]))
                inst.box[0][0],inst.box[1][1],inst.box[2][2]=boxdata[0:3]
                if len(boxdata)==9:
                    inst.box[0][1],inst.box[0][2],inst.box[1][0],inst.box[1][2],inst.box[2][0],inst.box[2][1]=boxdata[3:]
        return inst

    @classmethod
    def read_mol2(cls,filename=''):
        ''' Reads in a Sybyl MOL2 file into a Coordinates instance. 
            Note that this method only reads in
            MOLECULE, ATOM, and BOND sections.  '''
        inst=cls(filename)
        inst.format='mol2'
        ''' Length units in MOL2 are always Ångström '''
        inst.units['length']=_ANGSTROM_
        inst.units['velocity']=_ANGSTROM_+'/ps'
        if filename=='':
            return inst
        with open(filename,'r') as f:
            rawsections=f.read().split('@<TRIPOS>')[1:]
            sections={}
            for rs in rawsections:
                s=rs.split('\n')
                key=s[0].strip().lower()
                val=[a.strip() for a in s[1:] if len(a)>0]
                if key=='atom' or key=='bond':
                    val=StringIO('\n'.join(val))
                sections[key]=val
            inst.name=sections['molecule'][0]
            imetadat=list(map(int,sections['molecule'][1].strip().split()))
            inst.metadat['N']=inst.N=imetadat[0]
            inst.metadat['nBonds']=imetadat[1]
            inst.metadat['nSubs']=imetadat[2]
            inst.metadat['nFeatures']=imetadat[3]
            inst.metadat['nSets']=imetadat[4]
            inst.metadat['mol2type']=sections['molecule'][2]
            inst.metadat['mol2chargetype']=sections['molecule'][3]
            inst.D['atoms']=pd.read_csv(sections['atom'],sep='\s+',names=Coordinates.mol2_atom_colnames)
            inst.N=len(inst.D['atoms'])
            inst.D['bonds']=pd.read_csv(sections['bond'],sep='\s+',names=Coordinates.mol2_bond_colnames,dtype=Coordinates.mol2_bond_types)
            inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
        return inst

    def copy_coords(self,other):
        assert len(self.D['atoms'])==len(other.D['atoms']),f'Cannot copy -- atom count mismatch {len(self.D["atoms"])} vs {len(other.D["atoms"])}'
        otherfac=1.0
        if self.units['length']=='nm' and other.units['length']==_ANGSTROM_:
            otherfac=0.1
        elif self.units['length']==_ANGSTROM_ and other.units['length']=='nm':
            otherfac=10.0
        for c in ['posX','posY','posZ']:
            otherpos=other.D['atoms'][c].copy()
            otherpos*=otherfac
            self.D['atoms'][c]=otherpos
            
    def translate_coords(self,v=[0.,0.,0.]):
        self.D['atoms']['posX']+=v[0]
        self.D['atoms']['posY']+=v[1]
        self.D['atoms']['posZ']+=v[2]

    def geometric_center(self):
        a=self.D['atoms']
        return np.array([a.posX.mean(),a.posY.mean(),a.posZ.mean()])

    def rij(self,i,j,pbc=[1,1,1]):
        ''' compute distance between atoms i and j
            We assume that the DF row index is the
            globalIdx-1! '''
        ri=self.D['atoms'].iloc[i-1][['posX','posY','posZ']].values
        rj=self.D['atoms'].iloc[j-1][['posX','posY','posZ']].values
        rij=ri-rj
        for c in range(3):
            if pbc[c]:
                hbx=self.box[c][c]/2
                if rij[c]<-hbx:
                    rij[c]+=self.box[c][c]
                elif rij[c]>hbx:
                    rij[c]-=self.box[c][c]
        logging.info(f'i {i} j {j} rij {rij} mag {np.sqrt(rij.dot(rij))}')
        logging.info(f'pbc {pbc} box {self.box[0][0]} {self.box[1][1]} {self.box[2][2]}')
        return np.sqrt(rij.dot(rij))

    def calc_distance_matrix(self):
        M=np.zeros((self.N,self.N))
        for i in range(self.N-1):
            for j in range(i+1,self.N):
                d=self.rij(i,j)
                M[i][j]=M[j][i]=d
        self.distance_matrix=M

    def merge(self,other):
        if self.units['length']!=other.units['length']:
            raise Exception('Cannot merge coordinates in different units')
        ''' get atom index, bond index, and resnum index shifts '''
        idxshift=0 if 'atoms' not in self.D else len(self.D['atoms'])
        bdxshift=0 if 'bonds' not in self.D else len(self.D['bonds'])
        rdxshift=0 if 'atoms' not in self.D else self.D['atoms'].iloc[-1]['resNum']
        if 'atoms' in other.D:
            ''' shift residue indices in other before merging '''
            other.D['atoms']['resNum']+=rdxshift
        nOtherBonds=0
        if 'bonds' in other.D:
            ''' count number of bonds in other '''
            nOtherBonds=len(other.D['bonds'])
            ''' shift bond indices in other '''
            other.D['bonds']['globalIdx']+=bdxshift
        self.N+=len(other.D['atoms'])
        self.metadat['N']=self.N
        self.metadat['nBonds']+=nOtherBonds
        self._myconcat(other,directive='atoms',idxlabel=['globalIdx'],idxshift=idxshift)
        self._myconcat(other,directive='bonds',idxlabel=['ai','aj'],idxshift=idxshift)
        if 'bonds' in other.D:
            self.bondlist.update(other.D['bonds'])
        return (idxshift,bdxshift,rdxshift)

    def _myconcat(self,other,directive,idxlabel=[],idxshift=0):
        if not directive in other.D:
            return
        if directive in self.D:
            for i in idxlabel:
                other.D[directive][i]+=idxshift
            self.D[directive]=pd.concat((self.D[directive],other.D[directive]),ignore_index=True)
        else:
            self.D[directive]=other.D[directive]

    def write_sea(self,filename=''):
        if not 'sea-idx' in self.D['atoms'].columns:
            raise Exception('There is no sea-id in this atoms dataframe')
        if filename=='':
            raise Exception('Please provide a file name to write sea data')
        with open(filename,'w') as f:
            seaformatters = [ lambda x: f'{x:>7d}', lambda x: f'{x:>7d}' ]
            f.write(self.D['atoms'][['globalIdx','sea-idx']].to_string(header=False,index=False,formatters=seaformatters))

    def write_mol2(self,filename=''):
        if self.format!='mol2':
            raise Exception('This config instance is not in mol2 format')
        posscale=1.0
        if self.units['length']=='nm':
            posscale=10.0
        com=self.geometric_center()
        if filename!='':
            atomformatters = [
                lambda x: f'{x:>7d}',
                lambda x: f'{x:<8s}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:>9.4f}',
                lambda x: f'{x:<5s}',
                lambda x: f'{x:>3d}',
                lambda x: f' {x:<7s}',
                lambda x: f'{x:>9.4f}'
            ]
            bondformatters = [
                lambda x: f'{x:>6d}',
                lambda x: f'{x:>5d}',
                lambda x: f'{x:>5d}',
                lambda x: f'{str(x):>4s}'
            ]
            substructureformatters = [
                lambda x: f'{x:>6d}',
                lambda x: f'{x:<7s}',
                lambda x: f'{x:>6d}',
                lambda x: f'{x:<7s}'
            ]
            with open(filename,'w') as f:
                f.write('@<TRIPOS>MOLECULE\n')
                f.write(f'{self.name}\n')
                N=self.N
                # Infer the residue names and resids from the atom records
                rdf=self.D['atoms'][['resNum','resName']].copy().drop_duplicates()
                rdf['rootatom']=[1]*len(rdf)
                rdf['residue']=['RESIDUE']*len(rdf)
                nBonds=self.metadat.get('nBonds',0)
                nSubs=len(rdf)
                nFeatures=self.metadat.get('nFeatures',0)
                nSets=self.metadat.get('nSets',0)
                f.write('{:>6d}{:>6d}{:>3d}{:>3d}{:>3d}\n'.format(N,nBonds,nSubs, nFeatures,nSets))
                f.write(f"{self.metadat.get('mol2type','SMALL')}\n")
                f.write(f"{self.metadat.get('mol2chargetype','GASTEIGER')}\n")
                f.write('\n')
                f.write('@<TRIPOS>ATOM\n')
                if posscale!=1.0:
                    sdf=self.D['atoms'].copy()
                    pos=(sdf.loc[:,['posX','posY','posZ']]-com)*posscale
                    sdf.loc[:,['posX','posY','posZ']]=pos
                    f.write(sdf.to_string(columns=self.mol2_atom_colnames,header=False,index=False,formatters=atomformatters))
                else:
                    f.write(self.D['atoms'].to_string(columns=self.mol2_atom_colnames,header=False,index=False,formatters=atomformatters))
                f.write('\n')
                f.write('@<TRIPOS>BOND\n')
                f.write(self.D['bonds'].to_string(header=False,index=False,formatters=bondformatters))
                f.write('\n')
                ''' write substructure section '''
                f.write('@<TRIPOS>SUBSTRUCTURE\n')
                f.write(rdf.to_string(header=False,index=False,formatters=substructureformatters))
                
    def atomcount(self):
        return self.N

    def cap(self,capping_bonds=[],**kwargs):
        inst=deepcopy(self)
        ''' generate all capping bonds '''
        adf=inst.D['atoms']
        pairs=[]
        orders=[]
        deletes=[]
        for c in capping_bonds:
            ai,aj=c.pairnames
            o=c.bondorder
            idxi=adf[adf['atomName']==ai]['globalIdx'].values[0]
            # TODO: use iloc to get ri
            ri=adf[adf['atomName']==ai][['posX','posY','posZ']].values[0]
            idxj=adf[adf['atomName']==aj]['globalIdx'].values[0]
            # TODO: use iloc to get rj
            rj=adf[adf['atomName']==aj][['posX','posY','posZ']].values[0]
            pairs.append((idxi,idxj))
            orders.append(o)
            idxinidx=inst.bondlist.partners_of(idxi)
            idxjnidx=inst.bondlist.partners_of(idxj)
            inn=[k for k,v in zip(idxinidx,[adf[adf['globalIdx']==i]['atomName'].values[0] for i in idxinidx]) if v.startswith('H')]
            jnn=[k for k,v in zip(idxjnidx,[adf[adf['globalIdx']==i]['atomName'].values[0] for i in idxjnidx]) if v.startswith('H')]
            if len(inn)>0 and len(jnn)>0:
                jdists=[]
                for nj in jnn:
                    # TODO: use iloc -- assuming globalIdx [1...N]
                    rnj=adf[adf['globalIdx']==nj][['posX','posY','posZ']].values[0]
                    r=np.sqrt(((ri-rnj)**2).sum())
                    jdists.append((nj,r))
                jsac=sorted(jdists,key=lambda x: x[1])[0][0]
                deletes.append(jsac)
                idists=[]
                for ni in inn:
                    # TODO: use iloc -- assuming globalIdx [1...N]
                    rni=adf[adf['globalIdx']==ni][['posX','posY','posZ']].values[0]
                    r=np.sqrt(((rj-rni)**2).sum())
                    idists.append((ni,r))
                isac=sorted(idists,key=lambda x: x[1])[0][0]
                deletes.append(isac)
        # print('capping summary:')
        # print('pairs:',pairs)
        # print('orders:',orders)
        # print('deletes:',deletes)
        inst.add_bonds(pairs=pairs,orders=orders)
        inst.delete_atoms(idx=deletes)

        return inst

    def minimum_distance(self,other,self_excludes=[],other_excludes=[]):
        ''' computes the minimum distance between two configurations '''
        sp=self.D['atoms'][~self.D['atoms']['globalIdx'].isin(self_excludes)][['posX','posY','posZ']]
        op=other.D['atoms'][~other.D['atoms']['globalIdx'].isin(other_excludes)][['posX','posY','posZ']]
        minD=1.e9
        for i,srow in sp.iterrows():
            ri=srow.values
            for j,orow in op.iterrows():
                rj=orow.values
                rij=ri-rj
                # print(i,j,rij)
                D=np.sqrt(np.dot(rij,rij))
                if D<minD:
                    minD=D
        return minD

    def rotate(self,R):
        ''' multiplies rotation matrix R by position of each atom '''
        # print('rotate R.shape',R.shape)
        sp=self.D['atoms'][['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            # print(i,ri,ri.shape)
            newri=np.matmul(R,ri)
            # print('rot ri',ri,'newri',newri)
            self.D['atoms'].loc[i,'posX':'posZ']=newri
            # print(sp.loc[i,'posX':'posZ'])

    def translate(self,L):
        ''' translates all atom positions by L '''
        sp=self.D['atoms'][['posX','posY','posZ']]
        for i,srow in sp.iterrows():
            ri=srow.values
            newri=ri+L
            # print('tra ri',ri,'newri',newri)
            self.D['atoms'].loc[i,'posX':'posZ']=newri

    def bond_to(self,other,acc=None,don=None):
        # self.write_mol2(f'TMP-{self.name}-base.mol2')
        ''' creates a new bond from atom acc in self to atom don of other '''
        aadf=self.D['atoms']
        dadf=other.D['atoms']
        # get pre-merge indices of reactive atoms and any H's bound to them
        accidx=aadf[aadf['atomName']==acc]['globalIdx'].values[0]
        accr=aadf[aadf['atomName']==acc][['posX','posY','posZ']].values[0]
        donidx=dadf[dadf['atomName']==don]['globalIdx'].values[0]
        idxaccn=self.bondlist.partners_of(accidx)
        idxdonn=other.bondlist.partners_of(donidx)
        acch=[k for k,v in zip(idxaccn,[aadf[aadf['globalIdx']==i]['atomName'].values[0] for i in idxaccn]) if v.startswith('H')]
        donh=[k for k,v in zip(idxdonn,[dadf[dadf['globalIdx']==i]['atomName'].values[0] for i in idxdonn]) if v.startswith('H')]
        # Try out some bonds:  for each pair of H's (accH,donH), determine
        # the transformation matrix that will align the acc->accH vector 
        # and the donH->don vector, apply this transformation to all atoms in 
        # the donor, make a determination of degree of steric clash
        # find the pair of H's that gives zero steric clash!
        opt_idxstr='None'
        opt_oligidx=()
        overall_maximum=-11.e10
        for accH in acch:
            accHr=aadf[aadf['globalIdx']==accH][['posX','posY','posZ']].values[0]
            accb=accr-accHr
            accb*=1.0/np.linalg.norm(accb)
            for donH in donh:
                idxstr=f'TMP-{self.name}@{accidx}-{accH}+{other.name}@{donidx}-{donH}'
                oligidx=(accidx,accH,donidx,donH)
                donr=dadf[dadf['atomName']==don][['posX','posY','posZ']].values[0]
                donHr=dadf[dadf['globalIdx']==donH][['posX','posY','posZ']].values[0]
                delHr=accr-donHr
                donb=donHr-donr
                donb*=1.0/np.linalg.norm(donb)
                # accb and donb are the two vectors
                cp=np.cross(donb,accb)
                c=np.dot(donb,accb)
                # print(accb,donb,cp)
                v=np.array([[0,-cp[2],cp[1]],[cp[2],0,-cp[0]],[-cp[1],cp[0],0]])
                v2=np.dot(v,v)
                I=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
                # R is the rotation matrix that will rotate donb to align with accb
                R=I+v+v2/(1.+c)
                # rotate translate all donor atoms!
                other.rotate(R)
                donHr=dadf[dadf['globalIdx']==donH][['posX','posY','posZ']].values[0]
                delHr=accr-donHr
                other.translate(delHr)
                fileprefix=hashlib.shake_128(idxstr.encode("utf-8")).hexdigest(8)
                other.write_mol2(fileprefix+'.mol2')
                minD=self.minimum_distance(other,self_excludes=[accH],other_excludes=[donH])
                # print(minD)
                if minD>overall_maximum:
                    overall_maximum=minD
                    opt_idxstr=idxstr
                    opt_oligidx=oligidx
                # print(R)
        # print(f'best config is {opt_idxstr}')
        fileprefix=hashlib.shake_128(opt_idxstr.encode("utf-8")).hexdigest(8)
        take_me=Coordinates.read_mol2(fileprefix+'.mol2')
        other.copy_coords(take_me)
        # other.write_mol2(f'TMP-opt-{other.name}.mol2')
        # merge and grab the idxshift
        shifts=self.merge(other)
        # self.write_mol2(f'TMP-opt-{self.name}+{other.name}.mol2')
        idxshift=shifts[0]
        # add idxshift to don's original index -> you've got the new index!
        accidx,accHidx,donidx_unshifted,donHidx_unshifted=opt_oligidx
        donidx=donidx_unshifted+idxshift
        donHidx=donHidx_unshifted+idxshift
        # add_bond
        # print(f'preparing to add bond {accidx}-{donidx}...')
        self.add_bonds(pairs=[(accidx,donidx)])
        # delete hydrogens (remember to idxshift the one from the donor)
        # print(f'preparing to delete {accHidx} and {donHidx}...')
        self.delete_atoms(idx=[accHidx,donHidx])
        self.name+=f'+{other.name}'
        # this is a ready-to-minimize molecule!
        # exit()

    def add_bonds(self,pairs=[],orders=[]):
        if len(orders)==0:
            orders=[1]*len(pairs)
        ''' add bonds to a set of coordinates
            pairs:  list of 2-tuples of atom global indices '''
        bmi=self.D['bonds'].set_index(['ai','aj']).index
        h=self.mol2_bond_colnames
        for i,(b,o) in enumerate(zip(pairs,orders)):
            ai,aj=b
            # print(f'looking for {ai}-{aj}...')
            if not (ai,aj) in bmi and not (aj,ai) in bmi:
                data=[len(bmi)+i,ai,aj,o]
                # print('adding',data)
                bonddict={k:[v] for k,v in zip(h,data)}
                self.D['bonds']=pd.concat((self.D['bonds'],pd.DataFrame(bonddict)),ignore_index=True)
#        self.D['bonds'].sort_values(by=['ai','aj'],inplace=True)
        self.D['bonds'].globalIdx=self.D['bonds'].index
        self.bondlist=Bondlist.fromDataFrame(self.D['bonds'])
        if 'nBonds' in self.metadat:
            self.metadat['nBonds']=len(self.D['bonds'])

    def delete_atoms(self,idx=[],reindex=True):
        '''
        Deletes atoms whose global indices appear in the list idx.
        If parameter 'reindex' is true, then the global indices 
        are recalculated so that they are sequential starting at 1 with no
        gaps, and two new columns are added to self.DF:
          - 'oldGlobalIdx' contains the global index values before 
             the deletion.
          - 'globalIdxShift' is the change from the old to the new
             global index for each atom.
        '''
        # print('delete_atoms',idx)
        adf=self.D['atoms']
        indexes_to_drop=adf[adf.globalIdx.isin(idx)].index
        indexes_to_keep=set(range(adf.shape[0]))-set(indexes_to_drop)
        self.D['atoms']=adf.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            adf=self.D['atoms']
            oldGI=adf['globalIdx'].copy()
            adf['globalIdx']=adf.index+1
            mapper={k:v for k,v in zip(oldGI,adf['globalIdx'])}
        self.N-=len(idx)
        # print('mapper',mapper)
        ''' delete bonds '''
        # print('delete bonds containing',idx)
        if 'bonds' in self.D:
            d=self.D['bonds']
            # print('bonds before deletion:\n',d.to_string(index=False))
#            print(d.ai.isin(idx),d.aj.isin(idx))
            indexes_to_drop=d[(d.ai.isin(idx))|(d.aj.isin(idx))].index
            # print(self.D['bonds'].iloc[indexes_to_drop].to_string(index=False))
            indexes_to_keep=set(range(d.shape[0]))-set(indexes_to_drop)
            self.D['bonds']=d.take(list(indexes_to_keep)).reset_index(drop=True)
            # print('bonds after deletion:\n',self.D['bonds'].to_string(index=False))
            if reindex:
                d=self.D['bonds']
                d.ai=d.ai.map(mapper)
                d.aj=d.aj.map(mapper)
                d.globalIdx=d.index+1
                # print('bonds after reindexing:\n',self.D['bonds'].to_string(index=False))
            if 'nBonds' in self.metadat:
                self.metadat['nBonds']=len(self.D['bonds'])
            self.bondlist=Bondlist.fromDataFrame(self.D['bonds'])

    def write_gro(self,filename=''):
        ''' write coordinates in Gromacs format '''
        if filename=='':
            raise Exception('write_gro needs a filename')
        has_vel='velX' in self.D['atoms'].columns
        with open(filename,'w') as f:
            f.write(self.title+'\n')
            f.write(f'{self.N:>5d}\n')
            # TODO
            atomformatters = [
                lambda x: f'{x:>5d}',
                lambda x: f'{x:>5s}',
                lambda x: f'{x:>5s}',
                lambda x: f'{x:>5d}']+[lambda x: f'{x:>8.3f}']*6
            if has_vel:
                f.write(self.D['atoms'].to_string(columns=self.gro_colnames),header=False,index=False,formatters=atomformatters)
            else:
                f.write(self.D['atoms'].to_string(columns=self.gro_colnames[:-3]),header=False,index=False,formatters=atomformatters[:-3])
            f.write(f'{self.box[0][0]:10.5f}{self.box[1][1]:10.5f}{self.box[2][2]:10.5f}')
            # output off-diagonals only if at least one of them is non-zero
            x,y=self.box.nonzero()
            if not all(x==y):
                f.write(f'{self.box[0][1]:10.5f}{self.box[0][2]:10.5f}')
                f.write(f'{self.box[1][0]:10.5f}{self.box[1][2]:10.5f}')
                f.write(f'{self.box[2][0]:10.5f}{self.box[2][1]:10.5f}')
"""
MOL2 Format

SUBSTRUCTURE

@<TRIPOS>SUBSTRUCTURE
Each data record associated with this RTI consists of a single data line. The data
line contains the substructure ID, name, root atom of the substructure,
substructure type, dictionary type, chain type, subtype, number of inter
substructure bonds, SYBYL status bits, and user defined comment.
Format:
subst_id subst_name root_atom [subst_type [dict_type
[chain [sub_type [inter_bonds [status
[comment]]]]]]]
• subst_id (integer) = the ID number of the substructure. This is provided
for reference only and is not used by the MOL2 command when reading
the file.
• subst_name (string) = the name of the substructure.
• root_atom (integer) = the ID number of the substructures root atom.
• subst_type - string) = the substructure type: temp, perm, residue, group
or domain.
• dict_type (integer) = the type of dictionary associated with the
substructure.
• chain (string) = the chain to which the substructure belongs (ð 4 chars).
• sub_type (string) = the subtype of the chain.
• inter_bonds (integer) = the number of inter substructure bonds.
• status (string) = the internal SYBYL status bits. These should never be
set by the user. Valid status bit values are LEAF, ROOT, TYPECOL, DICT,
BACKWARD and BLOCK.
• comment (remaining strings on data line) = the comment for the
substructure.
Example:
1 ALA1 1 RESIDUE 1 A ALA 1 ROOT|DICT Comment here
The substructure has 1 as ID, ALA1 as name and atom 1 as root atom. It is
of type RESIDUE and the associated dictionary type is 1 (protein). It is part
of the A chain in the molecule and it is an ALAnine. There is only one inter
substructure bonds. The SYBYL status bits indicate it is the ROOT
substructure of the chain and it came from a dictionary. The comment reads
“Comment here”.
1 ALA1 1
Minimal representation of a substructure.


"""