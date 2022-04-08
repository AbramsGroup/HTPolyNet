'''
coordinates.py -- simple class for handling atom coordinates from gro and mol2 files
'''

import pandas as pd
import numpy as np
from io import StringIO
from shutil import copyfile

from HTPolyNet.bondlist import Bondlist
from HTPolyNet.ambertools import GAFFParameterize
from HTPolyNet.projectfilesystem import Library
from HTPolyNet.gromacs import grompp_and_mdrun

class Coordinates:
    gro_colnames = ['molNum', 'molName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    mol2_atom_colnames = ['globalIdx','atomName','posX','posY','posZ','type','molNum','molName','charge']
    mol2_bond_colnames = ['globalIdx','ai','aj','type']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_colnames, [int, int, int, str])}

    def __init__(self,name=''):
        self.name=name
        self.format=None
        self.title=''
        self.metadat={}
        self.N=0
        ''' dataframes: 'atoms' and 'bonds' '''
        self.D={}
        self.box=np.zeros((3,3))
        
    @classmethod
    def read_gro(cls,filename=''):
        inst=cls(filename)
        inst.format='gro'
        if filename!='':
            with open(filename,'r') as f:
                data=f.read().split('\n')
                while '' in data:
                    data.remove('')
                inst.title=data[0]
                inst.N=int(data[1])
                series={k:[] for k in cls.gro_colnames}
                for x in data[2:-1]:
                    series['molNum'].append(int(x[0:5].strip()))
                    series['molName'].append(x[5:10].strip())
                    series['atomName'].append(x[10:15].strip())
                    series['globalIdx'].append(int(x[15:20].strip()))
                    numbers=list(map(float,[y.strip() for y in x[20:].split()]))
                    if len(numbers)<4:
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
                inst.D['atoms']=pd.DataFrame(series)
                boxdataline=data[-1]
                n=10
                boxdata=list(map(float,[boxdataline[i:i+n].strip() for i in range(0,len(boxdataline),n)]))
                if len(boxdata)<4:
                    inst.box[0][0],inst.box[1][1],inst.box[2][2]=boxdata[0:3]
                if len(boxdata)==9:
                    inst.box[0][1],inst.box[0][2],inst.box[1][0],inst.box[1][2],inst.box[2][0],inst.box[2][1]=boxdata[3:]
        return inst

    @classmethod
    def read_mol2(cls,filename=''):
        inst=cls(filename)
        inst.format='mol2'
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
        for c in ['posX','posY','posZ']:
            self.D['atoms'][c]=other.D['atoms'][c].copy()

    def translate_coords(self,v=[0.,0.,0.]):
        self.D['atoms']['posX']+=v[0]
        self.D['atoms']['posY']+=v[1]
        self.D['atoms']['posZ']+=v[2]

    def geometric_center(self):
        a=self.D['atoms']
        return [a.posX.mean(),a.posY.mean(),a.posZ.mean()]

    def merge(self,other):
        idxshift=0 if 'atoms' not in self.D else len(self.D['atoms'])
        bdxshift=0 if 'bonds' not in self.D else len(self.D['bonds'])
        rdxshift=0 if 'atoms' not in self.D else self.D['atoms'].iloc[-1]['molNum']
        if 'atoms' in other.D:
            other.D['atoms']['molNum']+=rdxshift
        nOtherBonds=0
        if 'bonds' in other.D:
            nOtherBonds=len(other.D['bonds'])
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

    def write_mol2(self,filename=''):
        if self.format!='mol2':
            raise Exception('was this a gro file?')
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
            with open(filename,'w') as f:
                f.write('@<TRIPOS>MOLECULE\n')
                f.write(f'{self.name}\n')
                N=self.N
                nBonds=self.metadat.get('nBonds',0)
                nSubs=self.metadat.get('nSubs',0)
                nFeatures=self.metadat.get('nFeatures',0)
                nSets=self.metadat.get('nSets',0)
                f.write('{:>6d}{:>6d}{:>3d}{:>3d}{:>3d}\n'.format(N,nBonds,nSubs, nFeatures,nSets))
                f.write(f"{self.metadat.get('mol2type','SMALL')}\n")
                f.write(f"{self.metadat.get('mol2chargetype','GASTEIGER')}\n")
                f.write('\n')
                f.write('@<TRIPOS>ATOM\n')
                f.write(self.D['atoms'].to_string(header=False,index=False,formatters=atomformatters))
                f.write('\n')
                f.write('@<TRIPOS>BOND\n')
                f.write(self.D['bonds'].to_string(header=False,index=False,formatters=bondformatters))
                f.write('\n')
                
    def write_gro(self,filename=''):
        if filename!='':
            with open(filename,'w') as f:
                f.write(str(self))

    def atomcount(self):
        return self.N

    def cap(self,capping_bonds=[],**kwargs):
        ''' generate all capping bonds '''
        minimize=kwargs.get('minimize',False)
        adf=self.D['atoms']
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
            idxinidx=self.bondlist.partners_of(idxi)
            idxjnidx=self.bondlist.partners_of(idxj)
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
        self.add_bonds(pairs=pairs,orders=orders)
        self.delete_atoms(idx=deletes)

        return self

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
        self.write_mol2(f'TMP-{self.name}-base.mol2')
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
        overall_maximum=-11.e10
        for accH in acch:
            accHr=aadf[aadf['globalIdx']==accH][['posX','posY','posZ']].values[0]
            accb=accr-accHr
            accb*=1.0/np.linalg.norm(accb)
            for donH in donh:
                idxstr=f'TMP-{self.name}@{accidx}-{accH}+{other.name}@{donidx}-{donH}'
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
                other.write_mol2(idxstr+'.mol2')
                minD=self.minimum_distance(other,self_excludes=[accH],other_excludes=[donH])
                # print(minD)
                if minD>overall_maximum:
                     overall_maximum=minD
                     opt_idxstr=idxstr
                # print(R)
        print(f'best config is {opt_idxstr}')
        take_me=Coordinates.read_mol2(opt_idxstr+'.mol2')
        other.copy_coords(take_me)
        other.write_mol2(f'TMP-opt-{other.name}.mol2')
        # merge and grab the idxshift
        shifts=self.merge(other)
        self.write_mol2(f'TMP-opt-{self.name}+{other.name}.mol2')
        idxshift=shifts[0]
        # add idxshift to don's original index -> you've got the new index!
        opt_bondstrs=opt_idxstr.replace('TMP-','').split('+')
        acc_strs=opt_bondstrs[0].split('@')
        accidx,accHidx=list(map(int,acc_strs[1].split('-')))
        don_strs=opt_bondstrs[1].split('@')
        donidx,donHidx=list(map(int,don_strs[1].split('-')))
        donidx+=idxshift
        donHidx+=idxshift
        # add_bond
        print(f'preparing to add bond {accidx}-{donidx}...')
        self.add_bonds(pairs=[(accidx,donidx)])
        # delete hydrogens (remember to idxshift the one from the donor)
        print(f'preparing to delete {accHidx} and {donHidx}...')
        self.delete_atoms(idx=[accHidx,donHidx])
        # this is a ready-to-minimize molecule!
        # exit()
        pass

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

    def __str__(self):
        '''
        Generates contents of a *.gro file
        '''
        retstr=self.title+'\n'+f'{self.N:>5d}\n'
        has_vel='velX' in self.D['atoms'].columns
        for i,r in self.D['atoms'].iterrows():
            retstr+=f'{r.molNum:5d}{r.molName:>5s}{r.atomName:>5s}{r.globalIdx:5d}{r.posX:8.3f}{r.posY:8.3f}{r.posZ:8.3f}'
            if has_vel:
                retstr+=f'{r.velX:8.4f}{r.velY:8.4f}{r.velZ:8.4f}'
            retstr+='\n'
        retstr+=f'{self.box[0][0]:10.5f}{self.box[1][1]:10.5f}{self.box[2][2]:10.5f}'
        # output off-diagonals only if at least one of them is non-zero
        x,y=self.box.nonzero()
        if not all(x==y):
            retstr+=f'{self.box[0][1]:10.5f}{self.box[0][2]:10.5f}'
            retstr+=f'{self.box[1][0]:10.5f}{self.box[1][2]:10.5f}'
            retstr+=f'{self.box[2][0]:10.5f}{self.box[2][1]:10.5f}'
        retstr+='\n'
        return retstr


