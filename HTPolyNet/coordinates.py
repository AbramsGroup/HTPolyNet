'''
coordinates.py -- simple class for handling atom coordinates from gro and mol2 files
'''

import pandas as pd
import numpy as np
from io import StringIO

from HTPolyNet.bondlist import Bondlist
from HTPolyNet.ambertools import Command

class Coordinates:
    gro_colnames = ['molNum', 'molName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    mol2_atom_colnames = ['globalIdx','atomName','posX','posY','posZ','type','molNum','molName','charge']
    mol2_bond_colnames = ['globalIdx','ai','aj','type']
    mol2_bond_types = {k:v for k,v in zip(mol2_bond_colnames, [int, int, int, str])}

    def __init__(self,name=''):
        self.name=name
        self.format=None
        self.title=''
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
            metadat=sections['molecule'][1].strip().split()
            imetadat=list(map(int,metadat))
            inst.N=imetadat[0]
            inst.nBonds=imetadat[1]
            inst.nSubs=imetadat[2]
            inst.nFeatures=imetadat[3]
            inst.nSets=imetadat[4]
            inst.mol2type=sections['molecule'][2]
            inst.mol2chargetype=sections['molecule'][3]
            inst.D['atoms']=pd.read_csv(sections['atom'],sep='\s+',names=Coordinates.mol2_atom_colnames)
            inst.N=len(inst.D['atoms'])
            inst.D['bonds']=pd.read_csv(sections['bond'],sep='\s+',names=Coordinates.mol2_bond_colnames,dtype=Coordinates.mol2_bond_types)
            inst.bondlist=Bondlist.fromDataFrame(inst.D['bonds'])
        return inst

    def write_mol2(self,filename=''):
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
                f.write('{:>3d}{:>3d}{:>2d}{:>2d}{:>2d}\n'.format(self.N,self.nBonds,self.nSubs,self.nFeatures,self.nSets))
                f.write(f'{self.mol2type}\n')
                f.write(f'{self.mol2chargetype}\n')
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
        logf=kwargs.get('logf',None)
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
        print('capping summary:')
        print('pairs:',pairs)
        print('orders:',orders)
        print('deletes:',deletes)
        self.add_bonds(pairs=pairs,orders=orders)
        self.delete_atoms(idx=deletes)
        self.write_mol2(f'{self.name}-capped-unminimized.mol2')
        # cmd1 = 'obabel {}.mol2 -O {}.mol2 --minimize --sd --c 1e-5'
        c=Command(f'obabel {self.name}-capped-unminimized.mol2 --minimize --sd --c 1.e-5',O=f'{self.name}-capped-minimized.mol2')
        msg=c.run()
        if logf:
            logf(msg)
        self.read_mol2(f'{self.name}-capped-minimized.mol2')
        return self

    def add_bonds(self,pairs=[],orders=[]):
        ''' add bonds to a set of coordinates
            pairs:  list of 2-tuples of atom global indices '''
        bmi=self.D['bonds'].set_index(['ai','aj']).index
        h=self.mol2_bond_colnames
        for i,(b,o) in enumerate(zip(pairs,orders)):
            ai,aj=b
            if not (ai,aj) in bmi and not (aj,ai) in bmi:
                data=[len(bmi)+i,ai,aj,o]
                # print('adding',data)
                bonddict={k:[v] for k,v in zip(h,data)}
                self.D['bonds']=pd.concat((self.D['bonds'],pd.DataFrame(bonddict)),ignore_index=True)
#        self.D['bonds'].sort_values(by=['ai','aj'],inplace=True)
        self.bondlist=Bondlist.fromDataFrame(self.D['bonds'])

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
                # print('bonds after reindexing:\n',self.D['bonds'].to_string(index=False))
            self.nBonds=len(self.D['bonds'])
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


