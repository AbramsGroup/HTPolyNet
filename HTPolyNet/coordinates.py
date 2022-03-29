'''
coordinates.py -- simple classes for handling Gromacs commands and data
'''

import pandas as pd
import numpy as np

class Coordinates:
    colnames = ['molNum', 'molName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ', 'velX', 'velY', 'velZ']
    def __init__(self,name=''):
        self.name=name
        self.title=''
        self.N=0
        self.DF=None
        self.box=np.zeros((3,3))
        
    @classmethod
    def fromGroFile(cls,filename=''):
        inst=cls(filename)
        if filename!='':
            with open(filename,'r') as f:
                data=f.read().split('\n')
                while '' in data:
                    data.remove('')
                inst.title=data[0]
                inst.N=int(data[1])
                series={k:[] for k in cls.colnames}
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
                inst.DF=pd.DataFrame(series)
                boxdataline=data[-1]
                n=10
                boxdata=list(map(float,[boxdataline[i:i+n].strip() for i in range(0,len(boxdataline),n)]))
                if len(boxdata)<4:
                    inst.box[0][0],inst.box[1][1],inst.box[2][2]=boxdata[0:3]
                if len(boxdata)==9:
                    inst.box[0][1],inst.box[0][2],inst.box[1][0],inst.box[1][2],inst.box[2][0],inst.box[2][1]=boxdata[3:]
        return inst

    def write_gro(self,filename=''):
        if filename!='':
            with open(filename,'w') as f:
                f.write(str(self))

    def atomcount(self):
        return self.N

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
        indexes_to_drop=self.DF[self.DF.globalIdx.isin(idx)].index
        indexes_to_keep=set(range(self.DF.shape[0]))-set(indexes_to_drop)
        self.DF=self.DF.take(list(indexes_to_keep)).reset_index(drop=True)
        if reindex:
            oldGI=self.DF['globalIdx'].copy()
            self.DF['oldGlobalIdx']=oldGI
            self.DF['globalIdx']=self.DF.index+1
            self.DF['globalIdxShift']=self.DF['globalIdx']-oldGI
        self.N-=len(idx)
    
    def __str__(self):
        '''
        Generates contents of a *.gro file
        '''
        retstr=self.title+'\n'+f'{self.N:>5d}\n'
        has_vel='velX' in self.DF.columns
        for i,r in self.DF.iterrows():
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


