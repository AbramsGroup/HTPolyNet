import logging
import os
import yaml
import pandas as pd
from enum import Enum
class CPstate(Enum):
    fresh=0 # nothing there at all
    bondsearch_complete=1 # a bondsfile exists and bonds_are=='unrelaxed'
    relaxing=3 # 
    equilibrating=4
    finished=5

class CURECheckpoint:
    def __init__(self,iter,checkpoint_file='complete.yaml',filename_format='cure-{iter}:{radidx}-{stage}',bonds_file='bonds.csv'):
        # self.n_relaxed_bonds=0
        # self.n_unrelaxed_bonds=0
        self.state=CPstate.fresh
        self.iter=iter
        self.current_stage=0
        self.current_radidx=0
        self.checkpoint_file=checkpoint_file
        self.filename_format=filename_format
        self.bonds_file=bonds_file
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def read_checkpoint(self):
        self.current_stage=0
        self.current_radidx=0
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file,'r') as f:
                basedict=yaml.safe_load(f)
                self.top=os.path.basename(basedict['TOPOLOGY'])
                self.gro=os.path.basename(basedict['COORDINATES'])
                self.grx=os.path.basename(basedict['EXTRA_ATTRIBUTES'])
                self.current_stage=basedict['CURRENT_STAGE']
                self.current_radidix=basedict['CURRENT_RADIDX']
                bf=basedict.get('BONDSFILE',None)
                self.bonds_are=basedict.get('BONDS_ARE',None)
                if bf:
                    self.bondsfile=os.path.basename(bf)
                    self.read_bondsfile()

    def set_state(self):
        # TODO: fix this
        if self.bonds.shape[0]>0: # there are new bonds
            if self.bonds_are=='unrelaxed':
                self.state=CPstate.bondsearch_complete
                if self.current_stage>0:
                    self.state=CPstate.relaxing
            elif self.bonds_are=='relaxed':
                self.state=CPstate.finished
        else:
            self.state=CPstate.fresh

    def write_bondsfile(self):
        self.bonds.to_csv(self.bondsfile,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def read_bondsfile(self):
        assert os.path.exists(self.bondsfile),f'Error: {self.bondsfile} not found.'
        self.bonds=pd.read_csv(self.bondsfile,sep='\s+',header=0)

    def register_bonds(self,bonds,bonds_are='unrelaxed'):
        self.bonds=pd.DataFrame({'ai':[x[0] for x in bonds],'aj':[x[1] for x in bonds]})
        self.bonds_are=bonds_are
    
    def finished(self):
        if self.bonds_are=='relaxed':
            logging.info(f'CURE iteration {self.iter} marked complete:')
            logging.info(f'        topology    {self.top}')
            logging.info(f'        coordinates {self.gro}')
            logging.info(f'        extra_attr  {self.grx}')
            logging.info(f'        #newbonds   {len(self.bonds)}')
        return self.bonds_are=='relaxed'

    def write_checkpoint(self,radidx,stage,bonds,bonds_are='unrelaxed'):
        self.top=self.filename_format.format(iter,radidx,stage)+'.top'
        self.gro=self.filename_format.format(iter,radidx,stage)+'.gro'
        self.grx=self.filename_format.format(iter,radidx,stage)+'.grx'
        self.current_radidx=radidx
        self.current_stage=stage
        with open('complete.yaml','w') as f:
            f.write(f'ITERATION: {self.iter}\n')
            f.write(f'CURRENT_STAGE: {self.current_stage}\n')
            f.write(f'CURRENT_RADIDX: {self.current_radidx}\n')
            f.write(f'TOPOLOGY: {self.top}\n')
            f.write(f'COORDINATES: {self.gro}\n')
            f.write(f'EXTRA_ATTRIBUTES: {self.grx}\n')
            f.write(f'NUM_NEWBONDS: {num_newbonds}\n')
            f.close()
