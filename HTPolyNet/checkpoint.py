import logging
import os
import yaml
import pandas as pd
from enum import Enum
class CPstate(Enum):
    fresh=0 # nothing there at all
    generated_templates=1
    generated_initial_topology=2
    generated_initial_coordinates=3
    initial_equilibration=4
    bondsearch=5
#    bondsearch_complete=6
    update=7
    relax_prestage=8
    relax_poststage=9 
    equilibrate=10
    post_equilibration=11
    finished=12
    unknown=99
    def __str__(self):
        return self.name

class Checkpoint:
    def __init__(self,checkpoint_file='checkpoint.yaml',filename_format='cure-{iter}-{stage}',bonds_file='bonds.csv',n_stages=10):
        self.state=CPstate.fresh
        self.iter=0
        self.current_stage=0
        self.current_radidx=0
        self.checkpoint_file=checkpoint_file
        self.filename_format=filename_format
        self.n_stages=n_stages
        self.bonds_file=bonds_file
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def set_state(self,state):
        self.state=state

    def write_bondsfile(self):
        self.bonds.to_csv(self.bonds_file,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def read_bondsfile(self):
        assert os.path.exists(self.bonds_file),f'Error: {self.bonds_file} not found.'
        self.bonds=pd.read_csv(self.bonds_file,sep='\s+',header=0)

    def register_bonds(self,bonds,bonds_are='unrelaxed'):
        self.bonds=bonds
        self.bonds_are=bonds_are
        self.write_bondsfile()
    
    def finished(self):
        if self.bonds_are=='relaxed':
            logging.info(f'CURE iteration {self.iter} marked complete:')
            logging.info(f'        topology    {self.top}')
            logging.info(f'        coordinates {self.gro}')
            logging.info(f'        extra_attr  {self.grx}')
            logging.info(f'        #newbonds   {len(self.bonds)}')
        return self.bonds_are=='relaxed'

    def read_checkpoint(self,system):
        self.current_stage=0
        self.current_radidx=0
        self.state=CPstate.fresh
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file,'r') as f:
                basedict=yaml.safe_load(f)
                self.iter=basedict['ITERATION']
                self.state=CPstate[basedict['STATE']]
                self.top=os.path.basename(basedict['TOPOLOGY'])
                self.gro=os.path.basename(basedict['COORDINATES'])
                self.grx=os.path.basename(basedict['EXTRA_ATTRIBUTES'])
                self.current_stage=basedict['CURRENT_STAGE']
                self.current_radidx=basedict['CURRENT_RADIDX']
                bf=basedict.get('BONDSFILE',None)
                self.bonds_are=basedict.get('BONDS_ARE',None)
                if bf:
                    self.bonds_file=os.path.basename(bf)
                    self.read_bondsfile()
                system.set_system(CP=self)

    def write_checkpoint(self,system,state,prefix='checkpoint'):
        self.state=state
        self.top,self.gro,self.grx=[prefix+x for x in ['.top','.gro','.grx']]
        # system.TopoCoord.Topology.null_check(msg='writing checkpoint')
        system.register_system(CP=self)
        with open(self.checkpoint_file,'w') as f:
            f.write(f'ITERATION: {self.iter}\n')
            f.write(f'STATE: {str(self.state)}\n')
            f.write(f'CURRENT_STAGE: {self.current_stage}\n')
            f.write(f'CURRENT_RADIDX: {self.current_radidx}\n')
            f.write(f'TOPOLOGY: {self.top}\n')
            f.write(f'COORDINATES: {self.gro}\n')
            f.write(f'EXTRA_ATTRIBUTES: {self.grx}\n')
            if self.bonds.shape[0]>0:
                f.write(f'BONDS_ARE: {self.bonds_are}\n')
                f.write(f'BONDSFILE: {self.bonds_file}\n')
            f.close()
        self.write_bondsfile()