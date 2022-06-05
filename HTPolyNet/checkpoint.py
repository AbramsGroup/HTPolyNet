import os
import yaml
import pandas as pd
from enum import Enum
class CPstate(Enum):
    """Enumerated CURE state
    """
    fresh=0
    generated_templates=1
    generated_initial_topology=2
    generated_initial_coordinates=3
    initial_equilibration=4
    bondsearch=5
    drag=6
    update=7
    relax=8
    equilibrate=10
    post_equilibration=11
    postcure=12
    postcure_equilibration=13
    finished=25
    unknown=99
    def __str__(self):
        return self.name

class Checkpoint:
    """Manages checkpointing in the CURE algorithm
    """
    def __init__(self,checkpoint_file='checkpoint.yaml'):
        self.state=CPstate.fresh
        self.iter=0
        self.current_dragstage=0
        self.current_stage=0
        self.current_radidx=0
        self.radius=0.0
        self.checkpoint_file=checkpoint_file
        self.bonds_file=None
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def set_state(self,state):
        self.state=state

    def reset_for_next_iter(self,checkpoint_file='checkpoint.yaml'):
        self.iter+=1
        self.state=CPstate.fresh
        self.current_dragstage=0
        self.current_stage=0
        self.current_radidx=0
        self.radius=0.0
        self.checkpoint_file=checkpoint_file
        self.bonds_file=None
        self.bonds=pd.DataFrame()
        self.bonds_are='nonexistent!'

    def _write_bondsfile(self):
        self.bonds.to_csv(self.bonds_file,sep=' ',mode='w',index=False,header=True,doublequote=False)

    def _read_bondsfile(self):
        assert os.path.exists(self.bonds_file),f'Error: {self.bonds_file} not found.'
        self.bonds=pd.read_csv(self.bonds_file,sep='\s+',header=0)

    def register_bonds(self,bonds,pairs,bonds_are='unrelaxed'):
        self.bonds=bonds
        self.pairs=pairs
        self.bonds_are=bonds_are
    
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
            self.current_dragstage=basedict['CURRENT_DRAGSTAGE']
            self.current_stage=basedict['CURRENT_STAGE']
            self.current_radidx=basedict['CURRENT_RADIDX']
            self.radius=basedict['RADIUS']
            bf=basedict.get('BONDSFILE',None)
            assert bf,f'Error: BONDSFILE not found in {self.checkpoint_file}.'
            self.bonds_are=basedict.get('BONDS_ARE',None)
            self.bonds_file=os.path.basename(bf)
            self._read_bondsfile()
            system.set_system(CP=self)

    def write_checkpoint(self,system,state,prefix='checkpoint'):
        self.state=state
        self.top,self.gro,self.grx=[prefix+x for x in ['.top','.gro','.grx']]
        self.bonds_file=prefix+'-bonds.csv'
        system.register_system(CP=self)
        with open(self.checkpoint_file,'w') as f:
            f.write(f'ITERATION: {self.iter}\n')
            f.write(f'STATE: {str(self.state)}\n')
            f.write(f'CURRENT_DRAGSTAGE: {self.current_dragstage}\n')
            f.write(f'CURRENT_STAGE: {self.current_stage}\n')
            f.write(f'CURRENT_RADIDX: {self.current_radidx}\n')
            f.write(f'RADIUS: {self.radius}\n')
            f.write(f'TOPOLOGY: {self.top}\n')
            f.write(f'COORDINATES: {self.gro}\n')
            f.write(f'EXTRA_ATTRIBUTES: {self.grx}\n')
            if self.bonds.shape[0]>0:
                f.write(f'BONDS_ARE: {self.bonds_are}\n')
                f.write(f'BONDSFILE: {self.bonds_file}\n')
            f.close()
        self._write_bondsfile()

_CP_=Checkpoint()
def checkpoint_setup(checkpoint_file='checkpoint.yaml'):
    global _CP_
    _CP_.checkpoint_file=checkpoint_file

def checkpoint_write(system,state,prefix='checkpoint'):
    global _CP_
    _CP_.write_checkpoint(system,state,prefix=prefix)

def checkpoint_read(system):
    global _CP_
    _CP_.read_checkpoint(system)
