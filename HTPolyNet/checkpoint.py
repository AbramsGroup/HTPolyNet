import os
import yaml
import pandas as pd
from enum import Enum
import logging
from HTPolyNet.topocoord import TopoCoord

logger=logging.getLogger(__name__)



class Checkpoint:
    def __init__(self):
        self.checkpoint_file=''
        self.top=None
        self.gro=None
        self.grx=None
        self.mol2=None
        self.cwd=os.getcwd()
        self.stepschecked=[]


 
    def read_checkpoint(self): #,system):
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file,'r') as f:
                basedict=yaml.safe_load(f)
            self.top=os.path.basename(basedict['topology'])
            self.gro=os.path.basename(basedict['coordinates'])
            self.grx=os.path.basename(basedict['extra_attributes'])
            self.mol2=basedict.get('mol2_coordinates',None)
            if self.mol2:
                self.mol2=os.path.basename(basedict['mol2_coordinates'])
            self.cwd=os.path.commonpath([basedict['topology'],basedict['coordinates']])
            os.chdir(self.cwd)
            # self.iter=basedict['ITERATION']
            # self.state=CPstate[basedict['STATE']]
            # self.current_dragstage=basedict['CURRENT_DRAGSTAGE']
            # self.current_stage=basedict['CURRENT_STAGE']
            # self.current_radidx=basedict['CURRENT_RADIDX']
            # self.radius=basedict['RADIUS']
            # bf=basedict.get('BONDSFILE',None)
            # # assert bf,f'Error: BONDSFILE not found in {self.checkpoint_file}.'
            # self.bonds_are=basedict.get('BONDS_ARE',None)
            # if bf:
            #     self.bonds_file=os.path.basename(bf)
            #     self._read_bondsfile()
            # else:
                # self.bonds_file=None
            #system.set_system(CP=self)
        # else:
        #     logging.debug(f'read_checkpoint: no file, empty checkpoint')

    def write_checkpoint(self,TC:TopoCoord,stepname): #,state): #system,state,prefix='checkpoint'):
        # self.state=state
        self.stepschecked.append(stepname)
        self.top,self.gro,self.grx,self.mol2=[os.path.basename(TC.files[x]) for x in ['top','gro','grx','mol2']]
        if self.top:
            TC.write_top(self.top)
        if self.gro:
            TC.write_gro(self.gro)
        if self.grx:
            TC.write_grx_attributes(self.grx)
        if self.mol2:
            TC.write_mol2(self.mol2)
        self.cwd=os.getcwd()
        # self.bonds_file=prefix+'-bonds.csv'
        # system.register_system(CP=self)
        with open(self.checkpoint_file,'w') as f:
            # f.write(f'ITERATION: {self.iter}\n')
            # f.write(f'STATE: {str(self.state)}\n')
            # f.write(f'CURRENT_DRAGSTAGE: {self.current_dragstage}\n')
            # f.write(f'CURRENT_STAGE: {self.current_stage}\n')
            # f.write(f'CURRENT_RADIDX: {self.current_radidx}\n')
            # f.write(f'RADIUS: {self.radius}\n')
            f.write(f'stepschecked: {self.stepschecked}\n')
            f.write(f'cwd: {self.cwd}\n')
            if self.top:
                f.write(f'topology: {self.top}\n')
            if self.gro:
                f.write(f'coordinates: {self.gro}\n')
            if self.grx:
                f.write(f'extra_attributes: {self.grx}\n')
            if self.mol2:
                f.write(f'mol2_coordinates: {self.mol2}\n')
            # if self.bonds.shape[0]>0:
            #     f.write(f'BONDS_ARE: {self.bonds_are}\n')
            #     f.write(f'BONDSFILE: {self.bonds_file}\n')
            f.close()
        # self._write_bondsfile()

_CP_=Checkpoint()
def setup(checkpoint_file='checkpoint.yaml'):
    global _CP_
    if os.path.exists(checkpoint_file):
        _CP_.checkpoint_file=os.path.abspath(checkpoint_file)
    else:
        _CP_.checkpoint_file=os.path.join(os.getcwd(),checkpoint_file)
    logger.debug(f'checkpoints to {_CP_.checkpoint_file}')

def write(TC:TopoCoord,stepname):#,state:CPstate):
    global _CP_
    _CP_.write_checkpoint(TC,stepname)

def read():
    global _CP_
    _CP_.read_checkpoint()
    return TopoCoord(topfilename=_CP_.top,grofilename=_CP_.gro,grxfilename=_CP_.grx,mol2filename=_CP_.mol2)

def passed(stepname):
    return stepname in _CP_.stepschecked
