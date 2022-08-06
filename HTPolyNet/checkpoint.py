import os
import yaml
import logging
from HTPolyNet.topocoord import TopoCoord

logger=logging.getLogger(__name__)

class Checkpoint:
    reqd_attr=['stepschecked','substepschecked','currentstepname']
    reqd_files=['top','gro','grx']
    opt_files=['mol2']
    default_filename='checkpoint_state.yaml'
    def __init__(self,filename=default_filename):
        self.files={}
        for v in self.reqd_files:
            self.files[v]=''
        for v in self.opt_files:
            self.files[v]=''
        self.stepschecked=[]
        self.substepschecked=[]
        self.currentstepname='none'
        self.my_filename=filename
    
    def to_yaml(self):
        with open(self.my_filename,'w') as f:
            f.write(yaml.dump(self))

    @classmethod
    def from_yaml(cls,filename='checkpoint_state.yaml'):
        with open(filename,'r') as f:
            yaml_string=f.read()
        inst=yaml.load(yaml_string,Loader=yaml.Loader)
        return inst

    @classmethod
    def read(cls,TC:TopoCoord,filename):
        inst=cls.from_yaml(filename)
        inst.my_filename=os.path.abspath(filename)
        for ext in cls.reqd_files:
            TC.files[ext]=inst.files[ext]
        for ext in cls.opt_files:
            TC.files[ext]=inst.files[ext]
        return inst

    @classmethod
    def setup(cls,TC:TopoCoord,filename='checkpoint_state.yaml'):
        if not os.path.exists(filename):
            logger.debug('New checkpoint; no file')
            absfilename=os.path.join(os.getcwd(),filename)
            assert not os.path.exists(absfilename)
            return cls(filename=absfilename)
        return cls.read(TC,filename)

    def set(self,TC:TopoCoord,stepname):
        # logger.debug(f'{stepname} {TC.files}')
        self.stepschecked.append(stepname)
        self.substepschecked=[]
        self.currentstepname='none'
        self._filename_transfers(TC)
        self.to_yaml()

    def _filename_transfers(self,TC:TopoCoord):
        for v in self.reqd_files:
            self.files[v]=os.path.abspath(TC.files[v])
        for v in self.opt_files:
            if v in TC.files and TC.files[v]:
                self.files[v]=os.path.abspath(TC.files[v])

    def set_substep(self,TC:TopoCoord,currentstepname,substepname):
        self.currentstepname=currentstepname
        self.substepschecked.append(substepname)
        self._filename_transfers(TC)
        self.to_yaml()

    def passed(self,stepname):
        return stepname in self.stepschecked
    
    def is_currentstepname(self,stepname):
        return self.currentstepname==stepname

    def last_substep(self):
        return self.substepschecked[-1]

_CP_:Checkpoint=None
def setup(TC:TopoCoord,filename='checkpoint_state.yaml'):
    global _CP_
    _CP_=Checkpoint.setup(TC,filename)

def set(TC:TopoCoord,stepname):#,state:CPstate):
    global _CP_
    _CP_.set(TC,stepname)

def subset(TC:TopoCoord,currentstepname,substepname):
    global _CP_
    _CP_.set_substep(TC,currentstepname,substepname)

def passed(stepname):
    return _CP_.passed(stepname)

def is_currentstepname(testname):
    return _CP_.is_currentstepname(testname)

def last_substep():
    return _CP_.last_substep()