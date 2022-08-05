import os
import yaml
import logging
from HTPolyNet.topocoord import TopoCoord

logger=logging.getLogger(__name__)

class Checkpoint:
    reqd_attr=['stepschecked','substepschecked','currentstepname']
    reqd_files={'topology':'top','coordinates':'gro','extra_attributes':'grx'}
    opt_files={'mol2_coordinates':'mol2'}
    def __init__(self):
        self.checkpoint_file=''
        self.files={}
        self.stepschecked=[]
        self.substepschecked=[]
        self.currentstepname=None
 
    def read_checkpoint(self,TC:TopoCoord):
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file,'r') as f:
                basedict=yaml.safe_load(f)
            all_there=all([x in basedict for x in self.reqd_attr+list(self.reqd_files.keys())])
            assert all_there,f'{self.checkpoint_file} is missing one or more of {self.reqd_attr+list(self.reqd_files.keys())}'
            all_abs=all([basedict[filetype]==os.path.abspath(basedict[filetype]) for filetype in self.reqd_files.keys()])
            assert all_abs,f'{self.checkpoint_file} should contain absolute path names'
            self.stepschecked=basedict['stepschecked']
            self.substepschecked=basedict['substepschecked']
            self.currentstepname=basedict['currentstepname']
            for filetype,ext in self.reqd_files.items():
                self.files[filetype]=basedict[filetype]
                TC.files[ext]=basedict[filetype]
            for filetype,ext in self.opt_files.items():
                filename=basedict.get(filetype,None)
                if filename:
                    assert os.path.abspath(filename)==filename,f'{self.checkpoint_file} should contain absolute path names'
                    self.files[filetype]=filename
                    TC.files[ext]=basedict[filetype]
            return True
        return False

    def set_checkpoint(self,TC:TopoCoord,stepname):
        logger.debug(f'{stepname} {TC.files}')
        self.stepschecked.append(stepname)
        self.substepschecked=[]
        self.currentstepname='none'
        self.write(TC)

    def set_subcheckpoint(self,TC:TopoCoord,currentstepname,substepname):
        self.currentstepname=currentstepname
        self.substepschecked.append(substepname)
        self.write(TC)

    def write(self,TC:TopoCoord):
        for k,v in self.reqd_files.items():
            self.files[k]=os.path.abspath(TC.files[v])
        for k,v in self.opt_files.items():
            if v in TC.files and TC.files[v]:
                self.files[k]=os.path.abspath(TC.files[v])
        with open(self.checkpoint_file,'w') as f:
            f.write(f'stepschecked: {self.stepschecked}\n')
            f.write(f'substepschecked: {self.substepschecked}\n')
            f.write(f'currentstepname: {self.currentstepname}\n')
            for k,v in self.reqd_files.items():
                f.write(f'{k}: {self.files[k]}\n')
            for k,v in self.opt_files.items():
                if k in self.files:
                    f.write(f'{k}: {self.files[k]}\n')
            f.close()
    

_CP_=Checkpoint()
def setup(TC:TopoCoord,checkpoint_file='checkpoint.yaml'):
    global _CP_
    if os.path.exists(checkpoint_file):
        _CP_.checkpoint_file=os.path.abspath(checkpoint_file)
        logger.debug(f'Reading existing checkpoint file {checkpoint_file}')
        return _CP_.read_checkpoint(TC)
    else:
        logger.debug(f'Checkpoints to {_CP_.checkpoint_file}')
        _CP_.checkpoint_file=os.path.join(os.getcwd(),checkpoint_file)
        return False

def set(TC:TopoCoord,stepname):#,state:CPstate):
    global _CP_
    _CP_.set_checkpoint(TC,stepname)

def subset(TC:TopoCoord,currentstepname,substepname):
    global _CP_
    _CP_.set_subcheckpoint(TC,currentstepname,substepname)

# def read():
#     global _CP_
#     return _CP_.read_checkpoint()

def passed(stepname):
    return stepname in _CP_.stepschecked

def is_currentstepname(testname):
    return _CP_.currentstepname==testname

def last_substep():
    return _CP_.substepschecked[-1]