import os
import yaml
import logging
import functools
# from HTPolyNet.topocoord import TopoCoord
# import HTPolyNet.projectfilesystem as pfs

logger=logging.getLogger(__name__)

class Checkpoint:
    # reqd_attr=['stepschecked','substepschecked','currentstepname']
    # reqd_files=['top','gro','grx']
    # opt_files=['mol2']
    default_filename='checkpoint_state.yaml'
    def __init__(self,input_dict={}):
        self.my_abspath=os.getcwd()
        self.calls:list=input_dict.get('calls',[])
        self.results:dict=input_dict.get('results',{})
        self.narrative:list=input_dict.get('narrative',[])

    def to_yaml(self):
        with open(self.default_filename,'w') as f:
            f.write(yaml.dump(self))

    @classmethod
    def from_yaml(cls):
        try:
            with open(cls.default_filename,'r') as f:
                yaml_string=f.read()
            inst=yaml.load(yaml_string,Loader=yaml.Loader)
            # overwrite the absolute path in the file upon reading in
            inst.my_abspath=os.getcwd()
            return inst
        except FileNotFoundError:
            return cls()
    # @classmethod
    # def read(cls,TC:TopoCoord,filename):
    #     inst=cls.from_yaml(filename)
    #     inst.my_filename=filename
    #     for ext in cls.reqd_files:
    #         TC.files[ext]=os.path.abspath(inst.files[ext])
    #     for ext in cls.opt_files:
    #         TC.files[ext]=os.path.abspath(inst.files[ext])
    #     return inst

    # @classmethod
    # def setup(cls,TC:TopoCoord,filename='checkpoint_state.yaml'):
    #     if not os.path.exists(filename):
    #         logger.debug('New checkpoint; no file')
    #         absfilename=os.path.join(os.getcwd(),filename)
    #         assert not os.path.exists(absfilename)
    #         return cls(filename=filename)
    #     return cls.read(TC,filename)

#     def set(self,TC:TopoCoord,stepname):
#         # logger.debug(f'{stepname} {TC.files}')
#         self.stepschecked.append(stepname)
#         self.substepschecked=[]
#         self.currentstepname='none'
#         self._filename_transfers(TC)
#         self.to_yaml()

#     def _filename_transfers(self,TC:TopoCoord):
#         for v in self.reqd_files:
#             self.files[v]=pfs.proj_abspath(TC.files[v])
#         for v in self.opt_files:
#             if v in TC.files and TC.files[v]:
#                 self.files[v]=pfs.proj_abspath(TC.files[v])

#     def set_substep(self,TC:TopoCoord,currentstepname,substepname):
#         self.currentstepname=currentstepname
#         self.substepschecked.append(substepname)
#         self._filename_transfers(TC)
#         self.to_yaml()

#     def passed(self,stepname):
#         if stepname in self.stepschecked:
#             logger.debug(f'[ {stepname} ] checkpoint passed')
#             return True
#         return False
    
#     def is_currentstepname(self,stepname):
#         return self.currentstepname==stepname

#     def last_substep(self):
#         return self.substepschecked[-1]

# _CP_:Checkpoint=None
# def setup(TC:TopoCoord,filename='checkpoint_state.yaml'):
#     global _CP_
#     _CP_=Checkpoint.setup(TC,filename)

# def set(TC:TopoCoord,stepname):#,state:CPstate):
#     global _CP_
#     _CP_.set(TC,stepname)

# def subset(TC:TopoCoord,currentstepname,substepname):
#     global _CP_
#     _CP_.set_substep(TC,currentstepname,substepname)

# def passed(stepname):
#     return _CP_.passed(stepname)

# def is_currentstepname(testname):
#     return _CP_.is_currentstepname(testname)

# def last_substep():
#     return _CP_.last_substep()

_CP_=Checkpoint()
# class EnableCheckpointClass:
#     def __init__(self,func):
#         functools.update_wrapper(self,func)
#         self.func=func

#     def __call__(self,*args,**kwargs):
#         mywd=os.path.relpath(os.getcwd(),_CP_.my_abspath)
#         if len(_CP_.calls)>0 and (self.func.__name__,mywd) in _CP_.calls:
#             print(f'skipping {self.func.__name__}') 
#             return
#         result=self.func(*args,**kwargs)
#         _CP_.calls.append((self.func.__name__,mywd))
#         _CP_.results.append(result)
#         _CP_.narrative.append(f'Function {self.func.__name__} called with args {args} in {os.path.join(_CP_.my_abspath,mywd)} gave result {result}')
#         return result

def enableCheckpoint(method):
    @functools.wraps(method)
    def wrapper_method(self,*args,**kwargs):
        mywd=os.path.relpath(os.getcwd(),_CP_.my_abspath)
        if len(_CP_.calls)>0 and (method.__name__,mywd) in _CP_.calls:
            logger.info(f'Skipping {method.__name__} in {mywd}') 
            return
        result=method(self,*args,**kwargs)
        _CP_.calls.append((method.__name__,mywd))
        _CP_.results.update(result) # must be a dict
        _CP_.narrative.append(f'Method {method.__name__} called in {os.path.join(_CP_.my_abspath,mywd)} gave result {result}')
        _write_checkpoint()
        return result
    return wrapper_method

def _write_checkpoint():
    sv=os.getcwd()
    os.chdir(_CP_.my_abspath)
    _CP_.to_yaml()
    os.chdir(sv)

def read_checkpoint():
    global _CP_
    _CP_=Checkpoint.from_yaml()
    if len(_CP_.calls)>0:
        lwd=os.path.join(_CP_.my_abspath,_CP_.calls[-1][1])
    return {c:os.path.join(lwd,x) for c,x in _CP_.results.items()}