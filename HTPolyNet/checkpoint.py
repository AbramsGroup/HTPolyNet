"""

.. module:: checkpoint
   :synopsis: Implements a simple checkpointing scheme using a wrapper
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import yaml
import logging
import functools

logger=logging.getLogger(__name__)

class Checkpoint:
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

_CP_=Checkpoint()

def enableCheckpoint(method):
    """enableCheckpoint wraps any method so that every call is registered in a history of calls in a written checkpoint file

    :param method: name of method to be wrapped
    :type method: method
    :return: wrapped method
    :rtype: method
    """
    @functools.wraps(method)
    def wrapper_method(self,*args,**kwargs):
        ''' define working directory as current directory relative to the checkpoint '''
        mywd=os.path.relpath(os.getcwd(),_CP_.my_abspath)
        ''' if the method calling this from the current directory is already in the checkpoint history, 
        to not call the method, just return '''
        if len(_CP_.calls)>0 and (method.__name__,mywd) in _CP_.calls:
            logger.info(f'Skipping {method.__name__} in {mywd}') 
            return
        ''' call the method, save result '''
        result=method(self,*args,**kwargs)
        ''' register this method call in this directory and its results '''
        _CP_.calls.append((method.__name__,mywd))
        _CP_.results.update(result) # must be a dict
        _CP_.narrative.append(f'Method {method.__name__} called in {os.path.join(_CP_.my_abspath,mywd)} gave result {result}')
        ''' update the written checkpoint file '''
        _write_checkpoint()
        return result
    return wrapper_method

def _write_checkpoint():
    """_write_checkpoint writes the checkpoint file in its globally resolved location
    """
    sv=os.getcwd()
    os.chdir(_CP_.my_abspath)
    _CP_.to_yaml()
    os.chdir(sv)

def read_checkpoint():
    """read_checkpoint creates a new global Checkpoint object by reading from the default file

    :return: current results dictionary with any pathnames resolved as absolute
    :rtype: dict
    """
    global _CP_
    _CP_=Checkpoint.from_yaml()
    if len(_CP_.calls)>0:
        lwd=os.path.join(_CP_.my_abspath,_CP_.calls[-1][1])
    return {c:os.path.join(lwd,x) for c,x in _CP_.results.items()}