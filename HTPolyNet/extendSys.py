# THIS FILE IS NO LONGER USED
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:03:09 2020

@author: huang, abrams
"""
import os
from HTPolyNet.gromacs import GMXCommand

def insert_molecules(molInfo,boxSize,outName):
    if type(boxSize)==float:
        boxSize=[boxSize]*3
    assert len(boxSize)==3, f'Error: malformed boxsize {boxSize}'
    message=''
    for mol in molInfo:
        name = mol[1]
        num = mol[2]
        if os.path.isfile(f'{outName}.gro'):
            ''' final gro file exists; we must insert into it '''
            c=GMXCommand('insert-molecules',f=f'{outName}.gro',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=0.4)
        else:
            ''' no final gro file yet; make it '''
            c=GMXCommand('insert-molecules',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=0.4)
        message+=c.run()
    return message
                
def extendSys(monInfo,croInfo,boxSize,fileName):
    msg=''
    msg+=insert_molecules(monInfo,boxSize,fileName)
    msg+=insert_molecules(croInfo,boxSize,fileName)
    return msg
        
