# -*- coding: utf-8 -*-
"""
Created on Sun Oct 18 15:03:09 2020

@author: huang
"""
import subprocess
import os

class extendSys(object):
    def __init__(self, gmx):
        self.GMX = gmx
    
    def insertMol(self, molInfo, boxSize, outName):
        GMX = self.GMX
        if isinstance(boxSize) == list:
            if len(boxSize) == 3:
                pass
            else:
                raise ValueError(f'wrong boxSize: {boxSize}')
        else:
            boxSize = [boxSize, boxSize, boxSize]

        for mol in molInfo:
            name = mol[1]
            num = mol[2]
            if os.path.isfile('{}.gro'.format(outName)):
                cmd1 = '{} insert-molecules -f {}.gro -ci {}.gro -nmol {} -o {} -box {} {} {} -scale 0.4'.format(
                    GMX, outName, name, num, outName, boxSize[0], boxSize[1], boxSize[2])
                a1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = a1.communicate()
            else:
                cmd1 = '{} insert-molecules -ci {}.gro -nmol {} -o {} -box {} {} {} -scale 0.4'.format(
                    GMX, name, num, outName, boxSize[0], boxSize[1], boxSize[2])
                a1 = subprocess.Popen(cmd1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                out, err = a1.communicate()
                
    def extendSys(self, monInfo, croInfo, boxSize, fileName):
        print('Box size: ', boxSize)
        self.insertMol(monInfo, boxSize, fileName)
        self.insertMol(croInfo, boxSize, fileName)
        