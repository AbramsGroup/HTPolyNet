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
        for mol in molInfo:
            name = mol[1]
            num = mol[2]
            print('cwd: ', os.getcwd())
            print('extend outName: ', outName)
            if os.path.isfile('{}.gro'.format(outName)):
                cmd1 = '{} insert-molecules -f {}.gro -ci {}.gro -nmol {} -o {} -box {} {} {} -scale 0.4'.format(
                    GMX, outName, name, num, outName, boxSize, boxSize, boxSize)
                a1 = subprocess.call(cmd1, shell=True) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                #a1.wait()
            else:
                cmd1 = '{} insert-molecules -ci {}.gro -nmol {} -o {} -box {} {} {} -scale 0.4'.format(
                    GMX, name, num, outName, boxSize, boxSize, boxSize)
                a1 = subprocess.call(cmd1, shell=True) #, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                
    def extendSys(self, monInfo, croInfo, boxSize, fileName):
        print('Box size: ', boxSize)
        self.insertMol(monInfo, boxSize, fileName)
        self.insertMol(croInfo, boxSize, fileName)
        