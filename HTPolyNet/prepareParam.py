# -*- coding: utf-8 -*-
"""

@author: huang
"""
import subprocess
import os
import parmed
    
class prepareParam(object):
    def __init__(self):
        self.GMX = 'gmx'
        
    def PrepareFile(self, fileName, resName, outputName):
        print('fileName: ', fileName)
        if os.path.isfile('{}.gro'.format(outputName)) and os.path.isfile('{}.top'.format(outputName)):
            pass
        else:
            command1 = 'antechamber -fi mol2 -fo mol2 -c gas -at gaff -rn {} -i {} -o out.mol2 -pf Y -nc 0'.format(resName, fileName)
            command2 = 'parmchk2 -i out.mol2 -o out.frcmod -f mol2 -s gaff'
            command3 = 'tleap -f tleap.in'
            
            str1 = 'source leaprc.gaff\nSUS = loadmol2 out.mol2 \ncheck SUS\nloadamberparams out.frcmod \nsaveamberparm SUS out.top out.crd \nquit'
            with open('tleap.in', 'w') as f:
                f.write(str1)
            
            a1 = subprocess.Popen(command1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out1, err1 = a1.communicate()
            a2 = subprocess.Popen(command2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out2, err2 = a2.communicate()
            a3 = subprocess.Popen(command3, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out3, err3 = a3.communicate()
            
            file = parmed.load_file('out.top', xyz='out.crd')
            file.save('{}.gro'.format(outputName))
            file.save('{}.top'.format(outputName))
        
    

    
