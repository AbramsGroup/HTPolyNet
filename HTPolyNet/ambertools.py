# -*- coding: utf-8 -*-
"""

@author: huang
"""
import subprocess
import os
import parmed
    
class Parameterization(object):
    def __init__(self):
        self.GMX = 'gmx'
        
    def GAFFParameterize(self,inputMol2,outputPrefix,resName='UNK',force=False,extra_antechamber_params='',parmed_save_inline=True):
        message=f'GAFF is parameterizing {inputMol2}\n'
        groOut=f'{outputPrefix}.gro'
        topOut=f'{outputPrefix}.top'
        if os.path.isfile(groOut) and os.path.isfile(topOut) and not force:
            message+=f'   {groOut} and {topOut} already exist.\n'
            return message
        else:
            str1 = f'source leaprc.gaff\nSUS = loadmol2 {outputPrefix}.mol2 \ncheck SUS\nloadamberparams {outputPrefix}.frcmod \nsaveamberparm SUS {outputPrefix}-gaff.top {outputPrefix}.crd\nquit'
            with open('tleap.in', 'w') as f:
                f.write(str1)
            command1 = f'antechamber -j 4 -fi mol2 -fo mol2 -c gas -at gaff -rn {resName} -i {inputMol2} -o {outputPrefix}.mol2 -pf Y -nc 0 {extra_antechamber_params}'
            command2 = f'parmchk2 -i {outputPrefix}.mol2 -o {outputPrefix}.frcmod -f mol2 -s gaff'
            command3 = 'tleap -f tleap.in'
            
            commands=[command1,command2,command3]

            for c in commands:
                process=subprocess.Popen(c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
                out,err=process.communicate()
                if process.returncode!=0:
                    message=f'command "{c}" failed with returncode {process.returncode}:\n'
                    message+=out+'\n'+err+'\n'
                    raise subprocess.SubprocessError(message)
                else:
                    message+=f'command "{c}" successful.  Output:\n'
                    message+=out+'\n'
            
            # save the results of the antechamber/parmchk2/tleap sequence as Gromacs gro and top files
            try:
                file=parmed.load_file(f'{outputPrefix}-gaff.top', xyz=f'{outputPrefix}.crd')
                file.save(groOut)
                if parmed_save_inline:
                    file.save(topOut)
                else:
                    file.save(topOut,parameters=f'{outputPrefix}.itp',overwrite=True)
            except Exception as m:
                raise parmed.exceptions.GromacsError(m)
            return message
        
    

    
