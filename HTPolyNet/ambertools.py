# -*- coding: utf-8 -*-
"""

@author: huang
"""
import subprocess
import os
import parmed  # for converting Amber-style outputs of tleap to gromacs-format

def GAFFParameterize(inputMol2,outputPrefix,resName='UNK',force=False,extra_antechamber_params='',parmed_save_inline=True):
    message=f'GAFF is parameterizing {inputMol2}\n'
    groOut=f'{outputPrefix}.gro'
    topOut=f'{outputPrefix}.top'
    itpOut=f'{outputPrefix}.itp'
    if os.path.isfile(groOut) and os.path.isfile(topOut) and not force:
        message+=f'   {groOut} and {topOut} already exist,\n'
        message+=f'   and GAFFParameterize called with force=False.\n'
        return message
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.gaff\n')
        f.write(f'SUS = loadmol2 {inputMol2}\n')
        f.write(f'check SUS\n')
        f.write(f'loadamberparams {outputPrefix}.frcmod\n')
        f.write(f'saveamberparm SUS {outputPrefix}-tleap.top {outputPrefix}-tleap.crd\n')
        f.write('quit\n')
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
        file=parmed.load_file(f'{outputPrefix}-tleap.top', xyz=f'{outputPrefix}-tleap.crd')
        file.save(groOut)
        if parmed_save_inline:
            file.save(topOut)
        else:
            file.save(topOut,parameters=itpOut,overwrite=True)
    except Exception as m:
        raise parmed.exceptions.GromacsError(m)
    return message
        
    

    
