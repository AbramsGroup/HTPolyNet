# -*- coding: utf-8 -*-
"""

@author: huang
"""
import subprocess
import os
import parmed
from HTPolyNet.software import Command

def GAFFParameterize(inputPrefix,outputPrefix,parmed_save_inline=True,force=False,**kwargs):
    chargemodel=kwargs.get('chargemodel','gas')
    message=f'Ambertools: parameterizing {inputPrefix}\n'
    mol2in=f'{inputPrefix}.mol2'
    mol2out=f'{outputPrefix}.mol2'
    frcmodout=f'{outputPrefix}.frcmod'
    if mol2in==mol2out:
        message+=f'Warning: Antechamber will overwrite {mol2in}\n'
    groOut=f'{outputPrefix}.gro'
    topOut=f'{outputPrefix}.top'
    itpOut=f'{outputPrefix}.itp'
    if os.path.isfile(groOut) and os.path.isfile(topOut) and not force:
        message+=f'   {groOut} and {topOut} already exist,\n'
        message+=f'   and GAFFParameterize called with force=False.\n'
        return message
    c=Command('antechamber',j=4,fi='mol2',fo='mol2',c=chargemodel,at='gaff',i=mol2in,o=mol2out,pf='Y',nc=0,eq=1,pl=10)
    message+=c.run()
    c=Command('parmchk2',i=mol2out,o=frcmodout,f='mol2',s='gaff')
    message+=c.run()
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.gaff\n')
        f.write(f'SUS = loadmol2 {mol2out}\n')
        f.write(f'check SUS\n')
        f.write(f'loadamberparams {frcmodout}\n')
        f.write(f'saveamberparm SUS {outputPrefix}-tleap.top {outputPrefix}-tleap.crd\n')
        f.write('quit\n')
    c=Command('tleap',f='tleap.in')
    message+=c.run()
    # save the results of the antechamber/parmchk2/tleap sequence as Gromacs gro and top files
    try:
        file=parmed.load_file(f'{outputPrefix}-tleap.top', xyz=f'{outputPrefix}-tleap.crd')
        message+=f'Writing {groOut}\n'
        file.save(groOut)
        if parmed_save_inline:
            message+=f'Writing {topOut}\n'
            file.save(topOut)
        else:
            message+=f'Writing {topOut} and {itpOut}\n'
            file.save(topOut,parameters=itpOut,overwrite=True)
    except Exception as m:
        raise parmed.exceptions.GromacsError(m)
    return message
        
    

    
