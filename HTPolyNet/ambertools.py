# -*- coding: utf-8 -*-
"""

@author: huang
"""
import hashlib
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
    # tleap can't handle symbols in file names
    leapprefix=hashlib.shake_128(outputPrefix.encode("utf-8")).hexdigest(8)
    message+='Replacing string "{mol2out}" with hash "{leapprefix}" for leap input files.\n'
    Command(f'cp {mol2out} {leapprefix}.mol2').run()
    Command(f'cp {frcmodout} {leapprefix}.frcmod').run()
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.gaff\n')
        f.write(f'SUS = loadmol2 {leapprefix}.mol2\n')
        f.write(f'check SUS\n')
        f.write(f'loadamberparams {leapprefix}.frcmod\n')
        f.write(f'saveamberparm SUS {leapprefix}-tleap.top {leapprefix}-tleap.crd\n')
        f.write('quit\n')
    c=Command('tleap',f='tleap.in')
    message+=c.run()
    Command(f'mv {leapprefix}-tleap.top {outputPrefix}-tleap.top').run()
    Command(f'mv {leapprefix}-tleap.crd {outputPrefix}-tleap.crd').run()
    # save the results of the antechamber/parmchk2/tleap sequence as Gromacs gro and top files
    try:
        file=parmed.load_file(f'{outputPrefix}-tleap.top', xyz=f'{outputPrefix}-tleap.crd')
        message+=f'Writing {groOut}\n'
        file.save(groOut,overwrite=True)
        if parmed_save_inline:
            message+=f'Writing {topOut}\n'
            file.save(topOut,overwrite=True)
        else:
            message+=f'Writing {topOut} and {itpOut}\n'
            file.save(topOut,parameters=itpOut,overwrite=True)
    except Exception as m:
        raise parmed.exceptions.GromacsError(m)
    return message
        
    

    
