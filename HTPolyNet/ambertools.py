# -*- coding: utf-8 -*-
"""

@author: huang
"""
import logging
import hashlib
import shutil
import parmed
import os
from HTPolyNet.command import Command
from HTPolyNet.coordinates import Coordinates
logger=logging.getLogger(__name__)

def GAFFParameterize(inputPrefix,outputPrefix,input_structure_format='mol2',**kwargs):
    chargemethod=kwargs.get('charge_method','bcc')
    logger.info(f'AmberTools is parameterizing {inputPrefix}.{input_structure_format}')
    structin=f'{inputPrefix}.{input_structure_format}'
    mol2out=f'{outputPrefix}.mol2'
    frcmodout=f'{outputPrefix}.frcmod'
    new_structin=structin
    if structin==mol2out:
        # mol2pfx,mol2sfx=os.path.splitext(mol2in)
        new_structin=inputPrefix+f'-input.{input_structure_format}'
        logger.debug(f'Antechamber overwrites input {structin}; backing up to {new_structin}')
        shutil.copy(structin,new_structin)
    groOut=f'{outputPrefix}.gro'
    topOut=f'{outputPrefix}.top'
    itpOut=f'{outputPrefix}.itp'
    c=Command('antechamber',j=4,fi=input_structure_format,fo='mol2',c=chargemethod,at='gaff',i=new_structin,o=mol2out,pf='Y',nc=0,eq=1,pl=10)
    c.run(quiet=False)
    logger.debug(f'Antechamber generated {mol2out}')
    c=Command('parmchk2',i=mol2out,o=frcmodout,f='mol2',s='gaff')
    c.run(quiet=False)
    # Antechamber ignores SUBSTRUCTURES but we would like tleap to 
    # recognize them.  So we will simply use the antechamber-INPUT mol2 file
    # resName and resNum atom record fields over the antechamber-OUTPUT mol2 file
    # to generate the tleap-INPUT mol2 file.  Also, since tleap can't handle
    # wacky filenames, we'll hash the outputPrefix temporarily.
    if input_structure_format=='mol2':
        goodMol2=Coordinates.read_mol2(new_structin)
        acOutMol2=Coordinates.read_mol2(mol2out)
        adf=acOutMol2.A
        adf['resName']=goodMol2.A['resName'].copy()
        adf['resNum']=goodMol2.A['resNum'].copy()
        goodMol2.A=adf
    else:
        goodMol2=Coordinates.read_mol2(mol2out)
    leapprefix=hashlib.shake_128(outputPrefix.encode("utf-8")).hexdigest(16).replace('e','x')
    goodMol2.write_mol2(f'{leapprefix}.mol2')
    logger.debug(f'Replacing string "{outputPrefix}" with hash "{leapprefix}" for leap input files.')
    Command(f'cp {frcmodout} {leapprefix}.frcmod').run()
    with open(f'{inputPrefix}-tleap.in', 'w') as f:
        f.write(f'source leaprc.gaff\n')
        f.write(f'mymol = loadmol2 {leapprefix}.mol2\n')
        f.write(f'check mymol\n')
        f.write(f'loadamberparams {leapprefix}.frcmod\n')
        f.write(f'saveamberparm mymol {leapprefix}-tleap.top {leapprefix}-tleap.crd\n')
        f.write('quit\n')
    c=Command('tleap',f=f'{inputPrefix}-tleap.in')
    c.run(override=('Error!','Unspecified tleap error'))
    Command(f'rm -f {leapprefix}.frcmod').run()
    Command(f'mv {leapprefix}.mol2 {outputPrefix}.mol2').run()
    Command(f'mv {leapprefix}-tleap.top {outputPrefix}-tleap.top').run()
    Command(f'mv {leapprefix}-tleap.crd {outputPrefix}-tleap.crd').run()
    # save the results of the antechamber/parmchk2/tleap sequence as Gromacs gro and top files
    try:
        file=parmed.load_file(f'{outputPrefix}-tleap.top', xyz=f'{outputPrefix}-tleap.crd')
        logger.debug(f'Writing {groOut}, {topOut}, and {itpOut}')
        file.save(groOut,overwrite=True)
        file.save(topOut,parameters=itpOut,overwrite=True)
    except Exception as m:
        logger.error('Unspecified parmed error')
        raise parmed.exceptions.GromacsError(m)
        
    

    
