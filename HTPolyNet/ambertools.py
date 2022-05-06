# -*- coding: utf-8 -*-
"""

@author: huang
"""
import logging
import hashlib
import shutil
import parmed
from HTPolyNet.command import Command
from HTPolyNet.coordinates import Coordinates

def GAFFParameterize(inputPrefix,outputPrefix,**kwargs):
    chargemethod=kwargs.get('chargemethod','bcc')
    logging.info(f'Ambertools: parameterizing {inputPrefix}')
    mol2in=f'{inputPrefix}.mol2'
    mol2out=f'{outputPrefix}.mol2'
    frcmodout=f'{outputPrefix}.frcmod'
    if mol2in==mol2out:
        logging.info(f'Antechamber overwrites input {mol2in}; backing up to {mol2in}.bak')
        shutil.copy(f'{mol2in}',f'{mol2in}.bak')
    groOut=f'{outputPrefix}.gro'
    topOut=f'{outputPrefix}.top'
    itpOut=f'{outputPrefix}.itp'
    c=Command('antechamber',j=4,fi='mol2',fo='mol2',c=chargemethod,at='gaff',i=mol2in,o=mol2out,pf='Y',nc=0,eq=1,pl=10)
    c.run()
    logging.info(f'Antechamber generated {mol2out}')
    c=Command('parmchk2',i=mol2out,o=frcmodout,f='mol2',s='gaff')
    c.run()
    # Antechamber ignores SUBSTRUCTURES but we would like tleap to 
    # recognize them.  So we will simply use the antechamber-INPUT mol2 file
    # resName and resNum atom record fields over the antechamber-OUTPUT mol2 file
    # to generate the tleap-INPUT mol2 file.  Also, since tleap can't handle
    # wacky filenames, we'll hash the outputPrefix temporarily.
    goodMol2=Coordinates.read_mol2(mol2in)
    acOutMol2=Coordinates.read_mol2(mol2out)
    adf=acOutMol2.A
    adf['resName']=goodMol2.A['resName'].copy()
    adf['resNum']=goodMol2.A['resNum'].copy()
    goodMol2.A=adf
    leapprefix=hashlib.shake_128(outputPrefix.encode("utf-8")).hexdigest(8)
    goodMol2.write_mol2(f'{leapprefix}.mol2')
    logging.info(f'Replacing string "{outputPrefix}" with hash "{leapprefix}" for leap input files.')
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
        logging.info(f'Writing {groOut}')
        file.save(groOut,overwrite=True)
        logging.info(f'Writing {topOut} and {itpOut}')
        file.save(topOut,parameters=itpOut,overwrite=True)
    except Exception as m:
        logging.error('Unspecified parmed error')
        raise parmed.exceptions.GromacsError(m)
        
    

    
