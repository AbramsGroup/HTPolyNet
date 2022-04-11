'''
gromacs.py -- simple class for handling gromacs commands
'''
import logging
import os
import pandas as pd
from HTPolyNet.software import Command

def insert_molecules(monomers,composition,boxSize,outName,**kwargs):
    ''' launcher for `gmx insert-molecules`
        monomers:  dictionary of Monomer instances
        composition: dictionary keyed on monomer name with value monomer count
        boxSize:  3-element list of floats OR a single float (cubic box)
        outName:  output filename basename.  If {outName}.gro exists,
                  insertions are made into it.
    '''
    if type(boxSize)==int:
        boxSize=float(boxSize)
    if type(boxSize)==float:
        boxSize=[boxSize]*3
    scale=kwargs.get('scale',0.4) # our default vdw radius scaling
    for n,m in monomers.items():
        name = n+kwargs.get('basename_modifier','')
        num = composition[n]
        if os.path.isfile(f'{outName}.gro'):
            logging.info(f'gmx insert-molecules inserts into existing {outName}.gro')
            ''' final gro file exists; we must insert into it '''
            c=Command('gmx insert-molecules',f=f'{outName}.gro',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=scale)
        else:
            ''' no final gro file yet; make it '''
            c=Command('gmx insert-molecules',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=scale)
        out,err=c.run()
        out+=err
        logging.info(f'Output of "gmx insert-molecules"\n'+out)
        if 'Added' in out:
            outlines=out.split('\n')
            for ol in outlines:
                if 'Added' in ol:
                    tokens=ol.split()
                    numadded=int(tokens[1])
                    break
            if numadded!=num:
                logging.error(f'gmx insert molecules did not add enough {name}; only {numadded} out of {num} were placed.  Increase your boxsize.')
                raise Exception('need bigger box')

def grompp_and_mdrun(gro='',top='',out='',mdp='',boxSize=[],**kwargs):
    ''' launcher for grompp and mdrun
        gro: prefix for input coordinate file
        top: prefix for input topology
        out: prefix for desired output files
        boxsize: (optional) desired box size; triggers editconf before grompp
    '''
    if gro=='' or top=='' or out=='' or mdp=='':
        raise Exception('grompp_and_run requires gro, top, out, and mdp filename prefixes.')
    if len(boxSize)>0:
        c=Command('gmx editconf',f=f'{gro}.gro',o=gro,
                     box=' '.join([f'{x:.8f}' for x in boxSize]))
        c.run()
    maxwarn=kwargs.get('maxwarn',2)
    c=Command('gmx grompp',f=f'{mdp}.mdp',c=f'{gro}.gro',p=f'{top}.top',o=f'{out}.tpr',maxwarn=maxwarn)
    c.run()
    c=Command('gmx mdrun',deffnm=out)
    c.run()
    if os.path.exists(f'{out}.gro'):
        logging.info(f'grompp_and_run completed.  Check {gro}.gro.')
    else:
        logging.error(f'gmx mdrun ended prematurely; {gro}.gro not found.')
        raise Exception(f'gmx mdrun ended prematurely; {gro}.gro not found.')

def density_trace(edr='',**kwargs):
    if edr=='':
        raise Exception('density_trace requires and edr filename prefix.')
    with open('gmx.in','w') as f:
        f.write('22\n\n')
    c=Command(f'gmx energy -f {edr}.edr -o {edr}-density.dat -xvg none < gmx.in')
    c.run()
    density=pd.read_csv(f'{edr}-density.dat',sep='\s+',names=['time(ps)','density(kg/m^3)'])
    density['ra']=density['density(kg/m^3)'].expanding(1).mean()
    logging.info(density.to_string())