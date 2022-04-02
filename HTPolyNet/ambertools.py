# -*- coding: utf-8 -*-
"""

@author: huang
"""
import subprocess
import os
import parmed

class Command:
    def __init__(self,command,log=None,**options):
        self.command=command
        self.options=options
        self._openlog(log)
    def _openlog(self,filename):
        self.logio=None
        if filename:
            self.logio=open(filename,'w')
    def log(self,msg):
        if self.logio:
            self.logio.write(msg)
            self.logio.flush()
    def _closelog(self):
        if self.logio:
            self.logio.close()
    def run(self):
        c=f'{self.command} '+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        message=f'Issuing command "{c}"...\n'
        process=subprocess.Popen(c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        out,err=process.communicate()
        if process.returncode!=0:
            message=f'Command "{c}" failed with returncode {process.returncode}:\n'
            message+=out+'\n'+err+'\n'
            self.log(message)
            raise subprocess.SubprocessError(message)
        else:
            message+=f'Command "{c}" successful.\n'
            message+=out+'\n'
            self.log(message)
        self._closelog()
        return message

def GAFFParameterize(inputMol2,outputPrefix,parmed_save_inline=True,force=False):
    message=f'Ambertools: parameterizing {inputMol2}\n'
    groOut=f'{outputPrefix}.gro'
    topOut=f'{outputPrefix}.top'
    itpOut=f'{outputPrefix}.itp'
    if os.path.isfile(groOut) and os.path.isfile(topOut) and not force:
        message+=f'   {groOut} and {topOut} already exist,\n'
        message+=f'   and GAFFParameterize called with force=False.\n'
        return message
    user_Mol2=read_mol2(inputMol2)
    c=Command('antechamber',j=4,fi='mol2',fo='mol2',c='bcc',at='gaff',i=inputMol2,o=f'{outputPrefix}.mol2',pf='Y',nc=0,eq=1,pl=10)
    message+=c.run()
    ante_Mol2=read_mol2(f'{outputPrefix}.mol2')
    
    # TODO: antechamber will give atoms unique names that will persist, and they might not match what is in the input mol2 files
    # unfortunately, atom names in the input mol2 files are how users refer to particular atoms for describing reactions
    # and capping.  So, we need to map the user-specified names to current names.
    c=Command('parmchk2',i=f'{outputPrefix}.mol2',o=f'{outputPrefix}.frcmod',f='mol2',s='gaff')
    message+=c.run()
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.gaff\n')
        f.write(f'SUS = loadmol2 {inputMol2}\n')
        f.write(f'check SUS\n')
        f.write(f'loadamberparams {outputPrefix}.frcmod\n')
        f.write(f'saveamberparm SUS {outputPrefix}-tleap.top {outputPrefix}-tleap.crd\n')
        f.write('quit\n')
    c=Command('tleap',f='tleap.in')
    message+=c.run()
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
        
    

    
