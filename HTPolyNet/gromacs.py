'''
gromacs.py -- simple class for handling gromacs commands
'''
import subprocess
import os

class GMXCommand:
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
        c=f'gmx {self.command} '+' '.join([f'-{k} {v}' for k,v in self.options.items()])
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
    assert len(boxSize)==3, f'Error: malformed boxsize {boxSize}'
    scale=kwargs.get('scale',0.4) # our default vdw radius scaling
    message=''
    for n,m in monomers.items():
        name = n+kwargs.get('basename_modifier','')
        num = composition[n]
        if os.path.isfile(f'{outName}.gro'):
            ''' final gro file exists; we must insert into it '''
            c=GMXCommand('insert-molecules',f=f'{outName}.gro',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=scale)
        else:
            ''' no final gro file yet; make it '''
            c=GMXCommand('insert-molecules',ci=f'{name}.gro',nmol=num,o=outName,box=' '.join([f'{x:.8f}' for x in boxSize]),scale=scale)
        message+=c.run()
    return message

# def extendSys(monomers,composition,boxSize,fileName):
#     msg=''
#     msg+=insert_molecules(gro,num,boxSize,fileName)
#     # msg+=insert_molecules(croInfo,boxSize,fileName)
#     return msg

def grompp_and_mdrun(gro='',top='',out='',mdp='',boxSize=[],**kwargs):
    ''' launcher for grompp and mdrun
        gro: prefix for input coordinate file
        top: prefix for input topology
        out: prefix for desired output files
        boxsize: (optional) desired box size; triggers editconf before grompp
    '''
    if gro=='' or top=='' or out=='' or mdp=='':
        raise Exception('grompp_and_run requires gro, top, out, and mdp filename prefixes.')
    msg=''
    if len(boxSize)>0:
        c=GMXCommand('editconf',f=f'{gro}.gro',o=gro,
                     box=' '.join([f'{x:.8f}' for x in boxSize]))
        msg+=c.run()
    maxwarn=kwargs.get('maxwarn',2)
    c=GMXCommand('grompp',f=f'{mdp}.mdp',c=f'{gro}.gro',p=f'{top}.top',o=f'{out}.tpr',maxwarn=maxwarn)
    msg+=c.run()
    c=GMXCommand('mdrun',deffnm=out)
    msg+=c.run()
    assert os.path.exists(f'{out}.gro'), 'Error: mdrun failed.'
    return msg