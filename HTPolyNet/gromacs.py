'''
gromacs.py -- simple class for handling gromacs commands
'''
import subprocess

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
