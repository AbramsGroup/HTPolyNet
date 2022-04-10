''' Check for presence of required software '''
import logging
import subprocess

class Command:
    linelen=55
    def __init__(self,command,**options):
        self.command=command
        self.options=options
        self.c=f'{self.command} '+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        
    def run(self,override=()):
        logging.info(f'Issuing command "{self.c}"...')
        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        out,err=process.communicate()
        if process.returncode!=0:
            logging.error(f'Command "{self.c}" exited with returncode {process.returncode}')
            logging.error('stdout buffer follows\n'+'*'*self.linelen+'\n'+out+'\n'+'*'*self.linelen)
            logging.error('stderr buffer follows\n'+'*'*self.linelen+'\n'+err+'\n'+'*'*self.linelen)
            raise subprocess.SubprocessError(f'Command "{self.c}" failed with returncode {process.returncode}')
        else:
            logging.info(f'Command "{self.c}" exited with returncode {process.returncode}.')
            if len(override)==2:
                needle,msg=override
                if needle in out or needle in err:
                    logging.error('*'*self.linelen+'\n'+msg+'*'*self.linelen)
        return out,err

class Software:
    ambertools=['antechamber','tleap','parmchk2']
    commandsRequired=ambertools+['gmx']
    optionalCommands=['gmx_mpi','mdrun_mpi']
    def __init__(self):
        self.commands={}
        cnf=[]
        passes=True
        for c in Software.commandsRequired:
            CP=subprocess.run(['which',c],capture_output=True,text=True)
            if CP.returncode!=0:
                passes=False
                cnf.append(c)
            else:
                self.commands[c]=CP.stdout.strip()
        if not passes:
            logging.error(f'HTPolyNet cannot find command(s) {", ".join(cnf)}')
            raise FileNotFoundError(f'HTPolyNet cannot find command(s) {", ".join(cnf)}')
        cnf=[]
        passes=True
        for c in Software.optionalCommands:
            CP=subprocess.run(['which',c],capture_output=True,text=True)
            if CP.returncode!=0:
                passes=False
                cnf.append(c)
            else:
                self.commands[c]=CP.stdout.strip()
        if not passes:
            pass
        self.getVersions()

    def __str__(self):
        r='Commands available for HTPolyNet to use:\n'
        for c,p in self.commands.items():
            if not c in self.ambertools:
                verkey=c
            else:
                verkey='ambertools'
            r+=f'{c:>12s} (ver. {self.versions[verkey]:>6s}) at {p:<50s}\n'
        return r

    def getVersions(self):
        self.versions={}
        CP=subprocess.run(['gmx','-version'],capture_output=True,text=True)
        self.versions['gmx']=CP.stdout.split('\n')[0].split()[4].strip()
        CP=subprocess.run(['obabel','-V'],capture_output=True,text=True)
        self.versions['obabel']=CP.stdout.split()[2].strip()
        CP=subprocess.run(['antechamber','-h'],capture_output=True,text=True)
        l=CP.stdout.split('\n')[1].split()[3].strip().strip(':')
        self.versions['ambertools']=l
        if 'mdrun_mpi' in self.commands:
            CP=subprocess.run(['mdrun_mpi','-h'],capture_output=True,text=True)
            self.versions['mdrun_mpi']=CP.stdout.split('\n')[0].split()[4].strip()
        
    def info(self):
        print(str(self))


