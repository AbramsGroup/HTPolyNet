"""

.. module:: command
   :synopsis: Custom handling of calls to subprocess.Popen() 
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import subprocess
logger=logging.getLogger(__name__)

class Command:
    linelen=55
    def __init__(self,command,**options):
        self.command=command
        self.options=options
        self.c=f'{self.command} '+' '.join([f'-{k} {v}' for k,v in self.options.items()])
        
    def run(self,override=(),ignore_codes=[],quiet=True):
        if not quiet:
            logger.debug(f'{self.c}')
        process=subprocess.Popen(self.c,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
        out,err=process.communicate()
        if process.returncode!=0 and not process.returncode in ignore_codes:
            logger.error(f'Returncode: {process.returncode}')
            if len(out)>0:
                logger.error('stdout buffer follows\n'+'*'*self.linelen+'\n'+out+'\n'+'*'*self.linelen)
            if len(err)>0:
                logger.error('stderr buffer follows\n'+'*'*self.linelen+'\n'+err+'\n'+'*'*self.linelen)
            raise subprocess.SubprocessError(f'Command "{self.c}" failed with returncode {process.returncode}')
        else:
            if len(override)==2:
                needle,msg=override
                if needle in out or needle in err:
                    logger.error(f'Returncode: {process.returncode}, but another error was detected:')
                    logger.error(msg)
                    if len(out)>0:
                        logger.error('stdout buffer follows\n'+'*'*self.linelen+'\n'+out+'\n'+'*'*self.linelen)
                    if len(err)>0:
                        logger.error('stderr buffer follows\n'+'*'*self.linelen+'\n'+err+'\n'+'*'*self.linelen)
        return out,err