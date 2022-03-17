''' subdirectory tree for a run directory '''
import os
import sys
import pathlib
'''
root
    |-->proj0
    |-->proj1
    |...
    +-->projn
            |-->basic
            |-->mdp
            |-->systems
            |       |-->unrctSystems
            |       |-->rctSystems
            |       +-->typeSystems
            |
            +-->results
                    |-->init
                    |-->sim0
                    |-->sim1
                    |...
                    +-->simn
                            +-->init


'''

class ProjectFileSystem:
    def __init__(self,root='.',reProject=''):
        self.rootPath=root
        self.D=[]
        if reProject=='': # this is a fresh project
            i=0
            while(os.path.isdir(os.path.join(root,f'proj{i}'))):
                i+=1
            reProject=f'proj{i}'
        else:
            if not os.path.isdir(os.path.join(root,reProject)):
                raise FileNotFoundError(f'{reProject} not found')
        self.D.append(reProject)
        os.mkdir(os.path.join(self.D['root'],self.D['project']))

    def populate(self):
        toplevels=['basic','mdp','results','systems']
        for tl in toplevels:
            R[tl]=os.path.join(projPath,tl)
        syssubs=['unrctSystem','rctSystem','typeSystem']
        for ss in syssubs:
            R[ss]=os.path.join(R['systems'],ss)
        ressubs=['init']
        for rs in ressubs:
            R[rs]=os.path.join(R['results'],rs)
        self.D=R

