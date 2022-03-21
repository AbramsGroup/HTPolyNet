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
                            |-->init


'''
class ProjectFileSystem:
    def __init__(self,root='.',verbose=False,reProject=''):
        self.rootPath=pathlib.Path(root).resolve()
        self.cwd=self.rootPath
        self.verbose=verbose
        self.cd(self.rootPath)
        self._nextProjectDir(reProject='')
        self._setupProjectRoot()
        self.cd(self.rootPath)

    def cd(self,dest=''):
        os.chdir(dest)
        self.cwd=os.getcwd()
        if (self.verbose):
            print(f'cwd: {self.cwd}')

    def goToProjectRoot(self):
        self.cd(self.projPath)

    def goToProjectSubPath(self,toplevel):
        if toplevel in self.projSubPaths:
            self.cd(self.projSubPaths[toplevel])

    def __str__(self):
        return f'root {self.rootPath}: {self.D} cwd {self.cwd}'

    def _nextProjectDir(self,reProject=''):
        if reProject=='': # this is a fresh project
            i=0
            while(os.path.isdir(os.path.join(self.rootPath,f'proj{i}'))):
                i+=1
            reProject=f'proj{i}'
        else:
            if not os.path.isdir(os.path.join(self.rootPath,reProject)):
                raise FileNotFoundError(f'{reProject} not found')
        self.projPath=os.path.join(self.rootPath,reProject)
        os.mkdir(reProject)


    def _setupProjectRoot(self):
        self.goToProjectRoot()
        self.projSubPaths={}
        for tops in ['basic','mdp','systems','results']:
            self.projSubPaths[tops]=os.path.join(self.projPath,tops)
            os.mkdir(tops)
        self.basicPath=self.projSubPaths['basic']
        self.mdpPath=self.projSubPaths['mdp']
        self.systemsPath=self.projSubPaths['systems']
        self.resPath=self.projSubPaths['results']
        self.cd(self.systemsPath)
        self.systemsSubPaths={}
        for tops in ['rctSystems','unrctSystems','typeSystems']:
            self.systemsSubPaths[tops]=os.path.join(self.systemsPath,tops)
            os.mkdir(tops)
        self.unrctPath=self.systemsSubPaths['unrctSystems']
        self.rctPath=self.systemsSubPaths['rctSystems']
        self.typePath=self.systemsSubPaths['typeSystems']


if __name__=='__main__':
    pfs=ProjectFileSystem(verbose=True)

