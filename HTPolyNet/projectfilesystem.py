''' Classes for handling resource library and a runtime file system '''

import logging
import os
import pathlib
import importlib.resources

_MolecularDataFileTypes_=['mol2','gro','top','itp','sea']
_MoleculeClasses_=['monomers','oligomers']

class Library:
    ''' a library object -- default creation references the Library resource package. '''
    def __init__(self,libpackage='Library'):
        self.fullpaths={}
        try:
            # Important! This will allow the HTPolyNet runtime to know where the 
            # Library package lives on the system without the user needing to know
            tt=importlib.resources.files(libpackage)
        except:
            raise ImportError(f'Could not find library package {libpackage}.  Your HTPolyNet installation is corrupt.')
        for n in tt.iterdir():
            if os.path.isdir(n) and '__pycache__' not in str(n):
                bn=os.path.basename(n)
                self.fullpaths[bn]=n

    def dir(self,filename):
        prefix,ext=os.path.splitext(filename)
        try:
            ldir=os.path.join(self.fullpaths[ext],filename)
        except KeyError as msg:
            raise Exception(f'No top-level resource directory "{ext}"')
        if ext in _MolecularDataFileTypes_:
            for sd in _MoleculeClasses_:
                testfile=os.path.join(ldir,sd,filename)
                if os.path.exists(testfile):
                    return os.path.join(ldir,sd)


    def dir(self,typestr):
        return self.libdir.get(typestr,None)

    def __str__(self):
        s=''
        for k,v in self.libdir.items():
            s+=f'{k}: {v}\n'
        return s

    def info(self):
        print('\nHTPolyNet libraries:')
        for k,v in self.libdir.items():
            f=list(os.listdir(v))
            l=len(f)
            ess='' if l==1 else 's'
            c=':' if l>0 else '.'
            print(f'  {k} ({v}) has {len(f)} file{ess}{c}')
            if len(f)>0:
                maxlen=max([len(i) for i in f])
                ncol=70//(maxlen+2)
                outstr=r'  {:>'+str(maxlen)+r's}'
                for i in range(0,len(f),ncol):
                    print('    ',end='')
                    for j in range(i,i+ncol):
                        if j<len(f):
                            print(outstr.format(f[j]),end='')
                    print()


class ProjectFileSystem:
    def __init__(self,root='.',verbose=False,reProject=False):
        self.rootPath=pathlib.Path(root).resolve()
        self.cwd=self.rootPath
        self.verbose=verbose
        self.cd(self.rootPath)
        self._next_project_dir(reProject=reProject)
        self._setup_project_root()
        self.cd(self.rootPath)
        self.library=Library()

    def cd(self,dest=''):
        os.chdir(dest)
        self.cwd=os.getcwd()
        if (self.verbose):
            logging.info(f'cwd: {self.cwd}')

    def cdroot(self):
        self.cd(self.projPath)

    def cdrootsub(self,toplevel):
        if toplevel in self.projSubPaths:
            self.cd(self.projSubPaths[toplevel])

    def exist(self,name):
        ext=name.split('.')[-1]
        fullname=os.path.join(self.library.dir(ext),name)
        return os.path.exists(fullname)

    def fetch(self,name,destpath='.',overwrite='yes',altpath=None):
        ext=name.split('.')[-1]
        searchpath=[self.rootPath,self.projPath,self.library.dir(ext)]
        if altpath!=None:
            searchpath=[altpath]+searchpath
        destfname=os.path.join(destpath,name)
        for p in searchpath:
            afname=os.path.join(p,name)
            if os.path.exists(afname):
                if os.path.exists(destfname):
                    if overwrite=='yes':
                        logging.info(f'Fetching {afname} into {destpath} (overwriting {destfname})')
                        os.system(f'cp {afname} {destpath}')
                    elif overwrite=='if_newer':
                        logging.info(f'Fecthing newer {afname} into {destpath} (overwriting {destfname})')
                        if os.path.getctime(destfname) < os.path.getctime(afname):
                            os.system(f'cp {afname} {destpath}')
                    else:
                        logging.info(f'Refusing to fetch older {afname} into {destpath} because existing {destfname} is newer')
                else:
                    logging.info(f'Fetching {afname} into {destpath}')
                    os.system(f'cp {afname} {destpath}')
                return True
        else:
            return False

    def store(self,name,altdestpath=None,overwrite='yes'):
        ''' stores files in package library or in an alternate location'''
        ext=name.split('.')[-1]
        libdestpath=self.library.dir(ext)
        destfullname=os.path.join(libdestpath if not altdestpath else altdestpath,name)
        if os.path.exists(destfullname):
            if overwrite=='yes':
                logging.info(f'Storing {name} over {destfullname}')
                os.system(f'cp {name} {destfullname}')
            elif overwrite=='if_newer':
                logging.info(f'Storing newer {name} over {destfullname}')
                if os.path.getctime(destfullname) < os.path.getctime(name):
                    os.system(f'cp {name} {destfullname}')
            else:
                logging.info(f'Refusing to store older {name} over {destfullname}')
        else:
            logging.info(f'Storing {name} to {destfullname}')
            os.system(f'cp {name} {destfullname}')

    def __str__(self):
        return f'root {self.rootPath}: cwd {self.cwd}'

    def _next_project_dir(self,reProject=False):
        i=0
        lastprojdir=''
        currentprojdir=''
        while(os.path.isdir(os.path.join(self.rootPath,f'proj{i}'))):
            lastprojdir=f'proj{i}'
            i+=1
        if not reProject or lastprojdir=='': # this is a fresh project
            if lastprojdir=='':
                currentprojdir='proj0'
            else:
                currentprojdir=f'proj{i}'
            self.projPath=os.path.join(self.rootPath,currentprojdir)
            os.mkdir(currentprojdir)
        else:
            self.projPath=os.path.join(self.rootPath,lastprojdir)

    def _setup_project_root(self):
        self.cdroot()
        self.projSubPaths={}
        for tops in ['basic','mdp','systems','results']:
            self.projSubPaths[tops]=os.path.join(self.projPath,tops)
            if not os.path.isdir(self.projSubPaths[tops]):
                os.mkdir(tops)
        self.basicPath=self.projSubPaths['basic']
        self.mdpPath=self.projSubPaths['mdp']
        self.systemsPath=self.projSubPaths['systems']
        self.resPath=self.projSubPaths['results']
        self.cd(self.systemsPath)
        self.systemsSubPaths={}
        self.resultsSubPaths={}
        for tops in ['rctSystems','unrctSystems','typeSystems']:
            self.systemsSubPaths[tops]=os.path.join(self.systemsPath,tops)
            if not os.path.isdir(self.systemsSubPaths[tops]):
                os.mkdir(tops)
        self.unrctPath=self.systemsSubPaths['unrctSystems']
        self.rctPath=self.systemsSubPaths['rctSystems']
        self.typePath=self.systemsSubPaths['typeSystems']

    def next_results_dir(self):
        possibles=['init',*[f'step{i}' for i in range(30)]]
        for p in possibles:
            if not p in self.resultsSubPaths:
                break
        newpath=os.path.join(self.resPath,p)
        os.mkdir(newpath)
        self.resultsSubPaths[p]=newpath
        return newpath
        
if __name__=='__main__':
    pfs=ProjectFileSystem(verbose=True)

