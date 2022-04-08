''' Classes for handling a runtime file system '''
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
import os
import pathlib
import importlib.resources

_DefaultResourceTypes_=['cfg','mdp','mol2','top','itp','gro']
_DefaultResourcePackageDir_='Library'

class Library:
    def __init__(self,ResourceTypes=_DefaultResourceTypes_,basepath=_DefaultResourcePackageDir_):
        self.types=ResourceTypes
        self.basepath=basepath
        self.libdir={}
        if basepath==_DefaultResourcePackageDir_:
            tt=importlib.resources.files(basepath)
            for n in tt.iterdir():
                if os.path.isdir(n) and '__' not in str(n):
                    bn=str(n).split('/')[-1]
                    if bn in ResourceTypes:
                        self.libdir[bn]=n

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
            print(f'cwd: {self.cwd}')

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
        ''' default fetch behavior:  if filename exists in rootPath or projPath,
            take that one, otherwise look in the Library '''
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
                        print(f'Notice: copying {afname} over {destfname}')
                        os.system(f'cp {afname} {destpath}')
                    elif overwrite=='if_newer':
                        print(f'Notice: copying newer {afname} over {destfname}')
                        if os.path.getctime(destfname) < os.path.getctime(afname):
                            os.system(f'cp {afname} {destpath}')
                    else:
                        print(f'Notice: opting not to copy older {afname} over {destfname}')
                else:
                    print(f'Notice: copying {afname} into {destpath}')
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
                print(f'Notice: copying {name} over {destfullname}')
                os.system(f'cp {name} {destfullname}')
            elif overwrite=='if_newer':
                print(f'Notice: copying newer {name} over {destfullname}')
                if os.path.getctime(destfullname) < os.path.getctime(name):
                    os.system(f'cp {name} {destfullname}')
            else:
                print(f'Notice: opting not to copy older {name} over {destfullname}')
        else:
            print(f'Notice: copying {name} to {destfullname}')
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

