''' Classes for handling resource library and a runtime file system '''

from email.mime import base
import logging
import os
import pathlib
import importlib.resources

from HTPolyNet.command import Command
from Library import _LIBRARY_EXT_DIR_, which_ldir

_MolecularDataFileTypes_=['mol2','gro','top','itp','sea']
_mdf_last_required_=_MolecularDataFileTypes_.index('sea')

class RuntimeLibrary:
    ''' a library object -- default creation references the Library resource package. '''
    def __init__(self):
        self.alldirs=[]
        self.allfiles=[]
        self.designation='Empty'
    
    @classmethod
    def system(cls,libpackage='Library',verbose=False):
        inst=cls()
        inst.designation='System'
        try:
            tt=importlib.resources.files(libpackage)
        except:
            raise ImportError(f'Could not find package {libpackage}.  Your HTPolyNet installation is corrupt.')
        # make flat lists of absolute paths for all files
        for n in tt.iterdir():
            inst.allfiles.extend(inst._frec(n))
            inst.alldirs.extend(inst._drec(n))
        inst.ext_dirs=which_ldir
        basenames=[os.path.basename(x) for x in inst.allfiles]
        basenameset=set(basenames)
        for b in basenameset:
            if basenames.count(b)>1:
                logging.warning(f'Two or more files with basename {b} are detected in the system library.')
                logging.warning(f'checkout/checkin will resolve to {inst.allfiles[basenames.index(b)]}')
        logging.info(inst.info(verbose=verbose))
        return inst

    @classmethod
    def user(cls,pathname='.'):
        assert os.path.exists(pathname),f'Please ensure {pathname} exists.'
        tt=os.path.abspath(pathname)
        assert tt.is_dir(),f'Please ensure that {str(tt)} is a directory'
        inst=cls()
        inst.designation='User'
        inst.root=tt
        return inst

    def _frec(self,n):
        if n.is_file():
            return [n]
        elif n.is_dir():
            ret=[]
            for m in n.iterdir():
                ret.extend(self._frec(m))
            return ret
        else:
            return []

    def _drec(self,n):
        if n.is_dir():
            if '__pycache__' not in str(n):
                ret=[n]
                for m in n.iterdir():
                    ret.extend(self._drec(m))
                return ret
            else:
                return []
        else:
            return []

    def path(self,basefilename):
        ''' return the absolute path of the provided basename in the Library, or, if the basename is not in the library, return the path to the subdirectory it *should* go in.  Return False if no such path exists in the Library. '''
        assert not os.path.sep in basefilename
        matches = [x for x in self.allfiles if basefilename in str(x)]
        if len(matches)==1:
            return matches[0]
        elif len(matches)==0:
            prefix,ext=os.path.splitext(basefilename)
            ext=ext[1:]
            return self.ext_dirs(ext)
        else:
            logging.warning(f'More than one file {basefilename} found in library.  Only returning {matches[0]}')
            return matches[0]

    def checkin(self,basefilename,overwrite=False):
        if not os.path.exists(basefilename):
            logging.info(f'{basefilename} not found. No check-in performed.')
            return False
        ''' add the local file basefilename to the Library if it is not already there '''
        p=self.path(basefilename)
        if p and os.path.isfile(p):
            if overwrite:
                out,err=Command(f'cp -f {basefilename} {p}').run()
                return True
            else:
                logging.info(f'Check-in of {basefilename} to {self.designation} library is not necessary.')
                logging.info(f'    {p} is already there.')
        elif p and os.path.isdir(p):
            out,err=Command(f'cp {basefilename} {p}').run()
            self.allfiles.append(os.path.join(p,basefilename))
            return True
        logging.warning(f'Failed to check in {basefilename} -- I do not know where to put it.')
        return False

    def checkout(self,basefilename,nowarn=False):
        p=self.path(basefilename)
        if p and os.path.isfile(p):
            logging.info(f'Checking {basefilename} out of {self.designation} library into {os.getcwd()}')
            Command(f'cp {p} .').run()
            return True
        else:
            if not nowarn:
                logging.info(f'{basefilename} is not found in the {self.designation} library')
                if p and os.path.isdir(p):
                    logging.info(f'    but I expected to find it in {p}.')
                else:
                    root,ext=os.path.splitext(basefilename)
                    logging.info(f'    because I have no information on where {ext} files should go.')
            return False
    
    def exists(self,basefilename):
        p=self.path(basefilename)
        if p and os.path.isfile(p):
            return True
        return False

    def info(self,verbose=False):
        retstr=f'HTPolyNet Library {self.designation} Directories:\n'
        for d in self.alldirs:
            retstr+=f'   {d}\n'
        if verbose:
            retstr+=f'HTPolyNet Library {self.designation} Files:\n'
            for f in self.allfiles:
                retstr+=f'    {f}\n'
        return retstr

class ProjectFileSystem:
    def __init__(self,root='.',verbose=False,reProject=False,userlibrary=None):
        self.rootPath=pathlib.Path(root).resolve()
        self.cwd=self.rootPath
        self.verbose=verbose
        self.cd(self.rootPath)
        self._next_project_dir(reProject=reProject)
        self._setup_project_root()
        self.cd(self.rootPath)
        self.library=RuntimeLibrary.system()
        self.userlibrary=None
        if userlibrary:
            self.userlibrary=RuntimeLibrary.user(userlibrary)

    def cd(self,dest=''):
        os.chdir(dest)
        self.cwd=os.getcwd()
        if (self.verbose):
            logging.info(f'cwd: {self.cwd}')
        return self.cwd

    def cdroot(self):
        return self.cd(self.projPath)

    def cdrootsub(self,toplevel):
        if toplevel in self.projSubPaths:
            return self.cd(self.projSubPaths[toplevel])

    def exists(self,name):
        ''' look in user library first, if it exists '''
        if self.userlibrary:
            fullname=os.path.join(self.userlibrary.tt,name)
            if os.path.exists(fullname):
                return True
        return self.library.exists(name)

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
        for tops in ['rctSystems','unrctSystems']:
            self.systemsSubPaths[tops]=os.path.join(self.systemsPath,tops)
            if not os.path.isdir(self.systemsSubPaths[tops]):
                os.mkdir(tops)
        self.unrctPath=self.systemsSubPaths['unrctSystems']
        self.rctPath=self.systemsSubPaths['rctSystems']

_PFS_=None

def pfs_setup(root='.',verbose=False,reProject=False,userlibrary=None):
    global _PFS_
    _PFS_=ProjectFileSystem(root=root,verbose=verbose,reProject=reProject,userlibrary=userlibrary)

def checkout(filename):
    _PFS_.library.checkout(filename)

def exists(filename):
    return _PFS_.exists(filename)

def in_library(filename):
    in_system=all([exists(f'{filename}.{ex}') for ex in _MolecularDataFileTypes_[:_mdf_last_required_]])
    logging.info(f'in_library for {filename}: {in_system}')
    return in_system

def checkin(filename):
    _PFS_.library.checkin(filename)

def cd(pathstring):
    if pathstring=='root':
        return _PFS_.cdroot()
    cats=[_PFS_.projSubPaths,_PFS_.systemsSubPaths,_PFS_.resultsSubPaths]
    for cat in cats:
        if pathstring in cat:
            _PFS_.cd(cat[pathstring])
            return _PFS_.cwd
    raise Exception(f'Path {pathstring} cannot be located in the project file system')
    
def next_results_dir(maxstages=100):
    possibles=['init',*[f'step{i}' for i in range(maxstages)]]
    for p in possibles:
        if not p in _PFS_.resultsSubPaths:
            break
    newpath=os.path.join(_PFS_.resPath,p)
    os.mkdir(newpath)
    _PFS_.resultsSubPaths[p]=newpath
    return cd(newpath)

def local_data_searchpath():
    return [_PFS_.rootPath,_PFS_.projPath]

def info():
    if _PFS_.userlibrary:
        print(f'User library is {str(_PFS_.userlibrary.path)}')
    print(_PFS_.library.info())

if __name__=='__main__':
    pfs=ProjectFileSystem(verbose=True)

