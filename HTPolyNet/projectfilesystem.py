''' Classes for handling resource library and a runtime file system '''

import shutil
import logging
import os
import importlib.resources

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
        inst.package=libpackage
        try:
            tt=importlib.resources.files(inst.package)
        except:
            raise ImportError(f'Could not find package {libpackage}.  Your HTPolyNet installation is corrupt.')
        # make flat lists of absolute paths for all files
        for n in tt.iterdir():
            inst.allfiles.extend(inst._frec(n))
            inst.alldirs.extend(inst._drec(n))
        inst.root=os.path.commonpath(inst.alldirs)
        logging.info(inst.info(verbose=verbose))
        return inst

    @classmethod
    def user(cls,pathname='.'):
        assert os.path.exists(pathname),f'Cannot find {pathname} in {os.getcwd()}'
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

    def checkin(self,filename,overwrite=False):
        ''' filename must be a fully resolved pathname under the 
        system library.  We expect that the basename is in the current working directory'''
        basefilename=os.path.basename(filename)
        if not os.path.exists(basefilename):
            logging.info(f'{basefilename} not found in {os.getcwd()}. No check-in performed.')
            return False
        fullfilename=os.path.join(self.root,filename)
        if os.path.exists(fullfilename):
            if overwrite:
                shutil.copyfile(basefilename,fullfilename)
            else:
                logging.info(f'{filename} already exists in system library. No check-in performed.')
        else:
            shutil.copyfile(basefilename,fullfilename)
        return True

    def checkout(self,filename):
        basefilename=os.path.basename(filename)
        fullfilename=os.path.join(self.root,filename)
        if os.path.exists(fullfilename):
            shutil.copyfile(fullfilename,os.path.join(os.getcwd(),basefilename))
            return True
        else:
            logging.warning(f'Could not find {filename} in system library. No check-out performed.')
            return False
    
    def exists(self,filename):
        fullfilename=os.path.join(self.root,filename)
        return os.path.exists(fullfilename)

    def info(self,verbose=False):
        retstr=f'Library {self.designation} Directories:\n'
        for d in self.alldirs:
            retstr+=f'   {d}\n'
        if verbose:
            retstr+=f'Library {self.designation} Files:\n'
            for f in self.allfiles:
                retstr+=f'    {f}\n'
        return retstr

class ProjectFileSystem:
    def __init__(self,root='.',topdirs=['molecules','systems','plots'],verbose=False,reProject=False,userlibrary=None,mock=False):
        self.rootPath=os.path.abspath(root)
        os.chdir(self.rootPath)
        self.cwd=self.rootPath
        self.verbose=verbose
        if not mock:
            self._next_project_dir(reProject=reProject)
            self._setup_project_dir(topdirs=topdirs)
        self.library=RuntimeLibrary.system()
        self.userlibrary=None
        if userlibrary:
            self.userlibrary=RuntimeLibrary.user(userlibrary)

    def cdroot(self):
        os.chdir(self.rootPath)
        self.cwd=self.rootPath

    def cdproj(self):
        os.chdir(self.projPath)
        self.cwd=self.projPath

    def __str__(self):
        return f'root {self.rootPath}: cwd {self.cwd}'

    def _next_project_dir(self,reProject=False,prefix='proj-'):
        i=0
        lastprojdir=''
        currentprojdir=''
        while(os.path.isdir(os.path.join(self.rootPath,f'{prefix}{i}'))):
            lastprojdir=f'{prefix}{i}'
            i+=1
        if not reProject or lastprojdir=='': # this is a fresh project
            if lastprojdir=='':
                currentprojdir=f'{prefix}0'
            else:
                currentprojdir=f'{prefix}{i}'
            self.projPath=os.path.join(self.rootPath,currentprojdir)
            os.mkdir(currentprojdir)
        else:
            self.projPath=os.path.join(self.rootPath,lastprojdir)

    def _setup_project_dir(self,topdirs=['molecules','systems','plots']):
        os.chdir(self.projPath)
        self.projSubPaths={}
        for tops in topdirs:
            self.projSubPaths[tops]=os.path.join(self.projPath,tops)
            if not os.path.isdir(self.projSubPaths[tops]):
                os.mkdir(tops)

_PFS_=None

def pfs_setup(root='.',topdirs=['molecules','systems','plots'],
                verbose=False,reProject=False,userlibrary=None,mock=False):
    global _PFS_
    _PFS_=ProjectFileSystem(root=root,topdirs=topdirs,verbose=verbose,reProject=reProject,userlibrary=userlibrary,mock=mock)

def checkout(filename):
    if _PFS_.userlibrary and _PFS_.userlibrary.checkout(filename):
        return True
    return _PFS_.library.checkout(filename)

def exists(filename):
    if _PFS_.userlibrary and _PFS_.userlibrary.exists(filename):
        return True
    return _PFS_.library.exists(filename)

def checkin(filename,overwrite=False,priority='user'):
    if _PFS_.userlibrary and priority=='user':
        _PFS_.userlibrary.checkin(filename,overwrite=overwrite)
    else:
        _PFS_.library.checkin(filename,overwrite=overwrite)

def subpath(name):
    return _PFS_.projSubPaths[name]

def go_to(pathstr):
    """go_to Change the current working directory to "pathstr" which is assumed to be 
        under the project root.

    :param pathstr: pathname of directory relative to project root
    :type pathstr: _type_
    :return: absolute path of current working direcory
    :rtype: os.path
    """
    dirname=os.path.dirname(pathstr)
    if dirname=='':
        dirname=pathstr # assume this is a topdir
    assert dirname in _PFS_.projSubPaths,f'Error: cannot navigate using pathstring {pathstr}'
    _PFS_.cdproj()
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    os.chdir(dirname)
    basename=os.path.basename(pathstr)
    if basename!=pathstr:  # this is not a topdir
        if not os.path.exists(basename):
            logging.debug(f'PFS: making {basename}')
            os.mkdir(basename)
        os.chdir(basename)
    _PFS_.cwd=os.getcwd()
    return _PFS_.cwd

def root():
    return _PFS_.rootPath

def local_data_searchpath():
    return [_PFS_.rootPath,_PFS_.projPath]

def info():
    if _PFS_.userlibrary:
        print(f'User library is {str(_PFS_.userlibrary.path)}')
    print(_PFS_.library.info())

if __name__=='__main__':
    pfs=ProjectFileSystem(verbose=True)

