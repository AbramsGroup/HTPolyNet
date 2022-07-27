''' Classes for handling resource library and a runtime file system '''

import shutil
import logging
import os
import importlib.resources
logger=logging.getLogger(__name__)

class RuntimeLibrary:
    ''' a library object -- default creation references the Library resource package. '''
    def __init__(self):
        self.root=''
        self.subdirs=[]
        self.designation=''

    @classmethod
    def system(cls,libpackage='Library'):
        inst=cls()
        inst.designation='System'
        # inst.package=libpackage
        try:
            with importlib.resources.path(libpackage,'__init__.py') as f:
                inst.root=os.path.split(os.path.abspath(f))[0]
            subdirs=os.listdir(inst.root)
            for xxx in ['__init__.py','__pycache__','README.md']:
                subdirs.remove(xxx)
            inst.subdirs=[os.path.join(inst.root,xxx) for xxx in subdirs]
        except:
            raise ImportError(f'Could not find package {libpackage}.  Your HTPolyNet installation is corrupt.')
        logger.info(inst.info())
        return inst

    @classmethod
    def user(cls,pathname='.'):
        if not pathname:
            return None
        assert os.path.exists(pathname),f'Cannot find {pathname} in {os.getcwd()}'
        tt=os.path.abspath(pathname)
        assert os.path.isdir(tt),f'Please ensure that {str(tt)} is a directory'
        inst=cls()
        inst.designation='User'
        inst.root=tt
        for x in ['__init__.py', 'README.md', '__pycache__']:
            if x in tt:
                tt.remove(x)
        inst.subdirs=[x for x in tt if os.path.isdir(x)]
        logger.info(inst.info())
        return inst

    def checkin(self,filename,overwrite=False):
        ''' filename must be a fully resolved pathname under the 
        system library.  We expect that the basename is in the current working directory'''
        basefilename=os.path.basename(filename)
        if not os.path.exists(basefilename):
            logger.info(f'{basefilename} not found in {os.getcwd()}. No check-in performed.')
            return False
        fullfilename=os.path.join(self.root,filename)
        if os.path.exists(fullfilename):
            if overwrite:
                shutil.copyfile(basefilename,fullfilename)
            else:
                logger.info(f'{filename} already exists in system library. No check-in performed.')
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
            # logger.warning(f'Could not find {filename} in system library. No check-out performed.')
            return False
    
    def exists(self,filename):
        fullfilename=os.path.join(self.root,filename)
        return os.path.exists(fullfilename)

    def info(self):
        return f'{self.designation} library is {self.root}'

_SYSTEM_LIBRARY_=None
def lib_setup():
    global _SYSTEM_LIBRARY_
    if _SYSTEM_LIBRARY_==None:
        _SYSTEM_LIBRARY_=RuntimeLibrary.system()
    return _SYSTEM_LIBRARY_

class ProjectFileSystem:
    def __init__(self,root='.',topdirs=['molecules','systems','plots'],verbose=False,reProject=False,userlibrary=None,mock=False):
        self.library=lib_setup()
        self.userlibrary=None
        if userlibrary:
            self.userlibrary=RuntimeLibrary.user(userlibrary)
        self.rootPath=os.path.abspath(root)
        os.chdir(self.rootPath)
        self.cwd=self.rootPath
        self.verbose=verbose
        if not mock:
            self._next_project_dir(reProject=reProject)
            self._setup_project_dir(topdirs=topdirs)
        

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

_PFS_:ProjectFileSystem=None

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
            logger.debug(f'PFS: making {basename}')
            os.mkdir(basename)
        os.chdir(basename)
    _PFS_.cwd=os.getcwd()
    return _PFS_.cwd

def root():
    return _PFS_.rootPath

def proj():
    return _PFS_.projPath

def local_data_searchpath():
    return [_PFS_.rootPath,_PFS_.projPath]

def info():
    if _PFS_.userlibrary:
        print(_PFS_.userlibrary.info())
    print(_PFS_.library.info())

if __name__=='__main__':
    pfs=ProjectFileSystem(verbose=True)

