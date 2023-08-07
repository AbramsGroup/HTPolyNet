"""

.. module:: projectfilesystem
   :synopsis: handles project filesystems
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import shutil
import logging
import os
import glob
from . import resources
logger=logging.getLogger(__name__)
excludes=['__pycache__','__init__.py']
class RuntimeLibrary:
    ''' a library object -- default creation references the Library resource package. '''
    def __init__(self):
        self.root=''
        # self.subdirs=[]
        self.designation=''
        self.ResourcePaths={}
    @classmethod
    def system(cls):
        """system generates a RuntimeLibrary object corresponding to the installed Library subpackage

        :param libpackage: name of the HTPolyNet Library subpackage, defaults to 'Library'
        :type libpackage: str, optional
        :raises ImportError: if the system library subpackage is not found
        :return: a RuntimeLibrary object holding the system library subpackage
        :rtype: RuntimeLibrary
        """
        inst=cls()
        inst.designation='System'
        inst.root=os.path.dirname(resources.__file__)
        ResourceFullPaths=glob.glob(inst.root+'/*')
        for l in ResourceFullPaths:
            bn=os.path.basename(l)
            if not bn in excludes:
                inst.ResourcePaths[bn]=l
        logger.info(inst.info())
        return inst

    def get_example_depot_location(self):
        """get_example_depot_location reports the location of the examples in this RuntimeLibrary, if they exist

        :return: name of depot directory
        :rtype: str
        """
        assert 'example_depot' in self.ResourcePaths,f'No example depot found in {self.root} -- your installation is corrupt.'
        return self.ResourcePaths['example_depot']

    def get_example_names(self):
        """get_example_names returns the names of the example tarballs

        :return: list of names of example tarballs
        :rtype: list
        """
        depot=self.get_example_depot_location()
        assert os.path.exists(depot) and os.path.isdir(depot),f'Depot not found or not a directory -- your installation is corrupt.'
        owd=os.getcwd()
        os.chdir(depot)
        example_tarballs=glob.glob('*.tgz')
        basenames=[x.replace('.tgz','') for x in example_tarballs]
        os.chdir(owd)
        return basenames

    @classmethod
    def user(cls,pathname='.'):
        """user generates a new user-level RuntimeLibrary object

        :param pathname: where to find the library, defaults to '.'
        :type pathname: str, optional
        :return: a new RuntimeLibrary object
        :rtype: RuntimeLibrary
        """
        if not pathname:
            return None
        assert os.path.exists(pathname),f'Cannot find {pathname} in {os.getcwd()}'
        tt=os.path.abspath(pathname)
        assert os.path.isdir(tt),f'Please ensure that {str(tt)} is a directory'
        inst=cls()
        inst.designation='User'
        ResourceFullPaths=glob.glob(inst.root+'/*')
        for l in ResourceFullPaths:
            bn=os.path.basename(l)
            if not bn in excludes:
                inst.ResourcePaths[bn]=l
        return inst
    
    def checkin(self,filename,overwrite=False):
        """checkin check filename into this RuntimeLibrary

        :param filename: name of file to check in, ***relative to the top level directory in which the library is housed***
        :type filename: str
        :param overwrite: if True, overwrite the file in the RuntimeLibrary if it exists, defaults to False
        :type overwrite: bool, optional
        :return: False if filename not found in cwd; True if check-in was successful
        :rtype: bool
        """
        basefilename=os.path.basename(filename)
        if not os.path.exists(basefilename):
            logger.debug(f'{basefilename} not found in {os.getcwd()}. No check-in performed.')
            return False
        fullfilename=os.path.join(self.root,filename)
        if os.path.exists(fullfilename):
            if overwrite:
                shutil.copyfile(basefilename,fullfilename)
            else:
                logger.debug(f'{filename} already exists in system library. No check-in performed.')
        else:
            shutil.copyfile(basefilename,fullfilename)
        return True

    def checkout(self,filename,searchpath=[],altpath=[]):
        """checkout check the file named by filename out of this RuntimeLibrary and copy it to the cwd

        :param filename: name of file to checkout, relative to the top level directory in which the library is found
        :type filename: str
        :param searchpath: path of directories to search if the filename is not resolved, defaults to []
        :type searchpath: list, optional
        :param altpath: additional directories to add to the search path, defaults to []
        :type altpath: list, optional
        :return: True if file checkout was successful, False if not
        :rtype: bool
        """
        basefilename=os.path.basename(filename)
        fullfilename=os.path.join(self.root,filename)
        if os.path.exists(fullfilename):
            shutil.copyfile(fullfilename,os.path.join(os.getcwd(),basefilename))
            return True
        else:
            # logger.info(f'No {filename} found in libraries; checking local data searchpath {searchpath}')
            if altpath:
                searchpath.append(altpath)
                # logger.info(f'and alternative path {altpath}')
            for p in searchpath:
                # logger.debug(f'Searching {p}...')
                fullfilename=os.path.join(p,filename)
                if os.path.exists(fullfilename):
                    basefilename=os.path.basename(filename)
                    shutil.copyfile(fullfilename,basefilename)
                    return True
            # logger.debug(f'Could not find {filename} anywhere in user library!')
            return False
    
    def exists(self,filename):
        """exists checks to see if filename exists in this RuntimeLibrary

        :param filename: name of file to check for, relative to toplevel directory of RuntimeLibrary
        :type filename: str
        :return: True if filename exists in RuntimeLibrary, False otherwise
        :rtype: bool
        """
        fullfilename=os.path.join(self.root,filename)
        return os.path.exists(fullfilename)

    def info(self):
        """info returns a simple string describing this RuntimeLibrary

        :return: the string
        :rtype: str
        """
        return f'{self.designation} library is {self.root}'

_SYSTEM_LIBRARY_=None
def lib_setup():
    """lib_setup sets up the system RuntimeLibrary

    :return: system RuntimeLibrary object
    :rtype: RuntimeLibrary
    """     
    global _SYSTEM_LIBRARY_
    if _SYSTEM_LIBRARY_==None:
        _SYSTEM_LIBRARY_=RuntimeLibrary.system()
    return _SYSTEM_LIBRARY_
def system():
    """system returns the system RuntimeLibrary object

    :return: system RuntimeLibrary object
    :rtype: RuntimeLibrary
    """
    return _SYSTEM_LIBRARY_

class ProjectFileSystem:
    """ Handles all aspects of the creation and organization of a project filesystem
    """
    def __init__(self,root='.',topdirs=['molecules','systems','plots'],projdir='next',verbose=False,reProject=False,userlibrary=None,mock=False):
        """__init__ Generates a new ProjectFilesystem object

        :param root: path of root directory in which the project directory is to be housed, defaults to '.'
        :type root: str, optional
        :param topdirs: names of toplevel directories in the project directory, defaults to ['molecules','systems','plots']
        :type topdirs: list, optional
        :param projdir: name of the project directory, or 'next' if next available automatically generated name is to be used, defaults to 'next'
        :type projdir: str, optional
        :param verbose: flag to turn on verbose output, defaults to False
        :type verbose: bool, optional
        :param reProject: flag indicating if this is a restarted project, defaults to False
        :type reProject: bool, optional
        :param userlibrary: name of user library toplevel directory, relative to root, defaults to None
        :type userlibrary: str, optional
        :param mock: flag indicating if this is a mock call, defaults to False
        :type mock: bool, optional
        """
        self.library=lib_setup()
        self.userlibrary=None
        if userlibrary:
            self.userlibrary=RuntimeLibrary.user(userlibrary)
        self.rootPath=os.path.abspath(root)
        os.chdir(self.rootPath)
        self.cwd=self.rootPath
        self.verbose=verbose
        if not mock:
            self._next_project_dir(projdir=projdir,reProject=reProject)
            self._setup_project_dir(topdirs=topdirs)
        
    def cdroot(self):
        """cdroot change the cwd to the root (parent directory of project directory)
        """
        os.chdir(self.rootPath)
        self.cwd=self.rootPath

    def cdproj(self):
        """cdproj change the cwd to the toplevel project directory
        """
        os.chdir(self.projPath)
        self.cwd=self.projPath

    def go_to(self,subPath,make=False):
        """go_to change the cwd to the directory named by 'subPath'

        :param subPath: directory to change to, relative to project directory
        :type subPath: str
        """
        self.cdproj()
        if os.path.exists(subPath):
            os.chdir(subPath)
            self.cwd=os.getcwd()
        elif make:
            os.mkdir(subPath)
            os.chdir(subPath)
            self.cwd=os.getcwd()

    def __str__(self):
        return f'root {self.rootPath}: cwd {self.cwd}'

    def _next_project_dir(self,projdir='next',reProject=False,prefix='proj-'):
        """_next_project_dir sets the name of the current project directory and creates it if necessary

        :param projdir: name of next project directory, defaults to 'next'
        :type projdir: str, optional
        :param reProject: flag indicating if this is a restart, defaults to False
        :type reProject: bool, optional
        :param prefix: prefix for automatically generated project directories, defaults to 'proj-'
        :type prefix: str, optional
        """
        if not projdir=='next':  # explicit project directory is named
            # if os.path.exists(projdir) and not reProject:
            #     raise Exception(f'Project directory {projdir} exists but you did not indicate you wanted to restart with the "-restart" flag')
            self.projPath=os.path.join(self.rootPath,projdir)
            if os.path.exists(projdir):
                logger.info(f'Working in existing project {self.projPath}')
            else:
                os.mkdir(projdir)
                logger.info(f'Working in new project {self.projPath}')
        else:
            i=0
            lastprojdir=''
            currentprojdir=''
            while(os.path.isdir(os.path.join(self.rootPath,f'{prefix}{i}'))):
                lastprojdir=f'{prefix}{i}'
                logger.debug(f'{lastprojdir} exists')
                i+=1
            assert not os.path.exists(f'{prefix}{i}')
            if not reProject or lastprojdir=='': # this is a fresh project
                if lastprojdir=='':
                    currentprojdir=f'{prefix}0'
                else:
                    currentprojdir=f'{prefix}{i}'
                self.projPath=os.path.join(self.rootPath,currentprojdir)
                logger.info(f'New project in {self.projPath}')
                os.mkdir(currentprojdir)
            else:
                self.projPath=os.path.join(self.rootPath,lastprojdir)
                logger.info(f'Restarting project in {self.projPath} (latest project)')

    def _setup_project_dir(self,topdirs=['molecules','systems','plots']):
        """_setup_project_dir sets up the project directory after it is created by creating the requested top-level subdirectories

        :param topdirs: top-level subdirectories, defaults to ['molecules','systems','plots']
        :type topdirs: list, optional
        """
        os.chdir(self.projPath)
        self.projSubPaths={}
        for tops in topdirs:
            self.projSubPaths[tops]=os.path.join(self.projPath,tops)
            if not os.path.isdir(self.projSubPaths[tops]):
                os.mkdir(tops)

_PFS_:ProjectFileSystem=None

def pfs_setup(root='.',topdirs=['molecules','systems','plots'],projdir='next',
                verbose=False,reProject=False,userlibrary=None,mock=False):
    """pfs_setup sets up the global ProjectFileSystem object

    :param root: parent directory of this project file system, defaults to '.'
    :type root: str, optional
    :param topdirs: top-level subdirectories, defaults to ['molecules','systems','plots']
    :type topdirs: list, optional
    :param projdir: name of the project directory itself, defaults to 'next'
    :type projdir: str, optional
    :param verbose: flag indicating verbose output, defaults to False
    :type verbose: bool, optional
    :param reProject: flag indicating restart, defaults to False
    :type reProject: bool, optional
    :param userlibrary: name of user library relative to root, defaults to None
    :type userlibrary: str, optional
    :param mock: flag indicating this is a mock call, defaults to False
    :type mock: bool, optional
    """
    global _PFS_
    _PFS_=ProjectFileSystem(root=root,topdirs=topdirs,projdir=projdir,verbose=verbose,reProject=reProject,userlibrary=userlibrary,mock=mock)

def checkout(filename,altpath=[]):
    """checkout handles checking files out of the library owned by the global ProjectFileSystem object; tries to checkout from user-level local library first before trying to check out from the global system Library subpackage

    :param filename: name of file to check out relative to top-level of library
    :type filename: str
    :param altpath: list of alternate directories to search, defaults to []
    :type altpath: list, optional
    :return: True if checkout from user library successful, otherwise returns result of checkout from global system Library (also a bool)
    :rtype: bool
    """
    if _PFS_.userlibrary and _PFS_.userlibrary.checkout(filename,searchpath=[_PFS_.rootPath,_PFS_.projPath],altpath=altpath):
        return True
    return _PFS_.library.checkout(filename)

def fetch_molecule_files(mname):
    """fetch_molecule_files fetches all relevant molecule data files for molecule named 'mname'

    :param mname: name of molecule
    :type mname: str
    :return: list of filetypes found and fetched
    :rtype: list
    """
    ret_exts=[]
    dirname='molecules/parameterized'
    for e in ['mol2','pdb','gro','top','itp','grx']:
        prob_filename=os.path.join(dirname,f'{mname}.{e}')
        if exists(prob_filename):
            ret_exts.append(e)
            checkout(prob_filename)
    return ret_exts

def exists(filename):
    """exists checks for existence of filename in user-level local library, then in the system Library

    :param filename: name of file
    :type filename: str
    :return: True if found in either library
    :rtype: bool
    """
    # check user library first
    if _PFS_.userlibrary and _PFS_.userlibrary.exists(filename):
        return True
    # check system library
    return _PFS_.library.exists(filename)

def checkin(filename,overwrite=False,priority='user'):
    """checkin check the file named 'filename' into either the user-level local library or the global system Library package (the latter requires write priveleges wherever Python packages are stored if HTPolyNet is installed without -e or from PyPI)

    :param filename: name of file to check in 
    :type filename: str
    :param overwrite: flag indicates whether file should overwrite file of same name in library, defaults to False
    :type overwrite: bool, optional
    :param priority: string indicating which library to check into, defaults to 'user'
    :type priority: str, optional
    """
    if _PFS_.userlibrary and priority=='user':
        _PFS_.userlibrary.checkin(filename,overwrite=overwrite)
    else:
        _PFS_.library.checkin(filename,overwrite=overwrite)

def subpath(name):
    """subpath returns the path of the project subdirectory with name 'name'

    :param name: name of subdirectory
    :type name: str
    :return: path of subdirectory
    :rtype: os.path
    """
    return _PFS_.projSubPaths[name]

def go_proj():
    """go_proj change the current working directory to the project directory
    """
    _PFS_.cdproj()

def go_root():
    """go_root change the current working directory to the root directory
    """
    _PFS_.cdroot()

def go_to(pathstr):
    """go_to Change the current working directory to "pathstr" which is relative to the project root.

    :param pathstr: pathname of directory relative to project root
    :type pathstr: str
    :return: absolute path of current working direcory
    :rtype: os.path
    """
    _PFS_.cdproj()
    dirname=os.path.dirname(pathstr)
    if dirname=='':
        dirname=pathstr # assume this is a topdir
    assert dirname in _PFS_.projSubPaths,f'Error: cannot navigate using pathstring {pathstr}'
    reentry=True
    if not os.path.exists(dirname):
        os.mkdir(dirname)
        reentry=False
    os.chdir(dirname)
    basename=os.path.basename(pathstr)
    if basename!=pathstr:  # this is not a topdir
        if not os.path.exists(basename):
            logger.debug(f'PFS: making {basename}')
            os.mkdir(basename)
        os.chdir(basename)
    _PFS_.cwd=os.getcwd()
    return reentry

def root():
    """root returns name of root directory (parent of project directory)

    :return: name of root directory
    :rtype: os.path
    """
    return _PFS_.rootPath

def cwd():
    """cwd returns name of current working directory relative to the root path

    :return: name of current working directory relative to the root path
    :rtype: os.path
    """
    return os.path.relpath(os.getcwd(),start=_PFS_.rootPath)

def proj():
    """proj returns name of project directory

    :return: name of project directory
    :rtype: os.path
    """
    return _PFS_.projPath

def local_data_searchpath():
    """local_data_searchpath returns container containing root path and project path

    :return: list containing root path and project path
    :rtype: list
    """
    return [_PFS_.rootPath,_PFS_.projPath]

def info():
    """info prints some summary information about Libraries to the console
    """
    if _PFS_.userlibrary:
        print(_PFS_.userlibrary.info())
    print(_PFS_.library.info())

def proj_abspath(filename):
    """proj_abspath returns the path of the file named filename relative to the project directory

    :param filename: name of the probe file
    :type filename: str
    :return: name of probe file relative to the project directory
    :rtype: os.path
    """
    abf=os.path.abspath(filename)
    return os.path.relpath(abf,_PFS_.projPath)

