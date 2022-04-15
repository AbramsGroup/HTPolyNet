''' Classes for handling resource library and a runtime file system '''

import logging
import os
import pathlib
import importlib.resources

from HTPolyNet.software import Command

_MolecularDataFileTypes_=['mol2','gro','top','itp','sea']
_MoleculeClasses_=['monomers','oligomers']

class Library:
    ''' a library object -- default creation references the Library resource package. '''
    def __init__(self):
        self.alldirs=[]
        self.allfiles=[]
        self.molecule_class_subdirs={}
        self.designation='Empty'
    
    @classmethod
    def system(cls,libpackage='Library',verbose=False):
        inst=cls()
        inst.designation='System'
        try:
            tt=importlib.resources.files(libpackage)
        except:
            raise ImportError(f'Could not find package {libpackage}.  Your HTPolyNet installation is corrupt.')
        for n in tt.iterdir():
            inst.allfiles.extend(inst._frec(n))
            inst.alldirs.extend(inst._drec(n))
        for c in _MoleculeClasses_:
            dir=[x for x in inst.alldirs if c in str(x)]
            if len(dir)==0:
                raise Exception(f'No "molecules/{dir}" directory in {libpackage}')
            inst.molecule_class_subdirs[c]=dir[0]
        logging.info(inst.info(verbose=verbose))
        return inst

    @classmethod
    def user(cls,pathname='.'):
        inst=cls()
        inst.designation='User'
        tt=os.path.abspath(pathname)
        logging.warning('User-level data libraries are not yet implemented. Sorry!')
        if not tt.is_dir():
            logging.error('error')
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
            if n!='__pycache__':
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
        if os.path.sep in basefilename:
            logging.warning(f'basefilename {basefilename} seems to have path separators.  It should be a base name only.')
            basefilename=os.path.basename(basefilename)
            logging.warning(f'...interpreting as {basefilename}')
        matches = [x for x in self.allfiles if basefilename in str(x)]
        if len(matches)==0:
            ''' this file is not found, but we will return the path where it should go '''
            root,ext=os.path.splitext(basefilename)
            ext=ext.replace('.','')
            if not ext in _MolecularDataFileTypes_:
                ''' this is not a molecular data file '''
                matches = [x for x in self.alldirs if ext in str(x)]
                if len(matches)==0:
                    return False
                return matches[0]
            else:
                ''' this is a molecular data file '''
                ''' an oligomer will have @ in its file name '''
                if '@' in basefilename:
                    return self.molecule_class_subdirs['oligomers']
                else:
                    return self.molecule_class_subdirs['monomers']
        return matches[0]

    def checkin(self,basefilename,overwrite=False):
        if not os.path.exists(basefilename):
            logging.info(f'{basefilename} not found. No check-in performed.')
            return False
        ''' add the local file basefilename to the Library if it is not already there '''
        p=self.path(basefilename)
        if p.is_file():
            if overwrite:
                out,err=Command(f'cp -f {basefilename} {p}').run()
                return True
            else:
                logging.info(f'Check-in of {basefilename} to {self.designation} library is not necessary.')
                logging.info(f'    {p} is already there.')
        elif p.is_dir():
            out,err=Command(f'cp {basefilename} {p}').run()
            self.allfiles.append(os.path.join(p,basefilename))
            return True
        return False  # this would be the result of an unspecified error

    def checkout(self,basefilename,nowarn=False):
        p=self.path(basefilename)
        if p and p.is_file():
            logging.info(f'Checking {basefilename} out of {self.designation} library into {os.getcwd()}')
            Command(f'cp {p} .').run()
            return True
        else:
            if not nowarn:
                logging.info(f'{basefilename} is not found in the {self.designation} library')
                if p and p.is_dir():
                    logging.info(f'    but I expected to find it in {p}.')
                else:
                    root,ext=os.path.splitext(basefilename)
                    logging.info(f'    because I have no information on where {ext} files should go.')
            return False
    
    def exists(self,basefilename):
        p=self.path(basefilename)
        if p and p.is_file():
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
        self.library=Library.system()
        self.userlibrary=userlibrary
        if userlibrary:
            self.userlibrary=Library.user(userlibrary)

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

