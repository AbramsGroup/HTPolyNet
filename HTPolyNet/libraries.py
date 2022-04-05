'''
Sets up access to resources in the Library package directory
'''
import importlib.resources
import os

from sympy import EX 

_DefaultResourceTypes_=['cfg','mdp','mol2']
_DefaultResourcePackageDir_='Library'

class Library:
    def __init__(self,ResourceTypes=_DefaultResourceTypes_,basepath=_DefaultResourcePackageDir_):
        self.types=ResourceTypes
        self.basepath=basepath
        self.libdir={}
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

    def fetch(self,filename,destpath='.'):
        pref,ext=filename.split('.')
        if not ext in self.libdir:
            raise Exception(f'no dir {ext} in Library')
        afname=os.path.join(self.dir(ext),filename)
        if os.path.exists(afname):
            os.system(f'cp {afname} {destpath}')
        else:
            raise FileNotFoundError(f'{afname} in {self.dir(ext)} not found.')