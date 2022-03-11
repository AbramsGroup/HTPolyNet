'''
Sets up access to resources in the Templates directory
'''
import importlib.resources
import os 

_DefaultResourceTypes_=['cfg','Gromacs_mdp','mol2']
_DefaultResourcePackageDir_='Templates'

def IdentifyTemplateResourcePaths(ResourceTypes=_DefaultResourceTypes_,basedir=_DefaultResourcePackageDir_):
    td={}
    tt=importlib.resources.files(basedir)
    for n in tt.iterdir():
        if os.path.isdir(n) and '__' not in str(n):
            bn=str(n).split('/')[-1]
            if bn in ResourceTypes:
                td[bn]=n
    return td
