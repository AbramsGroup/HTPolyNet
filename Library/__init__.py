_LIBRARY_DIRS_=['cfg','mdp','molecules']
_LIBRARY_EXT_DIR_={'cfg':['json','yaml'],'mdp':['mdp'],'molecules':['mol2','gro','top','itp','sea']}
def which_ldir(ext):
    for d,elist in _LIBRARY_EXT_DIR_.items():
        if ext in elist:
            return d
    return None