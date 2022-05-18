# _LIBRARY_DIRS_=['cfg','mdp','molecules/inputs','molecules/parameterized']
# _LIBRARY_EXT_DIR_={'cfg':['json','yaml'],'mdp':['mdp'],'molecules/inputs':['mol2'],'molecules/parameterized':['mol2','gro','top','itp','sea']}
# _LIBRARY_DIRS_WRITEABLE_={'cfg':True,'mdp':True,'molecules/input':False,'molecules/parameterized':True}
# def which_ldir(ext,mode=None):
#     for d,elist in _LIBRARY_EXT_DIR_.items():
#         if ext in elist:
#             if mode=='for_checkin':
#                 if _LIBRARY_DIRS_WRITEABLE_[d]:
#                     return d
#                 else:
#                     return None
#             elif mode=='for_checkout':
#                 return d
#             else:
#                 return d
#     return None