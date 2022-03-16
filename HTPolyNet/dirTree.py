''' subdirectory tree for a run directory '''
import os
import sys

'''
root
   |----> proj0
   |----> proj1
   | ...

projn
    |----> basic
    |----> mdp
    |----> systems
    |            |----> unrctSystems
    |            |----> rctSystems
    |            |----> typeSystems
    +----> results
                 |----> init
                 |----> sim0
                 |----> sim1
                 ...

simn
   |----> init


'''

def dirTree(root='.',reProject='',nReplicas=1):
    R={}
    R['root']=root
    i=0
    if reProject=='': # this is a fresh project
        while(os.path.isdir(os.path.join(root,f'proj{i}'))):
            i+=1
        projPath=os.path.join(root,f'proj{i}')
    else:
        if not os.path.isdir(os.path.join(root,reProject)):
            raise FileNotFoundError(f'{reProject} not found')
        projPath=os.path.join(root,reProject)    
    R['project']=projPath
    toplevels=['basic','mdp','results','systems']
    for tl in toplevels:
        R[tl]=os.path.join(projPath,tl)
    syssubs=['unrctSystem','rctSystem','typeSystem']
    for ss in syssubs:
        R[ss]=os.path.join(R['systems'],ss)
    ressubs=['init']
    for rs in ressubs:
        R[rs]=os.path.join(R['results'],rs)
    return R

