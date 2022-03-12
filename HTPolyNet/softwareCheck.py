''' Check for presence of required software '''
import os
commandsRequired=['antechamber','tleap','parmchk2','gmx','obabel']

def CheckCommands():
    cnf=[]
    passes=True
    for c in commandsRequired:
        p=os.popen(f'which {c}').read()
        if p=='':
            passes=False
            cnf.append(c)
    if not passes:
        print('Error: HTPolyNet cannot find',', '.join(cnf))
        return False
    return True

def GetVersions():
    versions={}
    versions['gmx']=[s for s in os.popen(f'gmx -version').read().split('\n') if 'GROMACS version:' in s][0].split(':')[1].strip()
    versions['obabel']=os.popen('obabel -V').read().split()[2].strip()
    versions['ambertools']=[s for s in os.popen('antechamber -h').read().split('\n') if 'Welcome' in s][0].split()[3].strip().strip(':')
    return versions

