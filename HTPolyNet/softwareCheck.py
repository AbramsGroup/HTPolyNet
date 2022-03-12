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
     

