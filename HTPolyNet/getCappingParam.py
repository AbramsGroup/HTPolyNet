import os
import sys

#import HTPolyNet.configuration as configuration
import HTPolyNet.readTop2 as readTop2

def getCon(df_bonds, atIdx):
    df1 = df_bonds[(df_bonds.ai == atIdx) | (df_bonds.aj == atIdx)]
    atns = []
    for index, row in df1.iterrows():
        if row.ai == str(atIdx):
            atns.append(row.aj)
        elif row.aj == str(atIdx):
            atns.append(row.ai)
        else:
            raise TypeError('bond connection met error, please check!')
    return atns

def compareAType(atoms1, atoms2, atName, bonds):
    outLst = {}
    at1Type = atoms1.loc[(atoms1.atom == atName)].type.values[0]
    at1TypeUn = atoms2.loc[(atoms2.atom == atName)].type.values[0]
    if at1Type != at1TypeUn:
        atnRes = atoms1.loc[(atoms1.atom == atName)].residue.values[0]
        atnIdx = atoms2.loc[(atoms2.atom == atName)].nr.values[0]
        atnCharge = atoms2.loc[(atoms2.atom == atName)].charge.values[0]
        conAtns = getCon(bonds, atnIdx)
        outLst[atName] = {'charge': atnCharge, 'type': at1TypeUn}
        for atIdx in conAtns:
            row = atoms2.loc[(atoms2.nr == atIdx)]
            outLst[row.atom.values[0]] = {'charge': row.charge.values[0], 'type': row.type.values[0]}
            # outLst.append([atnRes, row.atom.values[0], row.type.values[0], row.charge.values[0]])
        return outLst
    else:
        return {}

def compareAtoms(topObj1, topObj2, cappingBonds):
    # topObj1: reactive part has been active
    # topObj2: reactive part doesn't active
    atoms1 = topObj1.atoms
    atoms2 = topObj2.atoms
    goalRes = str(atoms1['residue'].to_list()[0])
    outBonds = {}

    for bonds in cappingBonds:
        at1Name = bonds[1].strip()
        at2Name = bonds[2].strip()
        at1UnrctType = compareAType(atoms1, atoms2, at1Name, topObj2.bonds)
        at2UnrctType = compareAType(atoms1, atoms2, at2Name, topObj2.bonds)
        outBonds.update(at1UnrctType)
        outBonds.update(at2UnrctType)

    return goalRes, outBonds

def genUnrctMapping(inParam):
    unrctMap = {}
    for molPair in inParam.cappingMolPair:
        bonds = []
        mol1 = readTop2.initTop()
        mol2 = readTop2.initTop()
        mol1.setName('{}.top'.format(molPair[0].strip()), '{}.itp'.format(molPair[0].strip()))
        mol2.setName('{}.top'.format(molPair[1].strip()), '{}.itp'.format(molPair[1].strip()))
        mol1.genTopSession()
        mol2.genTopSession()

        for b in inParam.cappingBonds:
            if molPair[0].strip() in b[0].strip():
                bonds.append(b)
        resKey, resBonds = compareAtoms(mol1, mol2, bonds)
        unrctMap[resKey] = resBonds
    
    return unrctMap

if __name__ == '__main__':  # is this just for testing?
    oriPath = os.getcwd()
    sys.path.append(os.path.join(oriPath, 'mh-cl', 'src'))
    name = os.path.join(oriPath,'mh-cl', 'basic', 'options.txt-GMA-mST')
    a = configuration.configuration()
    a.setName(name)
    a.readParam()

    os.chdir(os.path.join(os.getcwd(), 'mh-cl', 'src', 'cappingFile', 'GMA-mST'))
    print(os.getcwd())
    
    '''
    finish getting the changing atom type
    TODO: 
    1. getting cappingBonds based on residue name
    2. put this in the work flow
    '''
    unrctMap = {}
    for molPair in a.cappingMolPair:
        bonds = []
        mol1 = readTop2.initTop()
        mol2 = readTop2.initTop()
        mol1.setName('{}.top'.format(molPair[0].strip()), '{}.itp'.format(molPair[0].strip()))
        mol2.setName('{}.top'.format(molPair[1].strip()), '{}.itp'.format(molPair[1].strip()))
        mol1.genTopSession()
        mol2.genTopSession()

        for b in a.cappingBonds:
            if molPair[0].strip() in b[0].strip():
                bonds.append(b)
        resKey, resBonds = compareAtoms(mol1, mol2, bonds)
        unrctMap[resKey] = resBonds
    print('unrctMap: ', unrctMap)
    
