# -*- coding: utf-8 -*-
"""
Created on Thu Jun 21 13:40:55 2018

@author: HuangMing
"""
# =============================================================================
# Read origin top file and record top info into molecules class using list form
# =============================================================================


import sys
import top as Top

def rmLine(lst, key): #remove the line start with ';'
    new_lst = lst.copy()
    for item in lst:
        if item.startswith(key, 0, 5):
            new_lst.remove(item)
    return new_lst

def SeparateLst(inLst, key1, key2, end=False):
    if end:
        idxStart = inLst.index(key1)
        lst = inLst[idxStart:-1]
    else:
        idxStart = inLst.index(key1)
        idxEnd = inLst.index(key2)
        lst = inLst[idxStart:idxEnd]
    return lst

def InitTOP(fileName):
    f = open(fileName, 'r')
    line = f.read().split('\n')
    
#    key1 = '[ defaults ]'
#    key2 = '[ atomtypes ]'
#    key3 = '[ moleculetype ]'
    key4 = '[ atoms ]'
    key5 = '[ bonds ]'
    key6 = '[ pairs ]'
    key7 = '[ angles ]'
    key8 = '[ dihedrals ]'
#    key9 = '[ system ]'
#    key10 = '[ molecules ]'
#    default = SeparateLst(line, key1, key2)
#    atomTypes = SeparateLst(line, key2, key3)
#    moleculeTypes = SeparateLst(line, key3, key4)
    atoms = rmLine(SeparateLst(line, key4, key5), ';')
    bonds = rmLine(SeparateLst(line, key5, key6), ';')
    pairs = rmLine(SeparateLst(line, key6, key7), ';')
    angles = rmLine(SeparateLst(line, key7, key8), ';')
    dihedrals = rmLine(SeparateLst(line, key8, '', end=True), ';')
#    dihedrals = rmLine(SeparateLst(line, key8, key9), ';')
#    system = SeparateLst(line, key9, key10)
#    molecules = SeparateLst(line, key10, '', end=True)
    
#    return default, atomTypes, moleculeTypes, atoms, bonds, pairs, angles, dihedrals, system, molecules
    return atoms, bonds, pairs, angles, dihedrals

def FindAtoms(idx, atomsList): #just design for one molecule
    for atom in atomsList:
        idx1 = atom.getlocalIndex()
        if idx == idx1:
            return atom
    print('Didnt find desired atom, please check idx!')
    sys.exit()
    
def Atoms(atomsList, lst): #lst, stands for atoms content
    aList = []
    if len(atomsList) != len(lst):
        print('molecules not match, please check')
        sys.exit()
    for i in range(len(lst)):
        item = lst[i]
        atom = atomsList[i]
        atomType = item[1]
        charge = item[6]
        mass = item[7]
        atom.setAtomType(atomType)
        atom.setCharge(charge)
        atom.setMass(mass)
        aList.append(atom)
    return aList

def Bonds(atomsList, lst):
    bondsLst = []
    for item in lst:
        atom1Idx = item[0]
        atom2Idx = item[1]
        info = '{:>6}{:>10}{:>14}'.format(item[2], item[3], item[4])
        bondsLst.append([atom1Idx, atom2Idx, info])
    return bondsLst

def Pairs(atomsList, lst):
    pairsLst = []
    for item in lst:
        atom1Idx = item[0]
        atom2Idx = item[1]
        info = item[2]
        pairsLst.append([atom1Idx, atom2Idx, info])
    return pairsLst

def Angles(atomsList, lst):
    anglesLst = []
    for item in lst:
        atom1Idx = item[0]
        atom2Idx = item[1]
        atom3Idx = item[2]
        info = '{:>6}{:>12}{:>11}'.format(item[3], item[4], item[5])
        anglesLst.append([atom1Idx, atom2Idx, atom3Idx, info])
    return anglesLst

def Dihs(atomsList, lst):
    dihsLst = []
    for item in lst:
        atom1Idx = item[0]
        atom2Idx = item[1]
        atom3Idx = item[2]        
        atom4Idx = item[3]
        info = '{:>6}{:>11}{:>11}{:>3}'.format(item[4], item[5], item[6],item[7])
        dihsLst.append([atom1Idx, atom2Idx, atom3Idx, atom4Idx, info])
    return dihsLst

def splitLst(lst):
    new = []
    for line in lst: 
        if line != '':
            new.append(line.split())
    return new

def TopInfoInput(atomsList, cont, top):
    defaults = cont[0]
    atomTypes = cont[1]
    moleculeTypes = cont[2]
    atoms = splitLst(cont[3])[1:]
    bonds = splitLst(cont[4])[1:]
    pairs = splitLst(cont[5])[1:]
    angles = splitLst(cont[6])[1:]
    dihedrals = splitLst(cont[7])[1:]
    system = cont[8]
    molecules = cont[9]
    
    atomsList = Atoms(atomsList, atoms)
    bondsList = Bonds(atomsList, bonds)
    pairsList = Pairs(atomsList, pairs)
    anglesList = Angles(atomsList, angles)
    dihsList = Dihs(atomsList, dihedrals)
    
    top.setDefaults(defaults)
    top.setAtomTypes(atomTypes)
    top.setMoleculeTypes(moleculeTypes)
    top.setAtoms(atomsList)
    top.setBonds(bondsList)
    top.setPairs(pairsList)
    top.setAngles(anglesList)
    top.setDihedrals(dihsList)
    top.setSystem(system)
    top.setMolecules(molecules)
    
    return top

def main(atoms, mol, topName):
    cont = InitTOP(topName)
    top = Top.TopInfo()
    top_new = TopInfoInput(atoms, cont, top)
    mol.setTop(top_new)
    return mol
    
atoms, bonds, pairs, angles, dihedrals = InitTOP('init.itp')