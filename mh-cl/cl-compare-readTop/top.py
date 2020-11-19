# -*- coding: utf-8 -*-
"""
Created on Sun Apr  8 13:17:58 2018

@author: HuangMing
"""

# =============================================================================
# Calculate position, using data from mol.getAtom().getPos()
# Others parameter using data from mol.getTop.getAtoms()
# =============================================================================
from copy import deepcopy

class TopInfo(object):
    def __init__(self):
        self.__defaults__ = []
        self.__atomTypes__ = []
        self.__moleculeType__ = []
        self.__atoms__ = []
        self.__bonds__ = []
        self.__pairs__ = []
        self.__angles__ = []
        self.__dihedrals__ = []
        self.__system__ = []
        self.__molecules__ = []

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)
        return result
    
    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        result.__dict__.update(deepcopy(self.__dict__))
        return result

    def setDefaults(self, lst):
        self.__defaults__ = lst
        
    def setAtomTypes(self, lst):
        self.__atomTypes__ = lst

    def setMoleculeTypes(self, lst):
        self.__moleculeType__ = lst
    
    def setAtoms(self, lst):
        self.__atoms__ = lst
        
    def setBonds(self, lst):
        self.__bonds__ = lst
        
    def setPairs(self, lst):
        self.__pairs__ = lst
        
    def setAngles(self, lst):
        self.__angles__ = lst
        
    def setDihedrals(self, lst):
        self.__dihedrals__ = lst
    
    def setSystem(self, lst):
        self.__system__ = lst
    
    def setMolecules(self, lst):
        self.__molecules__ = lst
    
    
    def getDefaults(self):
        return self.__defaults__
    
    def getAtomTypes(self):
        return self.__atomTypes__
    
    def getMoleculeTypes(self):
        return self.__moleculeType__
    
    def getAtoms(self):
        return self.__atoms__
    
    def getBonds(self):
        return self.__bonds__
    
    def getPairs(self):
        return self.__pairs__
    
    def getAngles(self):
        return self.__angles__
    
    def getDihedrals(self):
        return self.__dihedrals__
    
    def getSystem(self):
        return self.__system__
    
    def getMolecules(self):
        return self.__molecules__
    
# =============================================================================
#     append another top 
# =============================================================================
    def appendTop(self, top, systemInfo=True):
        defaults = top.getDefaults()
        atomTypes = top.getAtomTypes()
        molTypes = top.getMoleculeTypes()
        atoms = top.getAtoms()
        bonds = top.getBonds()
        pairs = top.getPairs()
        angles = top.getAngles()
        dihs = top.getDihedrals()
        system = top.getSystem()
        molecules = top.getMolecules()
        
        self.appendDefaults(defaults)
        self.appendAtomTypes(atomTypes)
        self.appendMoleculeType(molTypes)
        self.appendAtoms(atoms)
        self.appendBonds(bonds)
        self.appendPairs(pairs)
        self.appendAngles(angles)
        self.appendDihedrals(dihs)
        
        if systemInfo:
            self.setSystem(system)
            self.setMolecules(molecules)
        
    def appendDefaults(self, defaults):
        new = self.__defaults__.copy()
        for i in defaults:
            if i not in self.__defaults__:
                new.append(i)
        self.__defaults__ = new
        
    def appendAtomTypes(self, atomTypes):
        new = self.__atomTypes__.copy()
        for i in atomTypes:
            if i not in self.__atomTypes__:
                new.append(i)
        self.__atomTypes__ = new
        
    def appendMoleculeType(self, moleculeTypes):
        new = self.__moleculeType__.copy()
        for i in moleculeTypes:
            if i not in self.__moleculeType__:
                new.append(i)
        self.__moleculeType__ = new
        
    def appendAtoms(self, atoms):
        for a in atoms:
            self.__atoms__.append(a)
    
    def appendBonds(self, bonds):
        for b in bonds:
            self.__bonds__.append(b)
        
    def appendPairs(self, pairs):
        for p in pairs:
            self.__pairs__.append(p)
    
    def appendAngles(self, angles):
        for a in angles:
            self.__angles__.append(a)
    
    def appendDihedrals(self, dihs):
        for d in dihs:
            self.__dihedrals__.append(d)
        
# =============================================================================
#     Top function
# =============================================================================
    def findAtoms(self, idx, atomsList):
        for atom in atomsList:
            idx1 = atom.getlocalIndex()
            if type(idx) != str and type(idx) != int:
                idx = idx.getlocalIndex()
            if int(idx) == int(idx1):
                return atom
            
    def idx2Atom(self):
        atoms = self.__atoms__
        #bonds
        for b in self.__bonds__:
            idx1 = b[0]
            idx2 = b[1]
#            print(idx1, idx2)
            a1 = self.findAtoms(idx1, atoms)
            a2 = self.findAtoms(idx2, atoms)
            b[0] = a1
            b[1] = a2
        
        #pairs
        for p in self.__pairs__:
            idx1 = p[0]
            idx2 = p[1]
            a1 = self.findAtoms(idx1, atoms)
            a2 = self.findAtoms(idx2, atoms)
            p[0] = a1
            p[1] = a2
        
        #angles
        for a in self.__angles__:
            idx1 = a[0]
            idx2 = a[1]
            idx3 = a[2]
            a1 = self.findAtoms(idx1, atoms)
            a2 = self.findAtoms(idx2, atoms)
            a3 = self.findAtoms(idx3, atoms)
            a[0] = a1
            a[1] = a2
            a[2] = a3
        
        #dihedrals
        for d in self.__dihedrals__:
            idx1 = d[0]
            idx2 = d[1]
            idx3 = d[2]
            idx4 = d[3]
#            print(idx1, idx2, idx3, idx4)
            a1 = self.findAtoms(idx1, atoms)
            a2 = self.findAtoms(idx2, atoms)
            a3 = self.findAtoms(idx3, atoms)
            a4 = self.findAtoms(idx4, atoms)
            d[0] = a1
            d[1] = a2
            d[2] = a3
            d[3] = a4
                
    def lst2Str(self, lst):
        str1 = ''
        for i in lst:
            if i == '':
                continue
            str1 += i + '\n'
        return str1
    
    def atomsInfo(self):
        lst = ['[ atoms ]']
        for a in self.__atoms__:
            lst.append(a.outputInfo(gro=False, top=True))
        
        return lst
    
    def bondsInfo(self):
        lst = ['[ bonds ]']
        for b in self.__bonds__:
            idx1 = b[0].getGlobalIndex()
            idx2 = b[1].getGlobalIndex()
            str1 = '{:>7}{:>7}{}'.format(idx1, idx2, b[2])
            lst.append(str1)
            
        return lst
    
    def pairsInfo(self):
        lst = ['[ pairs ]']
        for p in self.__pairs__:
            idx1 = p[0].getGlobalIndex()
            idx2 = p[1].getGlobalIndex()
            str1 = '{:>7}{:>7}{:>6}'.format(idx1, idx2, p[2])
            lst.append(str1)
        
        return lst
    
    def anglesInfo(self):
        lst = ['[ angles ]']
        for a in self.__angles__:
            idx1 = a[0].getGlobalIndex()
            idx2 = a[1].getGlobalIndex()
            idx3 = a[2].getGlobalIndex()
            str1 = '{:>7}{:>7}{:>7}{}'.format(idx1, idx2, idx3, a[3])
            lst.append(str1)
        
        return lst
    
    def dihsInfo(self):
        lst = ['[ dihedrals ]']
        for d in self.__dihedrals__:
            idx1 = d[0].getGlobalIndex()
            idx2 = d[1].getGlobalIndex()
            idx3 = d[2].getGlobalIndex()
            idx4 = d[3].getGlobalIndex()
            str1 = '{:>7}{:>7}{:>7}{:>7}{}'.format(idx1, idx2, idx3, idx4, d[4])
            lst.append(str1)
        
        return lst
    def outputInfo(self, outputName):
        defaults = self.__defaults__
        atomType = self.__atomTypes__
        molType = self.__moleculeType__
        atomInfo = self.atomsInfo()
        bondInfo = self.bondsInfo()
        pairInfo = self.pairsInfo()
        angleInfo = self.anglesInfo()
        dihsInfo = self.dihsInfo()
        sysInfo = self.__system__
        molecules = self.__molecules__
        
        str1 = self.lst2Str(defaults)
        str2 = self.lst2Str(atomType)
        str3 = self.lst2Str(molType)
        str4 = self.lst2Str(atomInfo)
        str5 = self.lst2Str(bondInfo)
        str6 = self.lst2Str(pairInfo)
        str7 = self.lst2Str(angleInfo)
        str8 = self.lst2Str(dihsInfo)
        str9 = self.lst2Str(sysInfo)
        str10 = self.lst2Str(molecules)
        f = open(outputName, 'w')
        lst = [str1, str2, str3, str4, str5, str6, str7, str8, str9, str10]
        f.write(self.lst2Str(lst))
        f.close()
        
        
        