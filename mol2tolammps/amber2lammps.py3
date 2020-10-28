#!/usr/bin/python
# -*- coding: utf-8 -*-

# This is amber2lammps, a program written by Keir E. Novik to convert
# Amber files to Lammps files.
# Copyright 1999, 2000 Keir E. Novik; all rights reserved.
#
# Modified by Vikas Varshney, U Akron, 5 July 2005, as described in README
# Bug Fixed: Third argument in Dihedral Coeffs section is an integer - Ketan S Khare September 26, 2011
# Modified by Vikas Varshney, Oct 8, 2013 to include additional flags (Atomic_Number, Coulombic and van der Waals 1-4 factors which are included in newer vesions of .top and .crd files in amber12.
# Modified by Ketan S. Khare - September 15, 2016 - (1) Improve Roundoff for Coeffs and Charges especially with Python 2.7 and later 
#                (2) Fix Improper atom ordering.  If I is central atom of improper I,J,K,L, Amber order is J K I -L. Unless the 3rd or 4th atom in dihedral is the first atom in file (i.e. 0).  LAMMPS expected order is I, J, K, L
# Modified by Ketan S. Khare - December 10, 2019 Converted to Python 3. (1) Masses lines have a comment followed by the AMBER Atom Type (2) Impropers separated

# PYTHON 3+
# ============================================================

def Pop(S, I=-1):
    '''Pop item I from list'''

    X = S[I]
    del S[I]
    return X

# ============================================================

class Lammps:

    # --------------------------------------------------------

    def Dump(self):
        '''Write out contents of self (intended for debugging)'''

        Name_list = list(self.__dict__.keys())
        Name_list.sort()
        for Name in Name_list:
            print(Name + ':', self.__dict__[Name])

    # --------------------------------------------------------

    def Write_data(self, Basename, Item_list):
        '''Write the Lammps data to file (used by Write_Lammps)'''

        import os
        import sys

        Filename = 'Data.' + Basename

        Dir_list = os.listdir('.')
        i = 1
        while Filename in Dir_list:
            Filename = 'Data' + repr(i) + '.' + Basename
            i = i + 1
        del i

        print('Writing', Filename + '...', end=' ')
        sys.stdout.flush()

        try:
            F = open(Filename, 'w')
        except IOError as Detail:
            print('(error:', Detail[1] + '!)')
            return

        try:
            F.writelines(Item_list)
        except IOError as Detail:
            print('(error:', Detail[1] + '!)')
            F.close()
            return

        F.close()
        print('done.')

    # --------------------------------------------------------

    def Write_Lammps(self, Basename):
        '''Write the Lammps data file, ignoring blank sections'''

        import string
        L = []

        L.append('LAMMPS Data File for ' + self.name + '\n\n')

        L.append(repr(self.nAtoms) + ' atoms\n')
        L.append(repr(self.nBonds) + ' bonds\n')
        L.append(repr(self.nAngles) + ' angles\n')
        L.append(repr(self.nDihedrals) + ' dihedrals\n')
        L.append(repr(self.nImpropers) + ' impropers\n\n')

        L.append(repr(self.nAtom_Types) + ' atom types\n')
        if self.nBonds > 0:
            L.append(repr(self.nBond_Types) + ' bond types\n')
        if self.nAngles > 0:
            L.append(repr(self.nAngle_Types) + ' angle types\n')
        if self.nDihedrals > 0:
            L.append(repr(self.nDihedral_Types) + ' dihedral types\n')
        if self.nImpropers > 0:
            L.append(repr(self.nImproper_Types) + ' improper types\n')
        L.append('\n')

        L.append(repr(self.xlo) + ' ' + repr(self.xhi) + ' xlo xhi\n')
        L.append(repr(self.ylo) + ' ' + repr(self.yhi) + ' ylo yhi\n')
        L.append(repr(self.zlo) + ' ' + repr(self.zhi) + ' zlo zhi\n\n')

        if self.nAtom_Types != 0:
            L.append('Masses\n\n')
            for i in range(self.nAtom_Types):
                L.append(repr(i + 1) + ' ' + repr(self.Masses[i]) + ' # ' + self.AmberAtomTypes[i] + '\n')
            L.append('\n')

            L.append('Pair Coeffs # lj/cut/coul/long\n\n')
            for i in range(self.nAtom_Types):
                L.append(repr(i + 1))
                for j in range(len(self.Pair_Coeffs[0])):
                    L.append(' ' + repr(self.Pair_Coeffs[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nBonds != 0 and self.nBond_Types != 0:
            L.append('Bond Coeffs # harmonic\n\n')
            for i in range(self.nBond_Types):
                L.append(repr(i + 1))
                for j in range(len(self.Bond_Coeffs[0])):
                    L.append(' ' + repr(self.Bond_Coeffs[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nAngles != 0 and self.nAngle_Types != 0:
            L.append('Angle Coeffs # harmonic \n\n')
            for i in range(self.nAngle_Types):
                L.append(repr(i + 1))
                for j in range(len(self.Angle_Coeffs[0])):
                    L.append(' ' + repr(self.Angle_Coeffs[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nDihedrals != 0 and self.nDihedral_Types != 0:
            L.append('Dihedral Coeffs # charmm\n\n')
            for i in range(self.nDihedral_Types):
                L.append(repr(i + 1))
                for j in range(len(self.Dihedral_Coeffs[0])):
                    L.append(' ' + repr(self.Dihedral_Coeffs[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nImpropers != 0 and self.nImproper_Types != 0:
            L.append('Improper Coeffs # cvff\n\n')
            for i in range(self.nImproper_Types):
                L.append(repr(i + 1))
                for j in range(len(self.Improper_Coeffs[0])):
                    L.append(' ' + repr(self.Improper_Coeffs[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nAtoms != 0:
            L.append('Atoms # full\n\n')
            for i in range(self.nAtoms):
                L.append(repr(i + 1))
                for j in range(len(self.Atoms[0])):
                    L.append(' ' + repr(self.Atoms[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nBonds != 0 and self.nBond_Types != 0:
            L.append('Bonds\n\n')
            for i in range(self.nBonds):
                L.append(repr(i + 1))
                for j in range(len(self.Bonds[0])):
                    L.append(' ' + repr(self.Bonds[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nAngles != 0 and self.nAngle_Types != 0:
            L.append('Angles\n\n')
            for i in range(self.nAngles):
                L.append(repr(i + 1))
                for j in range(len(self.Angles[0])):
                    L.append(' ' + repr(self.Angles[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nDihedrals != 0 and self.nDihedral_Types != 0:
            L.append('Dihedrals\n\n')
            for i in range(self.nDihedrals):
                L.append(repr(i + 1))
                for j in range(len(self.Dihedrals[0])):
                    L.append(' ' + repr(self.Dihedrals[i][j]))
                L.append('\n')
            L.append('\n')

        if self.nImpropers != 0 and self.nImproper_Types != 0:
            L.append('Impropers\n\n')
            for i in range(self.nImpropers):
                L.append(repr(i + 1))
                for j in range(len(self.Impropers[0])):
                    L.append(' ' + repr(self.Impropers[i][j]))
                L.append('\n')
            L.append('\n')

        self.Write_data(Basename, L)


# ============================================================

class Amber:

    def __init__(self):
        '''Initialise the Amber class'''

        self.CRD_is_read = 0
        self.TOP_is_read = 0

    # --------------------------------------------------------

    def Dump(self):
        '''Write out contents of self (intended for debugging)'''

        Name_list = list(self.__dict__.keys())
        Name_list.sort()
        for Name in Name_list:
            print(Name + ':', self.__dict__[Name])

    # --------------------------------------------------------

    def Coerce_to_Lammps(self):
        '''Return the Amber data converted to Lammps format'''

        import math

        if self.CRD_is_read and self.TOP_is_read:
            l = Lammps()
            print('Converting...', end=' ')

            l.name = self.ITITL
            l.nAtoms = self.NATOM
            l.nBonds = self.NBONH + self.MBONA
            l.nAngles = self.NTHETH + self.MTHETA
            l.nDihedrals = self.NPHIH + self.MPHIA
            l.nImpropers = 0                            #AMBER merges Dihedrals and Impropers. Will be Separated Later.
            l.nAtom_Types = self.NTYPES
            l.nBond_Types = self.NUMBND
            l.nAngle_Types = self.NUMANG
            l.nDihedral_Types = self.NPTRA
            l.nImproper_Types = self.NPTRA

            Shift = 0
            if 'BOX' in self.__dict__:
                l.xlo = 0.0
                l.xhi = self.BOX[0]
                l.ylo = 0.0
                l.yhi = self.BOX[1]
                l.zlo = 0.0
                l.zhi = self.BOX[2]
                if l.xlo > min(self.X) or l.xhi < max(self.X) or l.ylo \
                    > min(self.Y) or l.yhi < max(self.Y) or l.zlo \
                    > min(self.Z) or l.zhi < max(self.Z):

                    #  Vikas Modification: Disabling Shifting. This means I am intend to send exact coordinates of each atom and let LAMMPS
                    #  take care of imaging into periodic image cells. If one wants to shift all atoms in the periodic box,
                    #  please uncomment the below 2 lines.

                    print('(warning: Currently not shifting the atoms to the periodic box)')
            else:                              
                #Adding a padding of 1.75 A to Box
                print('Warning: Guessing box size from crd. With pad of 1.75 A!', end=' ')
                l.xlo = min(self.X) - 1.75
                l.xhi = max(self.X) + 1.75
                l.ylo = min(self.Y) - 1.75
                l.yhi = max(self.Y) + 1.75
                l.zlo = min(self.Z) - 1.75
                l.zhi = max(self.Z) + 1.75

            # This doesn't check duplicate values

            l.Masses = []
            l.AmberAtomTypes = []
            for i in range(l.nAtom_Types):
                l.Masses.append(0)
                l.AmberAtomTypes.append(0)
            for i in range(self.NATOM):
                l.Masses[self.IAC[i] - 1] = self.AMASS[i]
                l.AmberAtomTypes[self.IAC[i] - 1] = self.ISYMBL[i]

            l.Pair_Coeffs = []
            for i in range(self.NTYPES):
                l.Pair_Coeffs.append([0, 0])
            for i in range(self.NTYPES):
                j = self.ICO[i * (self.NTYPES + 1)] - 1
                if self.CN1[j] == 0.0:
                    l.Pair_Coeffs[i][0] = 0.0
                else:
                    l.Pair_Coeffs[i][0] = round(0.25 * self.CN2[j] ** 2 / self.CN1[j], 4)
                if self.CN2[j] == 0.0:
                    l.Pair_Coeffs[i][1] = 0.0
                else:
                    l.Pair_Coeffs[i][1] = round((self.CN1[j] / self.CN2[j]) ** (1.0 / 6.0), 3)

            l.Bond_Coeffs = []
            for i in range(self.NUMBND):
                l.Bond_Coeffs.append([0, 0])
            for i in range(self.NUMBND):
                l.Bond_Coeffs[i][0] = self.RK[i]
                l.Bond_Coeffs[i][1] = self.REQ[i]

            l.Angle_Coeffs = []
            for i in range(self.NUMANG):
                l.Angle_Coeffs.append([0, 0])
            for i in range(self.NUMANG):
                l.Angle_Coeffs[i][0] = self.TK[i]
                l.Angle_Coeffs[i][1] = round(180 / math.pi * self.TEQ[i], 3)

            l.Dihedral_Coeffs = []
            for i in range(self.NPTRA):
                l.Dihedral_Coeffs.append([0, 0, 0, 0])
            for i in range(self.NPTRA):
                l.Dihedral_Coeffs[i][0] = round(self.PK[i],3)
                l.Dihedral_Coeffs[i][1] = int(self.PN[i])
                l.Dihedral_Coeffs[i][2] = int(180 / math.pi * self.PHASE[i])

            l.Improper_Coeffs = []  # Improper Coeffs section is a copy of Dihedral Coeffs. Will be corrected later.
            for i in range(self.NPTRA):
                l.Improper_Coeffs.append([0, 0, 0])
            for i in range(self.NPTRA):
                l.Improper_Coeffs[i][0] = round(self.PK[i],3)
                if self.PHASE[i] == 0:
                    l.Improper_Coeffs[i][1] = 1
                else:
                    l.Improper_Coeffs[i][1] = -1
                l.Improper_Coeffs[i][2] = int(self.PN[i])

            l.Atoms = []
            NetCharge = 0
            for i in range(self.NATOM):
                x = self.X[i]
                y = self.Y[i]
                z = self.Z[i]
                if Shift:
                    while x < l.xlo:
                        x = x + self.BOX[0]
                    while x > l.xhi:
                        x = x - self.BOX[0]
                    while y < l.ylo:
                        y = y + self.BOX[1]
                    while y > l.yhi:
                        y = y - self.BOX[1]
                    while z < l.zlo:
                        z = z + self.BOX[2]
                    while z > l.zhi:
                        z = z - self.BOX[2]
                l.Atoms.append([1,self.IAC[i],round(self.CHRG[i] / 18.2223, 6),x,y,z,])
                NetCharge = NetCharge + round(self.CHRG[i] / 18.2223, 6)

            l.Bonds = []
            for i in range(l.nBonds):
                l.Bonds.append([0, 0, 0])
            for i in range(self.NBONH):
                l.Bonds[i][0] = self.ICBH[i]
                l.Bonds[i][1] = abs(self.IBH[i]) // 3 + 1
                l.Bonds[i][2] = abs(self.JBH[i]) // 3 + 1
            for i in range(self.NBONA):
                l.Bonds[self.NBONH + i][0] = self.ICB[i]
                l.Bonds[self.NBONH + i][1] = abs(self.IB[i]) // 3 + 1
                l.Bonds[self.NBONH + i][2] = abs(self.JB[i]) // 3 + 1

            l.Angles = []
            for i in range(l.nAngles):
                l.Angles.append([0, 0, 0, 0])
            for i in range(self.NTHETH):
                l.Angles[i][0] = self.ICTH[i]
                l.Angles[i][1] = abs(self.ITH[i]) // 3 + 1
                l.Angles[i][2] = abs(self.JTH[i]) // 3 + 1
                l.Angles[i][3] = abs(self.KTH[i]) // 3 + 1
            for i in range(self.NTHETA):
                l.Angles[self.NTHETH + i][0] = self.ICT[i]
                l.Angles[self.NTHETH + i][1] = abs(self.IT[i]) // 3 + 1
                l.Angles[self.NTHETH + i][2] = abs(self.JT[i]) // 3 + 1
                l.Angles[self.NTHETH + i][3] = abs(self.KT[i]) // 3 + 1

            #Impropers are distinguised from dihedrals by making fourth atom index as negative
            for i in range(self.NPHIH):
                if self.LPH[i] < 0:
                    l.nImpropers += 1
            for i in range(self.NPHIA):
                if self.LP[i] < 0:
                    l.nImpropers += 1

            l.nDihedrals -= l.nImpropers

            l.Dihedrals = []
            for i in range(l.nDihedrals):
                l.Dihedrals.append([0, 0, 0, 0, 0])
            j = 0
            for i in range(self.NPHIH):
                if self.LPH[i] > 0:
                    l.Dihedrals[j][0] = self.ICPH[i]
                    l.Dihedrals[j][1] = abs(self.IPH[i]) // 3 + 1
                    l.Dihedrals[j][2] = abs(self.JPH[i]) // 3 + 1
                    l.Dihedrals[j][3] = abs(self.KPH[i]) // 3 + 1
                    l.Dihedrals[j][4] = abs(self.LPH[i]) // 3 + 1
                    j += 1
            for i in range(self.NPHIA):
                if self.LP[i] > 0:
                    l.Dihedrals[j][0] = self.ICP[i]
                    l.Dihedrals[j][1] = abs(self.IP[i]) // 3 + 1
                    l.Dihedrals[j][2] = abs(self.JP[i]) // 3 + 1
                    l.Dihedrals[j][3] = abs(self.KP[i]) // 3 + 1
                    l.Dihedrals[j][4] = abs(self.LP[i]) // 3 + 1
                    j += 1
            
            #Note difference in index order for impropers in AMBER versus LAMMPS
            l.Impropers = []
            for i in range(l.nImpropers):
                l.Impropers.append([0, 0, 0, 0, 0])
            j = 0
            for i in range(self.NPHIH):
                if self.LPH[i] < 0:
                    if self.JPH[i] != 0:
                        l.Impropers[j][0] = self.ICPH[i]
                        l.Impropers[j][1] = abs(self.KPH[i]) // 3 + 1
                        l.Impropers[j][2] = abs(self.IPH[i]) // 3 + 1
                        l.Impropers[j][3] = abs(self.JPH[i]) // 3 + 1
                        l.Impropers[j][4] = abs(self.LPH[i]) // 3 + 1
                    else:
                          # AMBER format will reverse order if 4th atom in dihedral is index 0, since -0 is 0
                        l.Impropers[j][0] = self.ICPH[i]
                        l.Impropers[j][1] = abs(self.JPH[i]) // 3 + 1
                        l.Impropers[j][2] = abs(self.LPH[i]) // 3 + 1
                        l.Impropers[j][3] = abs(self.KPH[i]) // 3 + 1
                        l.Impropers[j][4] = abs(self.IPH[i]) // 3 + 1
                    j += 1

            for i in range(self.NPHIA):
                if self.LP[i] < 0:
                    if self.JP[i] != 0:
                        l.Impropers[j][0] = self.ICP[i]
                        l.Impropers[j][1] = abs(self.KP[i]) // 3 + 1
                        l.Impropers[j][2] = abs(self.IP[i]) // 3 + 1
                        l.Impropers[j][3] = abs(self.JP[i]) // 3 + 1
                        l.Impropers[j][4] = abs(self.LP[i]) // 3 + 1
                    else:
                          # AMBER format reverse order if 4th atom in dihedral is index 0, since -0 is 0
                        l.Impropers[j][0] = self.ICP[i]
                        l.Impropers[j][1] = abs(self.JP[i]) // 3 + 1
                        l.Impropers[j][2] = abs(self.LP[i]) // 3 + 1
                        l.Impropers[j][3] = abs(self.KP[i]) // 3 + 1
                        l.Impropers[j][4] = abs(self.IP[i]) // 3 + 1
                    j += 1
                    
            #Remove unused dihedral and improper types from the Coeffs sections
            UseFlag=[]  
            for i in range(l.nDihedral_Types):
                UseFlag.append(0)
                for j in range(l.nDihedrals):
                    if(l.Dihedrals[j][0] == i+1):
                        UseFlag[i]=1
                        break
            
            i=0
            while(i<l.nDihedral_Types):
                if(UseFlag[i]==0):
                    del l.Dihedral_Coeffs[i]
                    del UseFlag[i]
                    l.nDihedral_Types-=1
                    for j in range(l.nDihedrals):
                        if(l.Dihedrals[j][0] > i):
                            l.Dihedrals[j][0]-=1
                else: 
                    i+=1                
                                      
            UseFlag=[]
            for i in range(l.nImproper_Types):
                UseFlag.append(0)
                for j in range(l.nImpropers):
                    if(l.Impropers[j][0] == i+1):
                        UseFlag[i]=1
                        break
            i=0
            while(i<l.nImproper_Types):
                if(UseFlag[i]==0):
                    del l.Improper_Coeffs[i]
                    del UseFlag[i]
                    l.nImproper_Types-=1
                    for j in range(l.nImpropers):
                        if(l.Impropers[j][0] > i):
                            l.Impropers[j][0]-=1
                else: 
                    i+=1   
            
            if(NetCharge!=0):
                print('Warning: Net Charge is ', round(NetCharge,6), '\n')
                
            print('done.')
            return l
        else:
            print('(Error: Not all the Amber data has been read!)')

    # --------------------------------------------------------

    def Read_data(self, Filename):
        '''Read the filename, returning a list of strings'''

        import string
        import sys

        print('Reading', Filename + '...', end=' ')
        sys.stdout.flush()

        try:
            F = open(Filename)
        except IOError as Detail:
            print('(error:', Detail[1] + '!)')
            return

        try:
            Lines = F.readlines()
        except IOError as Detail:
            print('(error:', Detail[1] + '!)')
            F.close()
            return

        F.close()

        # If the first line is empty, use the Basename

        if Filename[-4:] == '.crd':
            if (not Lines[0]):  # This line corresponds to TITLE name in CRD file
                Basename = Filename[:string.find(Filename, '.')]
                Item_list = [Basename]
                print('Warning: Title not present... Assigning Basename as Title')
            else:
                Item_list = []
        else:
            if (not Lines[3]):  # This line corresponds to TITLE name in TOPOLOGY file
                Basename = Filename[:string.find(Filename, '.')]
                Item_list = [Basename]
                print('Warning: Title not present... Assigning Basename as Title')
            else:
                Item_list = []

        for Line in Lines:
            if Line[0] != '%':  # Vikas' Modification: This condition ignores all the lines starting with % in the topology file.
                Item_list.extend(Line.split())

        return Item_list

    # --------------------------------------------------------

    def Read_CRD(self, Basename):
        '''Read the Amber coordinate/restart (.crd) file'''

        # The optional velocities and periodic box size are not yet parsed.

        Item_list = self.Read_data(Basename + '.crd')

        if Item_list == None:
            return
        elif len(Item_list) < 2:
            print('(error: File too short!)')
            return

        # Parse the data

        if 'ITITL' in self.__dict__:
            if Pop(Item_list, 0) != self.ITITL:
                print('(warning: ITITL differs!)', end=' ')
        else:
            self.ITITL = Pop(Item_list, 0)
        print("CRD Title is " + self.ITITL)  # Vikas Modification : Priting the Title

        if 'NATOM' in self.__dict__:
            if eval(Pop(Item_list, 0)) != self.NATOM:
                print('(error: NATOM differs!)')
                return
        else:
            self.NATOM = eval(Pop(Item_list, 0))
        print(repr(self.NATOM) + " atoms found!", end = " ")  # Vikas' Modification: Printing number of atoms just to make sure that the program is reading the correct value.

        # if len(Item_list) == 1 + 3 * self.NATOM:
        # Vikas' Modification: I changed the condition.

        if len(Item_list) % 3 != 0:
            self.TIME = eval(Pop(Item_list, 0))
        else:
            self.TIME = 0
        print("CRD Simulation time is: " + repr(self.TIME))  # Vikas' Modification : Printing simulation time, just to make sure that the program is readint the correct value.
        if len(Item_list) < 3 * self.NATOM:
            print('(error: File too short!)')
            return

        self.X = []
        self.Y = []
        self.Z = []
        for i in range(self.NATOM):
            self.X.append(eval(Pop(Item_list, 0)))
            self.Y.append(eval(Pop(Item_list, 0)))
            self.Z.append(eval(Pop(Item_list, 0)))

        if self.NATOM == 1 and len(Item_list):
            print('(warning: Ambiguity!)', end=' ')

        if len(Item_list) >= 3 * self.NATOM:
            self.VX = []
            self.VY = []
            self.VZ = []
            for i in range(self.NATOM):
                self.VX.append(eval(Pop(Item_list, 0)))
                self.VY.append(eval(Pop(Item_list, 0)))
                self.VZ.append(eval(Pop(Item_list, 0)))

        if len(Item_list) >= 3:
            self.BOX = []
            for i in range(3):
                self.BOX.append(eval(Pop(Item_list, 0)))

        if len(Item_list):
            print('(warning: File too large!)', end=' ')

        print('CRD file is done.')
        self.CRD_is_read = 1

    # --------------------------------------------------------

    def Read_TOP(self, Basename):
        '''Read the Amber parameter/topology (.top) file'''

        Item_list = self.Read_data(Basename + '.top')

        if Item_list == None:
            return
        elif len(Item_list) < 31:
            print('(error: File too short!)')
            return

       # Parse the data

        if 'ITITL' in self.__dict__:
            if Pop(Item_list, 0) != self.ITITL:
                print('(warning: ITITL differs!)')
        else:
            self.ITITL = Pop(Item_list, 0)
        print("TOP Title is " + self.ITITL)  # Printing Self Title

        if 'NATOM' in self.__dict__:
            if eval(Pop(Item_list, 0)) != self.NATOM:
                print('(error: NATOM differs!)')
                return
        else:
            self.NATOM = eval(Pop(Item_list, 0))
        print(repr(self.NATOM) + " atoms found!")  # Printing total number of atoms just to make sure that thing are going right
        self.NTYPES = eval(Pop(Item_list, 0))
        self.NBONH = eval(Pop(Item_list, 0))
        self.MBONA = eval(Pop(Item_list, 0))
        self.NTHETH = eval(Pop(Item_list, 0))
        self.MTHETA = eval(Pop(Item_list, 0))
        self.NPHIH = eval(Pop(Item_list, 0))
        self.MPHIA = eval(Pop(Item_list, 0))
        self.NHPARM = eval(Pop(Item_list, 0))
        self.NPARM = eval(Pop(Item_list, 0))
        self.NEXT = eval(Pop(Item_list, 0))
        self.NRES = eval(Pop(Item_list, 0))
        self.NBONA = eval(Pop(Item_list, 0))
        self.NTHETA = eval(Pop(Item_list, 0))
        self.NPHIA = eval(Pop(Item_list, 0))
        self.NUMBND = eval(Pop(Item_list, 0))
        self.NUMANG = eval(Pop(Item_list, 0))
        self.NPTRA = eval(Pop(Item_list, 0))
        self.NATYP = eval(Pop(Item_list, 0))
        self.NPHB = eval(Pop(Item_list, 0))
        self.IFPERT = eval(Pop(Item_list, 0))
        self.NBPER = eval(Pop(Item_list, 0))
        self.NGPER = eval(Pop(Item_list, 0))
        self.NDPER = eval(Pop(Item_list, 0))
        self.MBPER = eval(Pop(Item_list, 0))
        self.MGPER = eval(Pop(Item_list, 0))
        self.MDPER = eval(Pop(Item_list, 0))
        self.IFBOX = eval(Pop(Item_list, 0))
        self.NMXRS = eval(Pop(Item_list, 0))
        self.IFCAP = eval(Pop(Item_list, 0))

        # ....................................................

        if len(Item_list) < 5 * self.NATOM + self.NTYPES ** 2 + 2 \
            * (self.NRES + self.NUMBND + self.NUMANG) + 3 * self.NPTRA \
            + self.NATYP:
            print('(error: File too short!)')
            return -1

        self.IGRAPH = []
        Pop(Item_list, 0)

      # A little kludge is needed here, since the IGRAPH strings are
        # not separated by spaces if 4 characters in length.

        for i in range(self.NATOM):
            if len(Item_list[0]) > 4:
                Item_list.insert(1, (Item_list[0])[4:])
                Item_list.insert(1, (Item_list[0])[0:4])
                del Item_list[0]
            self.IGRAPH.append(Pop(Item_list, 0))

    # Vikas' Modification : In the following section, I am printing out each quantity which is currently being read from the topology file.

        print('Reading Charges...')
        self.CHRG = []
        for i in range(self.NATOM):
            self.CHRG.append(eval(Pop(Item_list, 0)))

        print('Reading Atomic Number...')
        self.ANUMBER = []
        for i in range(self.NATOM):
            self.ANUMBER.append(eval(Pop(Item_list, 0)))

        print('Reading Atomic Masses...')
        self.AMASS = []
        for i in range(self.NATOM):
            self.AMASS.append(eval(Pop(Item_list, 0)))

        print('Reading Atom Types...')
        self.IAC = []
        for i in range(self.NATOM):
            self.IAC.append(eval(Pop(Item_list, 0)))

        print('Reading Excluded Atoms...')
        self.NUMEX = []
        for i in range(self.NATOM):
            self.NUMEX.append(eval(Pop(Item_list, 0)))

        print('Reading Non-bonded Parameter Index...')
        self.ICO = []
        for i in range(self.NTYPES ** 2):
            self.ICO.append(eval(Pop(Item_list, 0)))

        print('Reading Residue Labels...')
        self.LABRES = []
        for i in range(self.NRES):
            self.LABRES.append(Pop(Item_list, 0))

        print('Reading Residues Starting Pointers...')
        self.IPRES = []
        for i in range(self.NRES):
            self.IPRES.append(eval(Pop(Item_list, 0)))

        print('Reading Bond Force Constants...')
        self.RK = []
        for i in range(self.NUMBND):
            self.RK.append(eval(Pop(Item_list, 0)))

        print('Reading Equilibrium Bond Values...')
        self.REQ = []
        for i in range(self.NUMBND):
            self.REQ.append(eval(Pop(Item_list, 0)))

        print('Reading Angle Force Constants...')
        self.TK = []
        for i in range(self.NUMANG):
            self.TK.append(eval(Pop(Item_list, 0)))

        print('Reading Equilibrium Angle Values...')
        self.TEQ = []
        for i in range(self.NUMANG):
            self.TEQ.append(eval(Pop(Item_list, 0)))

        print('Reading Dihedral Force Constants...')
        self.PK = []
        for i in range(self.NPTRA):
            self.PK.append(eval(Pop(Item_list, 0)))

        print('Reading Dihedral Periodicity...')
        self.PN = []
        for i in range(self.NPTRA):
            self.PN.append(eval(Pop(Item_list, 0)))

        print('Reading Dihedral Phase...')
        self.PHASE = []
        for i in range(self.NPTRA):
            self.PHASE.append(eval(Pop(Item_list, 0)))

        print('Reading 1-4 Electrostatic Scaling Factor...')
        self.SCEEFAC = []
        for i in range(self.NPTRA):
            self.SCEEFAC.append(eval(Pop(Item_list, 0)))

        print('Reading 1-4 Van der Waals Scaling Factor...')
        self.SCNBFAC = []
        for i in range(self.NPTRA):
            self.SCNBFAC.append(eval(Pop(Item_list, 0)))

        print('Reading Solty...')
        self.SOLTY = []
        for i in range(self.NATYP):
            self.SOLTY.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if len(Item_list) < 2 * self.NTYPES * (self.NTYPES + 1) / 2:
            print('(error: File too short!)')
            return -1

        print('Reading LJ A Coefficient...')
        self.CN1 = []
        for i in range(self.NTYPES * (self.NTYPES + 1) // 2):
            self.CN1.append(eval(Pop(Item_list, 0)))

        print('Reading LJ B Coefficient...')
        self.CN2 = []
        for i in range(self.NTYPES * (self.NTYPES + 1) // 2):
            self.CN2.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if len(Item_list) < 3 * (self.NBONH + self.NBONA) + 4 \
            * (self.NTHETH + self.NTHETA) + 5 * (self.NPHIH
                + self.NPHIA):
            print('(error: File too short!)')
            return -1

        print('Reading Bonds which include hydrogen...')
        self.IBH = []
        self.JBH = []
        self.ICBH = []
        for i in range(self.NBONH):
            self.IBH.append(eval(Pop(Item_list, 0)))
            self.JBH.append(eval(Pop(Item_list, 0)))
            self.ICBH.append(eval(Pop(Item_list, 0)))

        print('Reading Bonds which dont include hydrogen...')
        self.IB = []
        self.JB = []
        self.ICB = []
        for i in range(self.NBONA):
            self.IB.append(eval(Pop(Item_list, 0)))
            self.JB.append(eval(Pop(Item_list, 0)))
            self.ICB.append(eval(Pop(Item_list, 0)))

        print('Reading Angles which include hydrogen...')
        self.ITH = []
        self.JTH = []
        self.KTH = []
        self.ICTH = []
        for i in range(self.NTHETH):
            self.ITH.append(eval(Pop(Item_list, 0)))
            self.JTH.append(eval(Pop(Item_list, 0)))
            self.KTH.append(eval(Pop(Item_list, 0)))
            self.ICTH.append(eval(Pop(Item_list, 0)))

        print('Reading Angles which dont include hydrogen...')
        self.IT = []
        self.JT = []
        self.KT = []
        self.ICT = []
        for i in range(self.NTHETA):
            self.IT.append(eval(Pop(Item_list, 0)))
            self.JT.append(eval(Pop(Item_list, 0)))
            self.KT.append(eval(Pop(Item_list, 0)))
            self.ICT.append(eval(Pop(Item_list, 0)))

        print('Reading Dihedrals which include hydrogen...')
        self.IPH = []
        self.JPH = []
        self.KPH = []
        self.LPH = []
        self.ICPH = []
        for i in range(self.NPHIH):
            self.IPH.append(eval(Pop(Item_list, 0)))
            self.JPH.append(eval(Pop(Item_list, 0)))
            self.KPH.append(eval(Pop(Item_list, 0)))
            self.LPH.append(eval(Pop(Item_list, 0)))
            self.ICPH.append(eval(Pop(Item_list, 0)))

        print('Reading Dihedrals which dont include hydrogen...')
        self.IP = []
        self.JP = []
        self.KP = []
        self.LP = []
        self.ICP = []
        for i in range(self.NPHIA):
            self.IP.append(eval(Pop(Item_list, 0)))
            self.JP.append(eval(Pop(Item_list, 0)))
            self.KP.append(eval(Pop(Item_list, 0)))
            self.LP.append(eval(Pop(Item_list, 0)))
            self.ICP.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if len(Item_list) < self.NEXT + 3 * self.NPHB + 4 * self.NATOM:
            print('(error: File too short!)')
            return -1

        print('Reading Excluded Atom List...')
        self.NATEX = []
        for i in range(self.NEXT):
            self.NATEX.append(eval(Pop(Item_list, 0)))

        print('Reading H-Bond A Coefficient, corresponding to r**12 term for all possible types...')
        self.ASOL = []
        for i in range(self.NPHB):
            self.ASOL.append(eval(Pop(Item_list, 0)))

        print('Reading H-Bond B Coefficient, corresponding to r**10 term for all possible types...')
        self.BSOL = []
        for i in range(self.NPHB):
            self.BSOL.append(eval(Pop(Item_list, 0)))

        print('Reading H-Bond Cut...')  # I think it is not being used nowadays
        self.HBCUT = []
        for i in range(self.NPHB):
            self.HBCUT.append(eval(Pop(Item_list, 0)))

        print('Reading Amber Atom Types for each atom...')
        self.ISYMBL = []
        for i in range(self.NATOM):
            self.ISYMBL.append(Pop(Item_list, 0))

        print('Reading Tree Chain Classification...')
        self.ITREE = []
        for i in range(self.NATOM):
            self.ITREE.append(Pop(Item_list, 0))

        print('Reading Join Array: Tree joining information')  # Currently unused in Sander, an AMBER module
        self.JOIN = []
        for i in range(self.NATOM):
            self.JOIN.append(eval(Pop(Item_list, 0)))

        print('Reading IRotate...')  # Currently unused in Sander and Gibbs
        self.IROTAT = []
        for i in range(self.NATOM):
            self.IROTAT.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if self.IFBOX > 0:
            if len(Item_list) < 3:
                print('(error: File too short!)')
                return -1

            print('Reading final residue which is part of solute...')
            self.IPTRES = eval(Pop(Item_list, 0))
            print('Reading total number of molecules...')
            self.NSPM = eval(Pop(Item_list, 0))
            print('Reading first solvent moleule index...')
            self.NSPSOL = eval(Pop(Item_list, 0))

            if len(Item_list) < self.NSPM + 4:
                print('(error: File too short!)')
                return -1

            print('Reading atom per molecule...')
            self.NSP = []
            for i in range(self.NSPM):
                self.NSP.append(eval(Pop(Item_list, 0)))

            self.BETA = eval(Pop(Item_list, 0))

            print('Reading Box Dimensions...')
            if 'BOX' in self.__dict__:
                BOX = []
                for i in range(3):
                    BOX.append(eval(Pop(Item_list, 0)))
                for i in range(3):
                    if BOX[i] != self.BOX[i]:
                        print('(warning: BOX differs!)', end=' ')
                        break
                del BOX
            else:
                self.BOX = []
                for i in range(3):
                    self.BOX.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if self.IFCAP > 0:
            if len(Item_list) < 5:
                print('(error: File too short!)')
                return -1
            print('Reading ICAP variables::: For details, refer to online AMBER format manual')
            self.NATCAP = eval(Pop(Item_list, 0))
            self.CUTCAP = eval(Pop(Item_list, 0))
            self.XCAP = eval(Pop(Item_list, 0))
            self.YCAP = eval(Pop(Item_list, 0))
            self.ZCAP = eval(Pop(Item_list, 0))

        # ....................................................

        if self.IFPERT > 0:
            if len(Item_list) < 4 * self.NBPER + 5 * self.NGPER + 6 \
                * self.NDPER + self.NRES + 6 * self.NATOM:
                print('(error: File too short!)')
                return -1

            print('Reading perturb variables, 1. Bond, 2. Angles, 3. Dihedrals, etc etc.::: For details, refer to online AMBER format manual')
            self.IBPER = []
            self.JBPER = []
            for i in range(self.NBPER):
                self.IBPER.append(eval(Pop(Item_list, 0)))
                self.JBPER.append(eval(Pop(Item_list, 0)))

            self.ICBPER = []
            for i in range(2 * self.NBPER):
                self.ICBPER.append(eval(Pop(Item_list, 0)))

            self.ITPER = []
            self.JTPER = []
            self.KTPER = []
            for i in range(self.NGPER):
                self.ITPER.append(eval(Pop(Item_list, 0)))
                self.JTPER.append(eval(Pop(Item_list, 0)))
                self.KTPER.append(eval(Pop(Item_list, 0)))

            self.ICTPER = []
            for i in range(2 * self.NGPER):
                self.ICTPER.append(eval(Pop(Item_list, 0)))

            self.IPPER = []
            self.JPPER = []
            self.KPPER = []
            self.LPPER = []
            for i in range(self.NDPER):
                self.IPPER.append(eval(Pop(Item_list, 0)))
                self.JPPER.append(eval(Pop(Item_list, 0)))
                self.KPPER.append(eval(Pop(Item_list, 0)))
                self.LPPER.append(eval(Pop(Item_list, 0)))

            self.ICPPER = []
            for i in range(2 * self.NDPER):
                self.ICPPER.append(eval(Pop(Item_list, 0)))

            LABRES = []
            for i in range(self.NRES):
                LABRES.append(Pop(Item_list, 0))
            for i in range(self.NRES):
                if LABRES[i] != self.LABRES[i]:
                    print('(warning: BOX differs!)', end=' ')
                    break

            self.IGRPER = []
            for i in range(self.NATOM):
                self.IGRPER.append(eval(Pop(Item_list, 0)))

            self.ISMPER = []
            for i in range(self.NATOM):
                self.ISMPER.append(eval(Pop(Item_list, 0)))

            self.ALMPER = []
            for i in range(self.NATOM):
                self.ALMPER.append(eval(Pop(Item_list, 0)))

            self.IAPER = []
            for i in range(self.NATOM):
                self.IAPER.append(eval(Pop(Item_list, 0)))

            self.IACPER = []
            for i in range(self.NATOM):
                self.IACPER.append(eval(Pop(Item_list, 0)))

            self.CGPER = []
            for i in range(self.NATOM):
                self.CGPER.append(eval(Pop(Item_list, 0)))

        # ....................................................

        
        
        
        self.IPOL = 0
        if self.IPOL == 1:
            if len(Item_list) < self.NATOM:
                print('(error: File too short!)')
                return -1
            print('Reading Polarizability Data. For details, refer to online AMBER format manual')
            self.ATPOL = []
            for i in range(self.NATOM):
                self.ATPOL.append(eval(Pop(Item_list, 0)))

            if self.IFPERT == 1:
                if len(Item_list) < self.NATOM:
                    print('(error: File too short!)')
                    return -1
                self.ATPOL1 = []
                for i in range(self.NATOM):
                    self.ATPOL1.append(eval(Pop(Item_list, 0)))

        # ....................................................

        if len(Item_list):
            print('(warning: AMBER Topology file is too large! Usually Safe to ignore!)')

        print("TOP file is done.")
        self.TOP_is_read = 1


# ============================================================

def Find_Amber_files():
    '''Look for sets of Amber files to process'''

    # Date and existence checks not yet implemented

    import os
    import sys

    Basename_list = []

    # Extract basenames from command line

    for Name in sys.argv[1:]:
        if Name[-4:] == '.crd':
            Basename_list.append(Name[:-4])
        else:
            if Name[-4:] == '.top':
                Basename_list.append(Name[:-4])
            else:
                Basename_list.append(Name)

    # Remove duplicate basenames

    for Basename in Basename_list[:]:
        while Basename_list.count(Basename) > 1:
            Basename_list.remove(Basename)

    if Basename_list == []:
        print('Looking for Amber files...', end=' ')
        Dir_list = os.listdir('.')
        Dir_list.sort()
        for File in Dir_list:
            if File[-4:] == '.top':
                Basename = File[:-4]
                if Basename + '.crd' in Dir_list:
                    Basename_list.append(Basename)
        if Basename_list != []:
            print('found', end=' ')
            for i in range(len(Basename_list) - 1):
                print(Basename_list[i] + ',', end=' ')
            print(Basename_list[-1] + '\n')

    if Basename_list == []:
        print('none.\n')

    return Basename_list

# ============================================================

def Check_Python_Version():
    '''Python Version'''
    import sys
    if sys.version_info[0] < 3:
        raise Exception("Error: Must be using Python 3")

# ============================================================

def Convert_Amber_files():
    '''Handle the whole conversion process'''
    print("\nWelcome to amber2lammps, a program to convert AMBER files to LAMMPS format!\n")
    Check_Python_Version()
    Basename_list = Find_Amber_files()
    for Basename in Basename_list:
        a = Amber()
        a.Read_CRD(Basename)
        if a.CRD_is_read:
            a.Read_TOP(Basename)
            if a.TOP_is_read:
                l = a.Coerce_to_Lammps()
                l.Write_Lammps(Basename)
                del l
        del a
        print("Program: amber2lammps has exited gracefully.\n")


# ============================================================

Convert_Amber_files()
