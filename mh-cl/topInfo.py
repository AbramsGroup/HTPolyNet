# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:34:47 2020

@author: huangming
"""
import pandas as pd
from countTime import *
from copy import deepcopy
import numpy as np
import decimal
decimal.getcontext().prec = 6

class top(object):
    def __init__(self):
        self.default = pd.DataFrame(['1               2               yes             0.5     0.8333'])
        self.atomtypes = []
        self.bondtypes = []
        self.angletypes = []
        self.dihtypes = []
        self.imptypes = []
        self.moleculetype = []
        self.atoms = []
        self.bonds = []
        self.pairs = []
        self.angles = []
        self.dihedrals = []
        self.impropers = []
        self.system = []
        self.molecules = []
        self.dupDihTypeKey = []
        self.molNum = 0
        self.k = 1
        self.stepRelax = False

    def __copy__(self):
        cls = self.__class__
        result = cls.__new__(cls)
        result.__dict__.update(self.__dict__)

    def __deepcopy__(self, memo):
        cls = self.__class__
        result = cls.__new__(cls)
        memo[id(self)] = result
        result.__dict__.update(deepcopy(self.__dict__))
        return result

    def setInfo(self, info):
        self.atomtypes = info[1]
        self.moleculetype = info[2]
        self.bondtypes = info[3]
        self.angletypes = info[4]
        self.dihtypes = info[5]
        self.imptypes = info[6]
        self.atoms = info[7]
        self.bonds = info[8]
        self.pairs = info[9]
        self.angles = info[10]
        self.dihedrals = info[11]
        self.impropers = info[12]
        self.system = info[0][0]
        self.molecules = info[0][1]
        self.dupDihTypeKey = info[13]
    
    def subAtom2Atypes(self, a1, a2, a3, a4, df_atoms, length=2):
        # a1 = row.ai; a2 = row.aj; a3 = row.ak; a4 = row.al
        a1Type = df_atoms.loc[int(a1)-1, 'type']
        a2Type = df_atoms.loc[int(a2)-1, 'type']   
        a3Type = df_atoms.loc[int(a3)-1, 'type']
        a4Type = df_atoms.loc[int(a4)-1, 'type']
        key = '{}-{}-{}-{}'.format(a1Type, a2Type, a3Type, a4Type)
        return key
    
    def mergeRow(self, x, keys='aTypes'):
        if keys == 'aTypes':
            str1 = '{:<4}{:>12}{:>10}{:>10}{:>3}{:>14}{:>14}'.format(
                x['name'], x.bond_type, round(float(x.mass), 2), round(float(x.charge), 2),
                x.ptype, x.sigma, x.epsilon)
        elif keys == 'mtypes':
            str1 = '{:>7}{:>7}'.format(x.loc['name'], x.loc['nrexcl'])
        elif keys == 'bTypes':
            str1 = '{:>7}{:>7}{:>7}{:>11}{:>11}'.format(
                x.ai, x.aj, x.funct, round(float(x.c0), 2), float(x.c1))
        elif keys == 'angTypes':
            str1 = '{:>7}{:>7}{:>7}{:>7} {:>10}{:>11}'.format(
                x.ai, x.aj, x.ak, x.funct, x.c0, x.c1)
        elif keys == 'dihTypes':
            if len(x) != 8:
                str1 = ' '
            else:
                str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}{:>7}'.format(
                    x.ai, x.aj, x.ak, x.al, x.funct, round(float(x.c0), 2), x.c1, x.c2)
        elif keys == 'atoms':
            str1 = '{:>5}{:>11}{:>7}{:>7}{:>7}{:>6} {:>11}{:>11}'.format(
                    x.nr, x.type, x.resnr, x.residue, x.atom, x.cgnr, x.charge, x.mass)
        elif keys == 'bonds':
            if self.stepRelax:
                if x.new == '':
                    str1 = '{:>7}{:>7}{:>7}'.format(
                        x.ai, x.aj, x.funct)
                else:
                    dist = float(x.c0) + (float(x.dist) - float(x.c0)) * (1 - self.k)
                    str1 = '{:>7}{:>7}{:>7}{:>11}{:>11}'.format(
                            x.ai, x.aj, x.funct, round(dist, 2), float(x.c1) * self.k)
            else:
                str1 = '{:>7}{:>7}{:>7}'.format(
                    x.ai, x.aj, x.funct)

        elif keys == 'pairs':
            if self.stepRelax:
                if x.new == '':
                    str1 = '{:>7}{:>7}{:>7}'.format(
                        x.ai, x.aj, x.funct)
                else:
                    str1 = ' '
            else:
                str1 = '{:>7}{:>7}{:>7}'.format(
                        x.ai, x.aj, x.funct)
        elif keys == 'angles':
            if self.stepRelax:
                if x.new == '':
                    str1 = '{:>7}{:>7}{:>7}{:>7}'.format(
                        x.ai, x.aj, x.ak, x.funct)
                else:
                    str1 = ' '
            else:
                str1 = '{:>7}{:>7}{:>7}{:>7}'.format(
                    x.ai, x.aj, x.ak, x.funct)

        elif keys == 'dih':
            # if len(x) != 8:
            #     str1 = ' '
            # else:
            key = self.subAtom2Atypes(x.ai, x.aj, x.ak, x.al, self.atoms)
            # key = '{}-{}-{}-{}'.format(x.ai, x.aj, x.ak, x.al)
            if self.stepRelax:
                if x.new == '':
                    if key not in self.dupDihTypeKey:
                        str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}'.format(
                            x.ai, x.aj, x.ak, x.al, x.funct)
                    else:
                        str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}{:>7}'.format(
                            x.ai, x.aj, x.ak, x.al, x.funct, round(float(x.c0), 2), x.c1, x.c2)
                else:
                    str1 = ' '
            else:
                if key not in self.dupDihTypeKey:
                    str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}'.format(
                            x.ai, x.aj, x.ak, x.al, x.funct)
                else:
                    str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}{:>7}'.format(
                            x.ai, x.aj, x.ak, x.al, x.funct, round(float(x.c0), 2), x.c1, x.c2)
        elif keys == 'imp':
            str1 = ' '

        return str1

    
    def avgCharge(self, incharge, add=True):
        c = incharge / len(self.atoms.charge)
        with open('chargeInfo.txt', 'a') as f:
            f.write('residue charge: {}\n'.format(c))

        for i in range(len(self.atoms.charge)):
            self.atoms.charge[i] = decimal.Decimal(self.atoms.charge[i]) - c

    def addCharge(self, incharge):
        c = max(self.atoms.charge)
        for i in range(len(self.atoms.charge)):
            if self.atoms.charge[i] == c:
                self.atoms.charge[i] = decimal.Decimal(self.atoms.charge[i]) - decimal.Decimal(incharge)
                break

    def setChargeDicimal(self, row):
        row.charge = decimal.Decimal(row.charge)
        # row.charge = round(float(row.charge), 8)
        return row

    def countCharge(self):
        charges = decimal.Decimal('0')
        for index, row in self.atoms.iterrows():
            charges += decimal.Decimal(row.charge)

        return charges

    def checkCharge(self):

        self.atoms = self.atoms.apply(lambda x: self.setChargeDicimal(x), axis=1)
        charges = self.countCharge()

        with open('chargeInfo.txt', 'a') as f:
            f.write('ori charge: {}\n'.format(charges))
            f.write('\t{}'.format(self.atoms.charge[:10]))
        self.avgCharge(charges)

        charges = self.countCharge()
        with open('chargeInfo.txt', 'a') as f:
            f.write('new charge: {}\n'.format(charges))
            f.write('\t{}'.format(self.atoms.charge[:10]))

        self.addCharge(charges)
        charges = self.countCharge()
        with open('chargeInfo.txt', 'a') as f:
            f.write('final charge: {}\n'.format(charges))
            f.write('\t{}'.format(self.atoms.charge[:10]))

    def addTopRow(self, df, inStr):
        df.loc[-1] = inStr
        df.index = df.index + 1
        df = df.sort_index().reset_index(drop=True)
        return df
    
    def addBondTypes(self, lst):
        cNames = ['ai', 'aj', 'func', 'c0', 'c1']
        tmp = pd.DataFrame([lst], columns=cNames)
        self.bondtypes = self.bondtypes.append(tmp, sort=False).reset_index(drop=True)

    
    def addAngleTypes(self, lst):
        cNames = ['ai', 'aj', 'ak', 'func', 'c0', 'c1']
        tmp = pd.DataFrame([lst], columns=cNames)
        self.angletypes = self.angletypes.append(tmp, sort=False).reset_index(drop=True)

    def addDihTypes(self, lst):
        cNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        tmp = pd.DataFrame([lst], columns=cNames)
        self.dihtypes = self.dihtypes.append(tmp, sort=False).reset_index(drop=True)

    @countTime
    def addBonds(self, pairs):
        # cNames = ['ai', 'aj', 'funct', 'new']
        cNames = ['ai', 'aj', 'funct', 'c0', 'c1', 'new', 'dist']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        # b_tmp['funct'] = '1'
        self.bonds['c0'] = ''
        self.bonds['c1'] = ''
        self.bonds['new'] = ''
        self.bonds['dist'] = ''
        self.bonds = pd.concat([self.bonds, b_tmp])

        # for p in pairs:
        #     a = p.copy()
        #     a.append('1')
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.bonds = self.bonds.append(b_tmp, sort=False).reset_index(drop=True)

    @countTime
    def addPairs(self, pairs):
        # cNames = ['ai', 'aj', 'funct']
        cNames = ['ai', 'aj', 'funct', 'new']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        self.pairs['new'] = ''
        self.pairs = pd.concat([self.pairs, b_tmp])
        # for p in pairs:
        #     a = p.copy()
        #     a.append('1')
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.pairs = self.pairs.append(b_tmp, sort=False).reset_index(drop=True)

    @countTime
    def addAngles(self, pairs):
        cNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1', 'new']
        # cNames = ['ai', 'aj', 'ak']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        self.angles['new'] = ''
        self.angles = pd.concat([self.angles, b_tmp])

    @countTime
    def addDih(self, pairs):
        cNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2', 'new']
        # cNames = ['ai', 'aj', 'ak', 'al']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        self.dihedrals['new'] = ''
        self.dihedrals = pd.concat([self.dihedrals, b_tmp])

        # for p in pairs:
        #     a = p.copy()
        #     a.append('9')
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.dihedrals = self.dihedrals.append(b_tmp, sort=False).reset_index(drop=True)
    
    def addImp(self, pairs):
        # cNames = ['ai', 'aj', 'ak', 'al', 'funct']
        cNames = ['ai', 'aj', 'ak', 'al']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        b_tmp['funct'] = '4'
        self.dihedrals = pd.concat([self.dihedrals, b_tmp])

        # for p in pairs:
        #     a = p.copy()
        #     a.append('4')
        #     b_tmp = pd.DataFrame([pairs], columns=cNames)
        #     self.dihedrals = self.dihedrals.append(b_tmp, sort=False).reset_index(drop=True)
    
    def outTop(self, df, outName):
        with open('{}.top'.format(outName), 'w') as f:
            for index, row in df.iterrows():
                f.write('{}\n'.format(row['0']))
        
    def outDf(self, outName, k=1, simple=False, stepRelax=False):
        self.k = k
        self.stepRelax = stepRelax
        self.bonds.reset_index(drop=True)
        self.angles.reset_index(drop=True)
        self.pairs.reset_index(drop=True)
        self.dihedrals.reset_index(drop=True)

        df_atypes_str = self.atomtypes.apply(lambda x: self.mergeRow(x, keys='aTypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_btypes_str = self.bondtypes.apply(lambda x: self.mergeRow(x, keys='bTypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_angTypes_str = self.angletypes.apply(lambda x: self.mergeRow(x, keys='angTypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_dihTypes_str = self.dihtypes.apply(lambda x: self.mergeRow(x, keys='dihTypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_impTypes_str = self.imptypes.apply(lambda x: self.mergeRow(x, keys='dihTypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_atoms_str = self.atoms.apply(lambda x: self.mergeRow(x, keys='atoms'), axis=1).to_frame().rename(columns={0: '0'})
        df_bonds_str = self.bonds.apply(lambda x: self.mergeRow(x, keys='bonds'), axis=1).to_frame().rename(columns={0: '0'})
        df_pairs_str = self.pairs.apply(lambda x: self.mergeRow(x, keys='pairs'), axis=1).to_frame().rename(columns={0: '0'})
        df_angles_str = self.angles.apply(lambda x: self.mergeRow(x, keys='angles'), axis=1).to_frame().rename(columns={0: '0'})
        df_dih_str = self.dihedrals.apply(lambda x: self.mergeRow(x, keys='dih'), axis=1).to_frame().rename(columns={0: '0'})
        df_imp_str = self.impropers.apply(lambda x: self.mergeRow(x, keys='imp'), axis=1).to_frame().rename(columns={0: '0'})
        df_default = self.default
        df_default = df_default.rename(columns={0: '0'})
        df_molType = self.moleculetype.apply(lambda x: self.mergeRow(x, keys='mtypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_sys = self.system
        df_mol = self.molecules
        df_itp = pd.DataFrame(['#include "{}.itp"'.format(outName)], columns=['0'])
        # Top file section: default, atomtype, bondtype, angletype, dihtype, molecular type, 
        # atom, bond, pair, angle, dihedral, system, molecules
        df_lst0 = [df_default, df_itp, df_sys, df_mol]
        str_top_tmp = ['[ defaults ]', '; Include', '[ system ]', '[ molecules ]']
        if simple:
            df_lst1 = [df_atypes_str, df_btypes_str, df_molType, df_atoms_str, df_bonds_str]
            str_itp_tmp = ['[ atomtypes ]', '[ bondtypes ]', '[ moleculetype ]', '[ atoms ]', '[ bonds ]']

        else:
            df_lst1 = [df_atypes_str, df_btypes_str, df_angTypes_str, df_dihTypes_str, df_impTypes_str, df_molType,
                       df_atoms_str, df_bonds_str, df_pairs_str, df_angles_str, df_dih_str, df_imp_str]
            str_itp_tmp = ['[ atomtypes ]', '[ bondtypes ]', '[ angletypes ]', '[ dihedraltypes ]', '[ dihedraltypes ]',
                           '[ moleculetype ]', '[ atoms ]', '[ bonds ]', '[ pairs ]', '[ angles ]', '[ dihedrals ]',
                           '[ dihedrals ]',
                           ]
        df0 = []; df = []

        for i in range(len(str_top_tmp)):
            s = str_top_tmp[i]
            df_tmp = self.addTopRow(df_lst0[i], s)
            df0.append(df_tmp)
            
        for i in range(len(str_itp_tmp)):
            s = str_itp_tmp[i]
            df_tmp = self.addTopRow(df_lst1[i], s)
            df.append(df_tmp)
        
        df_out_top = pd.concat(df0)
        df_out_itp = pd.concat(df) #TODO: need to drop empty line

        self.outTop(df_out_top, outName) # Not use to_csv, since it cannot handle " well, it will always generate two quotation marks
        df_out_itp.to_csv('{}.itp'.format(outName), mode = 'w', index=False, header=None)
        return df_out_top
