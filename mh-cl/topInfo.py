# -*- coding: utf-8 -*-
"""
Created on Mon Sep 28 15:34:47 2020

@author: huangming
"""
import pandas as pd
from decimal import Decimal
from countTime import *

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
                x.ai, x.aj, x.funct, x.c0, round(float(x.c1), 2))
        elif keys == 'angTypes':
            str1 = '{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}'.format(
                x.ai, x.aj, x.ak, x.funct, x.c0, round(float(x.c1), 2))
        elif keys == 'dihTypes':
            str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}{:>7}'.format(
                x.ai, x.aj, x.ak, x.al, x.funct, round(float(x.c0), 2), round(float(x.c1), 2), x.c2)
        elif keys == 'atoms':
            str1 = '{:>5}{:>11}{:>7}{:>7}{:>7}{:>6}{:>11}{:>11}'.format(
                    x.nr, x.type, x.resnr, x.residue, x.atom, x.cgnr, x.charge, x.mass)
        elif keys == 'bonds':
            str1 = '{:>7}{:>7}{:>7}'.format(
                    x.ai, x.aj, x.funct)
        elif keys == 'pairs':
            str1 = '{:>7}{:>7}{:>7}'.format(
                    x.ai, x.aj, x.funct)
        elif keys == 'angles':
            str1 = '{:>7}{:>7}{:>7}{:>7}'.format(
                    x.ai, x.aj, x.ak, x.funct)
        elif keys == 'dih':
            key = self.subAtom2Atypes(x.ai, x.aj, x.ak, x.al, self.atoms)
            # key = '{}-{}-{}-{}'.format(x.ai, x.aj, x.ak, x.al)
            if key not in self.dupDihTypeKey:
                str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}'.format(
                        x.ai, x.aj, x.ak, x.al, x.funct)
            else:
                str1 = '{:>7}{:>7}{:>7}{:>7}{:>7}{:>11}{:>11}{:>7}'.format(
                        x.ai, x.aj, x.ak, x.al, x.funct, round(float(x.c0), 2), round(float(x.c1), 4), x.c2)
        return str1
    
    def addCharge(self, incharge):
        c = max(self.atoms.charge)
        for i in range(len(self.atoms.charge)):
            if self.atoms.charge[i] == c:
                # print('old: {}, new: {}'.format(self.atoms.charge[i], str(Decimal(self.atoms.charge[i]) + Decimal(incharge))))
                self.atoms.charge[i] = str(Decimal(self.atoms.charge[i]) - Decimal(incharge))
                break
        
    def setChargeDicimal(self, row):
        row.charge = str(round(float(row.charge), 4))
        return row
        
    def checkCharge(self):
        self.atoms = self.atoms.apply(lambda x: self.setChargeDicimal(x), axis=1)
        charges = 0
        for index, row in self.atoms.iterrows():
            charges += Decimal(row.charge)
        
        self.addCharge(charges)
        
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
        # cNames = ['ai', 'aj', 'funct']
        cNames = ['ai', 'aj']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        b_tmp['funct'] = '1'
        self.bonds = pd.concat([self.bonds, b_tmp])
        # for p in pairs:
        #     a = p.copy()
        #     a.append('1')
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.bonds = self.bonds.append(b_tmp, sort=False).reset_index(drop=True)

    @countTime
    def addPairs(self, pairs):
        # cNames = ['ai', 'aj', 'funct']
        cNames = ['ai', 'aj']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        b_tmp['funct'] = '1'
        self.pairs = pd.concat([self.pairs, b_tmp])
        # for p in pairs:
        #     a = p.copy()
        #     a.append('1')
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.pairs = self.pairs.append(b_tmp, sort=False).reset_index(drop=True)

    @countTime
    def addAngles(self, pairs):
        # cNames = ['ai', 'aj', 'ak', 'funct']
        cNames = ['ai', 'aj', 'ak']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        b_tmp['funct'] = '1'
        self.angles = pd.concat([self.angles, b_tmp])

        # for p in pairs:
        #     a = p.copy()
        #     a.append('1') # This is the func
        #     b_tmp = pd.DataFrame([a], columns=cNames)
        #     self.angles = self.angles.append(b_tmp, sort=False).reset_index(drop=True)

    @countTime
    def addDih(self, pairs):
        # cNames = ['ai', 'aj', 'ak', 'al', 'funct']
        cNames = ['ai', 'aj', 'ak', 'al']
        b_tmp = pd.DataFrame(pairs, columns=cNames)
        b_tmp['funct'] = '9'
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
        
    def outDf(self, outName):
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
        df_imp_str = self.impropers.apply(lambda x: self.mergeRow(x, keys='dih'), axis=1).to_frame().rename(columns={0: '0'})
        df_default = self.default
        df_default = df_default.rename(columns={0: '0'})
        df_molType = self.moleculetype.apply(lambda x: self.mergeRow(x, keys='mtypes'), axis=1).to_frame().rename(columns={0: '0'})
        df_sys = self.system
        df_mol = self.molecules
        df_itp = pd.DataFrame(['#include "{}.itp"'.format(outName)], columns=['0'])
        # Top file section: default, atomtype, bondtype, angletype, dihtype, molecular type, 
        # atom, bond, pair, angle, dihedral, system, molecules
        df_lst0 = [df_default, df_itp, df_sys, df_mol]
        df_lst1 = [df_atypes_str, df_btypes_str, df_angTypes_str, df_dihTypes_str, df_impTypes_str, df_molType,
                   df_atoms_str, df_bonds_str, df_pairs_str, df_angles_str, df_dih_str, df_imp_str]
        
        df0 = []; df = []
        str_top_tmp = ['[ defaults ]', '; Include', '[ system ]', '[ molecules ]']
        str_itp_tmp = ['[ atomtypes ]', '[ bondtypes ]', '[ angletypes ]', '[ dihedraltypes ]', '[ dihedraltypes ]', 
                   '[ moleculetype ]', '[ atoms ]', '[ bonds ]', '[ pairs ]', '[ angles ]', '[ dihedrals ]', '[ dihedrals ]', 
                   ]
        for i in range(len(str_top_tmp)):
            s = str_top_tmp[i]
            df_tmp = self.addTopRow(df_lst0[i], s)
            df0.append(df_tmp)
            
        for i in range(len(str_itp_tmp)):
            s = str_itp_tmp[i]
            df_tmp = self.addTopRow(df_lst1[i], s)
            df.append(df_tmp)
        
        df_out_top = pd.concat(df0)
        df_out_itp = pd.concat(df)
        self.outTop(df_out_top, outName) # Not use to_csv, since it cannot handle " well, it will always generate two quotation marks
        df_out_itp.to_csv('{}.itp'.format(outName), mode = 'w', index=False, header=None)
        return df_out_top