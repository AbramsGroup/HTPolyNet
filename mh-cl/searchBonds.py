# -*- coding: utf-8 -*-
"""
Searching potential bonds

@author: huang
"""
import pandas as pd
import numpy as np
import time
from countTime import *
from multiprocessing import Process, Pool
from functools import partial

import random
import sys

class searchBonds(object):
    def __init__(self, basicParameters, generatedBonds, inGro, inTop):
        self.monInfo = basicParameters.monInfo
        self.croInfo = basicParameters.croInfo
        self.cutoff = basicParameters.cutoff
        self.bondsRatio = basicParameters.bondsRatio
        self.boxSize = basicParameters.boxSize  # one dimension box length
        self.maxBonds = basicParameters.maxBonds
        self.rctInfo = basicParameters.rctInfo
        self.rctMap = self.genRctMap(basicParameters.rctInfo)
        self.generatedBonds = generatedBonds
        self.gro = inGro  # gro object
        self.top = inTop  # top object
        self.rctMon = []
        self.rctCro = []
        self.mol = []

        # links-cell properties
        self.maxCellId = []
        self.cellId = []
        self.div_box = []

    def genRctMap(self, rctInfo):
        rctMap = {}
        for rct in rctInfo:
            if rct[0][0] not in rctMap.keys():
                rctMap[rct[0][0]] = [[rct[1][0], rct[2][0]]]
            elif rct[0][0] in rctMap.keys():
                if [rct[1][0], rct[2][0]] not in rctMap[rct[0][0]]:
                    rctMap[rct[0][0]].append([rct[1][0], rct[2][0]])

            if rct[1][0] not in rctMap.keys():
                rctMap[rct[1][0]] = [[rct[0][0], rct[2][0]]]
            elif rct[1][0] in rctMap.keys():
                if [rct[0][0], rct[2][0]] not in rctMap[rct[1][0]]:
                    rctMap[rct[1][0]].append([rct[0][0], rct[2][0]])
            else:
                print('undecide rctions. Please check!')

        return rctMap

    def calDist(self, aPos, bPos, pbc=True):
        boxSize = [self.boxSize, self.boxSize, self.boxSize]
        xlen = 0.5 * float(boxSize[0])
        ylen = 0.5 * float(boxSize[1])
        zlen = 0.5 * float(boxSize[2])

        x1 = float(aPos[0]) - float(bPos[0])
        y1 = float(aPos[1]) - float(bPos[1])
        z1 = float(aPos[2]) - float(bPos[2])
        if pbc:
            if x1 > xlen:
                x1 = x1 - float(boxSize[0])
            elif x1 <= -xlen:
                x1 = x1 + float(boxSize[0])

            elif y1 > ylen:
                y1 = y1 - float(boxSize[1])
            elif y1 <= -ylen:
                y1 = y1 + float(boxSize[1])

            elif z1 > zlen:
                z1 = z1 - float(boxSize[2])
            elif z1 <= -zlen:
                z1 = z1 + float(boxSize[2])

        x2 = np.power(x1, 2)
        y2 = np.power(y1, 2)
        z2 = np.power(z1, 2)

        dist = np.sqrt(x2 + y2 + z2)
        return dist

    def filterRctPairs(self, x, df_mon):
        lst = []
        for index, row in df_mon.iterrows():
            if x.rctNum > 0 and row.rctNum > 0:
                dist = self.calDist([x.posX, x.posY, x.posZ], [row.posX, row.posY, row.posZ])
                if dist < self.cutoff and dist > 0.05:
                    lst.append([x.globalIdx, row.globalIdx, round(dist, 2), x.rctNum, row.rctNum, x.molCon])
        return lst

    def getId(self, x, maxId):
        outList = []
        tmpLst = [-1, 0, 1]
        for i in tmpLst:
            tmp = int(x) + i
            # PBC condition applied
            if tmp < 0:
                tmp = maxId
            elif tmp > maxId:
                tmp = 0
            else:
                pass
            outList.append(tmp)
        return outList

    def appendDist(self, x, atoms, boxSize):
        x0 = float(atoms.posX)
        y0 = float(atoms.posY)
        z0 = float(atoms.posZ)
        x1 = float(x.posX)
        y1 = float(x.posY)
        z1 = float(x.posZ)
        x['dist'] = self.calDist([x0, y0, z0], [x1, y1, z1])
        return x

    def condFilter(self, atom, atomsDf):
        # Reaction define
        rctMolInfo = self.rctMap[atom.molName]
        rctName = []; rctPct = []
        for r in rctMolInfo:
            if r[0] not in rctName:
                rctName.append(r[0])
                rctPct.append(r[1])

        tmpDf = atomsDf.loc[atomsDf.molName.isin(rctName)]
        for i in range(len(rctName)):
            tmpDf.loc[(tmpDf.molName == rctName[i]), 'rctPct'] = rctPct[i]

        return tmpDf


    def getPairs(self, atom, atomsDf):  # collect atoms based on cell id. itself and adjacent cell
        # atomsDf contains all atoms
        # maxCellId used for pbc condition
        cell0 = atom.cellId
        # print('cell0: ', cell0)
        maxCellId = self.maxCellId
        tmpLst = [-1, 0, 1]
        xList = self.getId(cell0[0], maxCellId[0])
        yList = self.getId(cell0[1], maxCellId[1])
        zList = self.getId(cell0[2], maxCellId[2])

        cellSum = []
        for i in xList:
            for ii in yList:
                for iii in zList:
                    id = ''.join([str(i), str(ii), str(iii)])
                    cellSum.append(id)

        df1 = atomsDf.loc[atomsDf.cellId.isin(cellSum)]
        df2 = df1.loc[df1.rctNum > 0]

        # Using needed condition to filter the atomsDf
        df2 = self.condFilter(atom, df2)

        # obtain the distance between the atoms and distance between atoms and the potential atoms within neighbor cell
        df3 = df2.apply(lambda x: self.appendDist(x, atom, self.boxSize), axis=1)
        df3 = df3.iloc[1:] # remove atoms which connect to itself
        pd.to_numeric(df3.dist)
        df4 = df3.loc[(df3.dist < float(self.cutoff)) & (df3.molNum != atom.molNum)]
        df_out1 = []

        if len(df4) == 0:
            return []
        else:
            for index, row in df4.iterrows():
                df_out1.append([atom.globalIdx, row.globalIdx, round(row.dist, 2), atom.rctNum, row.rctNum,
                                atom.molCon, random.random(), row.rctPct])
            return df_out1

    @countTime
    def getAllMolNames(self):
        atomNames = []
        rctFunc = self.rctInfo
        for i in rctFunc:
            if i[0] not in atomNames:
                atomNames.append(i[0][0])
            elif i[1] not in atomNames:
                atomNames.append(i[1][0])

        return atomNames

    def openList(self, inList):
        out = []
        for i in inList:
            if i == []:
                continue
            else:
                for ii in i:
                    out.append(ii)
        return out

    def parallel_getPairs(self, df, df_sum):
        res = df.apply(lambda x: self.getPairs(x, df_sum), axis=1)
        return res

    @countTime
    def checkHydrogen(self, atomDf):
        bonds = self.top.bonds
        atoms = self.top.atoms
        filterIdx = []
        for index, value in atomDf.iterrows():
            atIdx = str(value.globalIdx)
            atCon = bonds.loc[(bonds.ai == atIdx) | (bonds.aj == atIdx)]
            for idx2, jj in atCon.iterrows():
                a1 = atoms.loc[atoms.nr == jj.ai]
                a2 = atoms.loc[atoms.nr == jj.aj]
                if 'H' in a1.atom.to_list()[0] or 'H' in a2.atom.to_list()[0]:
                    filterIdx.append(atIdx)
        outDf = atomDf.loc[atomDf.globalIdx.isin(filterIdx)]
        return outDf


    @countTime
    def getRctDf(self):
        rctFunc = self.rctInfo
        atoms = self.gro.df_atoms
        top = self.top
        rctDf = []
        if 'molCon' in atoms.keys():
            pass
        else:
            atoms['molCon'] = '0'

        top.atoms['molCon'] = atoms['molCon']
        top.atoms['molNum'] = atoms['molNum']

        atomNames = self.getAllMolNames()
        df_tmp0 = atoms.loc[(atoms.rct == 'True') & (atoms.molName.isin(atomNames))]
        df_tmp0 = self.assignAtoms(df_tmp0)

        df_tmp = self.checkHydrogen(df_tmp0) # This check just to confirm selected atom can react
        # ##### START PARALLEL
        print('start parallel searching!!')
        p = Pool(processes=4) #TODO: should be able to tune based on the number of cell and free CPU cores
        dfSplit = np.array_split(df_tmp, 4)
        results = p.map(partial(self.parallel_getPairs, df_sum=df_tmp), dfSplit)
        p.close()
        p.join()
        parts = pd.concat(results, axis=0)
        # ##### END PARALLEL
        # Non-parallel method to search bonds
        # lst_tmp = df_tmp.apply(lambda x: self.getPairs(x, df_tmp), axis=1)

        lst_tmp = parts
        lst_tmp2 = self.openList(lst_tmp)
        for ii in lst_tmp2:
            if len(ii) == 0:
                continue
            else:
                rctDf.append(ii)

        names = ['acro', 'amon', 'dist', 'rctNum1', 'rctNum2', 'molCon', 'p', 'rctP']
        df = pd.DataFrame(rctDf, columns=names)
        return df

    def checkRepeat(self, lst, inLst):
        for i in lst:
            if set(i) == set(inLst):
                return False
        return True

    def checkConSingleMolecule(self, row, atom):
        if row.molNum == atom.molNum:
            return False
        else:
            return True

    def checkCircuit(self, df_rctAtoms, row):
        mol1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molNum'].to_list()
        mol2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molNum'].to_list()
        molCon1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molCon'].to_list()[0].split(',')
        molCon2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molCon'].to_list()[0].split(',')

        if mol1[0] in molCon2 or mol2[0] in molCon1 or mol1[0] == mol2[0]:
            return False
        else:
            return True

    def checkAtomsRepeat(self, lst, inAtoms):
        for a in inAtoms:
            if a in lst:
                return False
            else:
                continue
        return True

    def updateMolCon(self, atomsDf, a1, a2):
        mol1 = atomsDf[atomsDf.globalIdx == str(a1)].molNum.values[0]
        mol2 = atomsDf[atomsDf.globalIdx == str(a2)].molNum.values[0]

        con1 = list(atomsDf[atomsDf.molNum == mol1].molCon)[0]
        con2 = list(atomsDf[atomsDf.molNum == mol2].molCon)[0]

        tmp = '{},{}'.format(con1, mol2)
        atomsDf.loc[atomsDf['molNum'] == mol1, 'molCon'] = tmp

        tmp = '{},{}'.format(con2, mol1)
        atomsDf.loc[atomsDf['molNum'] == mol2, 'molCon'] = tmp

    @countTime
    def finalRctPairs(self, df_pairs):
        '''
        1. check repeat
        2. check circuit
        3. check kinetic ratio
        '''
        atomsDf = self.gro.df_atoms
        # Remove repeat rows
        lst = [];
        rowList = []
        pcriteria = 0
        if len(df_pairs) <= 1:
            pcriteria = 1

        for index, row in df_pairs.iterrows():
            if self.checkRepeat(lst, [row.acro, row.amon]):
                if pcriteria == 1:
                    k = 1
                else:
                    k = float(row.rctP.strip('%'))/100

                if float(row.p) > 1 - k:
                    # lst.append([row.acro, row.amon, random.random(), row.p]);
                    rowList.append(row)
            else:
                continue
        df_tmp1 = pd.DataFrame(rowList)
        # check circuit connection. Molecules cannot connect to the same molecules
        rowList = [];
        atomsList = []
        for index, row in df_tmp1.iterrows():
            if self.checkCircuit(atomsDf, row):
                if self.checkAtomsRepeat(atomsList, [row.acro, row.amon]):
                    atomsList.append(row.acro);
                    atomsList.append(row.amon);
                    rowList.append(row)
                    self.updateMolCon(atomsDf, row.acro, row.amon)
            else:
                continue
        df_tmp2 = pd.DataFrame(rowList)

        return df_tmp2

    def idx2Mol(self, pairs):
        # pairs is a dataframe
        atoms = [];
        mols = []
        for index, row in pairs.iterrows():
            atoms.append(row.acro)
            atoms.append(row.amon)

        groDf = self.gro.df_atoms
        for a in atoms:
            index = list(groDf[groDf.globalIdx == str(a)].molNum)[0]
            mols.append(index)

        self.mol = mols

    # links-cell algorithm
    def genCell(self, parts):
        parts = parts - 1
        xdiv = 0
        if self.cutoff > float(self.boxSize) * 0.5:
            print('cutoff should smaller than the half of the box size, please modify it.')
            sys.exit()

        while (xdiv <= self.cutoff):
            print('parts: ', parts)
            if parts <= 1:
                print('Unusual occasion occurs! Please check the cutoff you set and rerun the script!')
                sys.exit()

            x = np.linspace(0, float(self.boxSize), parts + 1)
            y = np.linspace(0, float(self.boxSize), parts + 1)
            z = np.linspace(0, float(self.boxSize), parts + 1)
            xdiv = x[1] - x[0]
            ydiv = y[1] - y[0]
            zdiv = z[1] - z[0]
            parts -= 1

        cell_id = [];
        div_box = [xdiv, ydiv, zdiv]
        com = [0.5 * (x[1] - x[0]),
               0.5 * (y[1] - y[0]),
               0.5 * (z[1] - z[0])]
        xNum = 0;
        yNum = 0;
        zNum = 0
        xMax = 0;
        yMax = 0;
        zMax = 0

        for i in range(parts):  # x dir
            xCom = com[0] + xdiv * xNum
            yMax = yNum;
            zMax = zNum
            xNum += 1;
            yNum = 0;
            zNum = 0
            for j in range(parts):  # y dir
                yCom = com[1] + ydiv * yNum
                yNum += 1;
                zNum = 0
                for k in range(parts):  # z dir
                    zCom = com[2] + zdiv * zNum
                    zNum += 1
                    id = ''.join([str(i), str(j), str(k)])
                    info = [id, [xCom, yCom, zCom]]  # [[id], [center coord]]
                    cell_id.append(info)
            xMax = xNum
            maxCellId = [xMax - 1, yMax - 1, zMax - 1]  # Since it starts from 0
        self.cellId = cell_id
        self.div_box = div_box
        self.maxCellId = maxCellId

    def searchCell(self, row):
        cell_id = self.cellId
        box_div = self.div_box
        xDiv = 0.5 * box_div[0]
        yDiv = 0.5 * box_div[1]
        zDiv = 0.5 * box_div[2]
        row['cellId'] = '000'
        row['comCoord'] = '0'
        coord = [float(row.posX), float(row.posY), float(row.posZ)]
        for c in cell_id:
            xlow, xhigh = [c[1][0] - xDiv, c[1][0] + xDiv]
            ylow, yhigh = [c[1][1] - yDiv, c[1][1] + yDiv]
            zlow, zhigh = [c[1][2] - zDiv, c[1][2] + zDiv]
            if xlow <= coord[0] <= xhigh and ylow <= coord[1] <= yhigh and zlow <= coord[2] <= zhigh:
                row['cellId'] = c[0]
                row['comCoord'] = c[1]
                return row
            else:
                continue
        return row

    @countTime
    def assignAtoms(self, df_atoms):
        df1 = df_atoms.apply(lambda x: self.searchCell(x), axis=1)
        return df1

    @countTime
    def main(self):
        pairs = [];
        count = 2
        parts = 8
        # print('Generate cells number on one dimension: ', parts)
        self.genCell(parts)
        while (len(pairs) == 0):
            for i in range(count):  # max trial times
                df_pairs = self.getRctDf()
                if len(df_pairs) > 0:
                    break
            a1 = self.finalRctPairs(df_pairs)
            if len(a1) == 0:
                self.cutoff += 0.5
                if self.cutoff > 0.5 * float(self.boxSize):
                    break
            pairs = a1

        if len(pairs) == 0:
            return [], self.mol

        self.idx2Mol(pairs)
        print('{} bonds are going to be generated'.format(len(pairs)))
        print('Following bonds will be formed: \n')
        for index, value in pairs.iterrows():
            print('\t', value.acro, '\t', value.amon, '\t',
                  round(value.p, 2), '\t', value.rctP)
        return a1, self.mol
