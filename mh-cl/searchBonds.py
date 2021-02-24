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
    def __init__(self, basicParameters, generatedBonds, inGro, inTop, conv, desBonds):
        self.monInfo = basicParameters.monInfo
        self.croInfo = basicParameters.croInfo
        self.cutoff = basicParameters.cutoff
        self.bondsRatio = basicParameters.bondsRatio
        self.boxSize = basicParameters.boxSize  # one dimension box length
        self.maxBonds = basicParameters.maxBonds
        self.rctInfo = basicParameters.rctInfo
        self.HTProcess = basicParameters.HTProcess
        self.rctMap = self.genRctMap(basicParameters.rctInfo)
        self.generatedBonds = generatedBonds
        self.gro = inGro  # gro object
        self.top = inTop  # top object
        self.rctMon = []
        self.rctCro = []
        self.mol = []
        self.conv = conv
        self.desBonds = desBonds
        self.dumpPairs = {}

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
        # dump pairs define, filter short circuit
        try:
            dumpList = self.dumpPairs[atom.globalIdx]
        except:
            dumpList = []
        df1 = atomsDf.loc[~atomsDf.molNum.isin(dumpList)]
        # Reaction define
        rctMolInfo = self.rctMap[atom.molName]
        rctName = []; rctPct = []

        for r in rctMolInfo:
            if r[0] not in rctName:
                rctName.append(r[0])
                rctPct.append(r[1])

        tmpDf = atomsDf.loc[df1.molName.isin(rctName)]
        for i in range(len(rctName)):
            tmpDf.loc[(tmpDf.molName == rctName[i]), 'rctPct'] = rctPct[i]

        return tmpDf

    def getPairs(self, atom, atomsDf):  # collect atoms based on cell id. itself and adjacent cell
        # atomsDf contains all atoms
        # maxCellId used for pbc condition
        cell0 = atom.cellId
        maxCellId = self.maxCellId
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
        df1.to_csv('df1.csv')
        df2 = df1.loc[df1.rctNum > 0]

        # Using needed condition to filter the atomsDf
        # Set chemical reaction
        df2 = self.condFilter(atom, df2)

        # obtain the distance between the atoms and distance between atoms and the potential atoms within neighbor cell
        df3 = df2.apply(lambda x: self.appendDist(x, atom, self.boxSize), axis=1)
        df3 = df3.iloc[1:] # remove atoms which connect to itself
        pd.to_numeric(df3.dist)
        df3.sort_values(by='dist', inplace=True)
        df3.to_csv('df3_sorted.csv')
        print('cutoff: ', self.cutoff)
        df4 = df3.loc[(df3.dist < float(self.cutoff)) & (df3.molNum != atom.molNum)]
        df4.to_csv('df4.csv')
        df_out1 = []

        if len(df4) == 0:
            print('No pairs found')
            return [[]]
        else:
            print('Pairs found')
            for index, row in df4.iterrows():
                df_out1.append([atom.globalIdx, row.globalIdx, round(row.dist, 2), atom.rctNum, row.rctNum,
                                atom.molCon, random.random(), row.rctPct])
            return df_out1

    @countTime
    def getAllMolNames(self):
        atomNames = []
        rctFunc = self.rctInfo
        for i in rctFunc:
            if i[0][0] not in atomNames:
                atomNames.append(i[0][0])
            elif i[1][0] not in atomNames:
                atomNames.append(i[1][0])
        return atomNames

    def openList(self, inList):
        out = []
        print('inList type: ', type(inList))
        print('inList: ', inList)
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
        print('df_tmp: ', df_tmp.head())
        results = p.map(partial(self.parallel_getPairs, df_sum=df_tmp), dfSplit)
        p.close()
        p.join()
        tmpLst = []
        for l in results:
            print('l type: ', type(l))
            print('l: ', l)
            if l.empty:
                tmpLst.append(pd.Series(dtype=object))
            else:
                tmpLst.append(l)

        parts = pd.concat(tmpLst, axis=0)
        # ##### END PARALLEL
        # Non-parallel method to search bonds
        # parts = df_tmp.apply(lambda x: self.getPairs(x, df_tmp), axis=1)

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
            print('generate circuit, mol {}'.format(mol1[0]))
            keys = self.dumpPairs.keys()
            if row.acro in keys:
                self.dumpPairs[row.acro].append(row.amon)
            else:
                self.dumpPairs[row.acro] = [row.amon]

            if row.amon in keys:
                self.dumpPairs[row.amon].append(row.acro)
            else:
                self.dumpPairs[row.acro] = [row.acro]
            return False
        else:
            return True

    def checkAtomsRepeat(self, inAtoms, lst):
        for a in lst:
            if a in inAtoms:
                return False
            else:
                continue
        return True

    def checkHT(self, at1Idx, at2Idx):
        if self.HTProcess == 'False':
            return True
        else:
            atomsDf = self.gro.df_atoms
            at1Prop = atomsDf.loc[atomsDf.globalIdx == at1Idx].prop.values[0]
            at2Prop = atomsDf.loc[atomsDf.globalIdx == at2Idx].prop.values[0]
            # print('at1Prop: ', at1Prop.values[0])
            # print('at2Prop: ', at2Prop.values[0])
            if at1Prop == at2Prop:
                return False
            else:
                return True

    def checkKineticRatio(self, pcriteria, row):
        if pcriteria == 1:
            k = 1
        else:
            k = float(row.rctP.strip('%')) / 100

        if float(row.p) > 1 - k:
            return True
        else:
            return False

    def updateMolCon(self, atomsDf, a1, a2):
        mol1 = atomsDf[atomsDf.globalIdx == str(a1)].molNum.values[0]
        mol2 = atomsDf[atomsDf.globalIdx == str(a2)].molNum.values[0]

        con1 = list(atomsDf[atomsDf.molNum == mol1].molCon)[0]
        con2 = list(atomsDf[atomsDf.molNum == mol2].molCon)[0]

        tmp = '{},{}'.format(con1, mol2)
        atomsDf.loc[atomsDf['molNum'] == mol1, 'molCon'] = tmp

        tmp = '{},{}'.format(con2, mol1)
        atomsDf.loc[atomsDf['molNum'] == mol2, 'molCon'] = tmp

    def setRctP(self, df_pairs):
        for index, value in df_pairs.iterrows():
            df_pairs.loc[index, 'p'] = random.random()

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
        rowList0 = []
        rowList1 = []
        for index, row in df_pairs.iterrows():
            if self.checkHT(row.acro, row.amon):
                if self.checkRepeat(lst, [row.acro, row.amon]):
                    rowList0.append(row)
            else:
                continue

        df_tmp0 = pd.DataFrame(rowList0)
        df_tmp0.to_csv('all_bonds_within_cutoff.csv')

        print('{} bonds goes to check probability'.format(len(df_tmp0)))
        pcriteria = 0
        if len(df_tmp0) <= 2:
            pcriteria = 1

        cc = 0
        while(len(rowList1) == 0 and cc < 5):
            for index, row in df_tmp0.iterrows():
                if self.checkKineticRatio(pcriteria, row):
                    rowList1.append(row)

            if self.conv < 0.7:
                if len(rowList1) < 0.1 * self.desBonds:
                    self.setRctP(df_tmp0)
                    rowList = []
                    cc += 1
                    continue
                else:
                    break
            elif 0.7 <= self.conv <= 0.9:
                if len(rowList1) < 0.05 * self.desBonds:
                    self.setRctP(df_tmp0)
                    rowList = []
                    cc += 1
                    continue
                else:
                    break

            else:
                break

        df_tmp1 = pd.DataFrame(rowList1)
        df_tmp2 = self.updateAllCon(df_pairs)

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
        if self.cutoff > float(self.boxSize) * 0.6:
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

    def updateAllCon(self, df_tmp):
        atomsDf = self.gro.df_atoms
        rowList = []
        atomsList = []
        for index, row in df_tmp.iterrows():
            if self.checkCircuit(atomsDf, row):
                if self.checkAtomsRepeat(atomsList, [row.acro, row.amon]):
                    atomsList.append(row.acro);
                    atomsList.append(row.amon);
                    rowList.append(row)
                    self.updateMolCon(atomsDf, row.acro, row.amon)
            else:
                continue

        df_out = pd.DataFrame(rowList)
        df_out.to_csv('tmp2.csv')
        return df_out

    def collectBonds(self, count):
        pairs = []
        a1 = []
        count = 0
        while (len(pairs) == 0 and count < 10):
            while (len(pairs) == 0):
                df_pairs = self.getRctDf()
                df_pairs.to_csv('all_bonds_within_cutoff.csv')
                # df_pairs.to_csv('final_bonds.csv')
                print('df_pairs: ', df_pairs)
                if self.conv < 0.7:
                    if len(df_pairs) == 0 or len(df_pairs) < 0.2 * self.desBonds:
                        self.cutoff += 0.1
                        if self.cutoff > 0.5 * float(self.boxSize):
                            break
                        else:
                            continue
                    else:
                        break

                elif 0.7 <= self.conv <= 0.9:
                    if len(df_pairs) == 0 or len(df_pairs) < 0.2 * self.desBonds:
                        self.cutoff += 0.1
                        if self.cutoff > 0.6 * float(self.boxSize):
                            break
                        else:
                            continue
                    else:
                        break
                else:
                    if len(df_pairs) == 0:
                        self.cutoff += 0.1
                        if self.cutoff > 0.6 * float(self.boxSize):
                            break
                        else:
                            continue
                    else:
                        break

            df_pairs = self.finalRctPairs(df_pairs)  # TODO: update mol connection should not be here
            pairs = df_pairs
            count += 1
            if len(pairs) == 0:
                self.cutoff += 0.2
        return pairs

    @countTime
    def main(self):
        count = 5
        parts = 8
        self.genCell(parts)
        pairs = self.collectBonds(count)
        if len(pairs) == 0:
            return pairs, self.mol, self.cutoff

        print('pairs: ', pairs)
        self.idx2Mol(pairs)
        print('{} bonds are going to be generated'.format(len(pairs)))
        print('Following bonds will be formed: \n')
        for index, value in pairs.iterrows():
            print('\t', value.acro, '\t', value.amon, '\t',
                  round(value.p, 2), '\t', value.rctP)
        return pairs, self.mol, self.cutoff
