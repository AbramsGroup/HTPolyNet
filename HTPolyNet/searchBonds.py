# -*- coding: utf-8 -*-
"""
Searching potential bonds

@author: huang
"""
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import networkx as nx
import time
from HTPolyNet.countTime import *
from multiprocessing import Process, Pool
from functools import partial
from copy import deepcopy
import random
import sys

class searchBonds(object):
    def __init__(self, cpu, basicParameters, generatedBonds, inGro, inTop, 
                    conv, desBonds, chains, boxLimit):
        self.cpu = cpu

        self.monInfo = basicParameters.monInfo
        self.croInfo = basicParameters.croInfo
        self.cutoff = basicParameters.cutoff
        self.bondsRatio = basicParameters.bondsRatio
        self.boxSize = basicParameters.boxSize  # three dimensions of box length
        self.maxBonds = basicParameters.maxBonds
        self.rctInfo = basicParameters.rctInfo
        self.HTProcess = basicParameters.HTProcess
        self.rctMap = self.genRctMap(basicParameters.rctInfo)
        self.generatedBonds = generatedBonds
        self.gro = inGro  # gro object
        self.top = inTop  # top object

        self.rctAtoms = rctAtoms = self.gro.df_atoms.loc[(self.gro.df_atoms.rct == 'True') | (self.gro.df_atoms.rct == 'False')].globalIdx.to_list()
        self.rctBonds = self.genRctBondsMap(inGro.df_atoms, inTop.bonds)

        self.start = '0' # use for circle detect
        self.tmpBonds = {}

        self.chains = chains

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

        # layer limit
        self.boxLimit = boxLimit
        self.dimLimit = 1
        self.layerDir = basicParameters.layerDir

    def genRctBondsMap(self, atoms, bonds):
        rctAtoms = self.rctAtoms
        rctBonds = {}
        for index, row in bonds.iterrows():
            if row.ai in rctAtoms and row.aj in rctAtoms:
                if row.ai not in rctBonds.keys():
                    rctBonds[row.ai] = [row.aj]
                else:
                    rctBonds[row.ai].append(row.aj)

                if row.aj not in rctBonds.keys():
                    rctBonds[row.aj] = [row.ai]
                else:
                    rctBonds[row.aj].append(row.ai)

        return rctBonds

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
        boxSize = self.boxSize
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

            if y1 > ylen:
                y1 = y1 - float(boxSize[1])
            elif y1 <= -ylen:
                y1 = y1 + float(boxSize[1])

            if z1 > zlen:
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

        tmpDf = df1.loc[df1.molName.isin(rctName)]
        for i in range(len(rctName)):
            if not tmpDf.empty:
                tmpDf.loc[(tmpDf.molName == rctName[i]), 'rctPct'] = rctPct[i]

        return tmpDf

    def convChainsNum2Name(self, inChain, atomsDf):
        chain = inChain.split(',')
        names = []
        if chain == ['']:
            return ['']
        else:
            for ele in chain:
                n = atomsDf.loc[atomsDf['molNum'] == ele, 'molName'].values[0]
                names.append(n)

            return names

    def calAtnDistance(self, atn, atomDf):
        if self.boxSize[0] == self.boxSize[1] == self.boxSize[2]:
            if self.boxLimit != 1:
                # won't generate bond across y boundary
                atn1 = np.array([atn.posX, atn.posZ]).astype(float)
                atn2Lst = [atomDf.posX.to_list(), atomDf.posZ.to_list()]
                atn2Lst = np.asarray(atn2Lst).T.astype(float)
                delta = np.abs(atn1 - atn2Lst)
                delta = np.where(delta > 0.5 * float(self.boxSize[0]), delta - float(self.boxSize[0]), delta) # assume same box size for all dimension
                
                deltaY = np.abs(np.array(atn.posY).astype(float) - np.asarray(atomDf.posY.to_list()).T.astype(float))
                deltaSum = []
                for i in range(len(deltaY)):
                    tmp = np.insert(delta[i], 1, deltaY[i])
                    deltaSum.append(tmp)
                
                delta = np.asarray(deltaSum)
                dist = np.sqrt((delta ** 2).sum(axis=-1))
                atomDf['dist'] = dist
            else:
                atn1 = np.array([atn.posX, atn.posY, atn.posZ]).astype(float)
                atn2Lst = [atomDf.posX.to_list(), atomDf.posY.to_list(), atomDf.posZ.to_list()]
                atn2Lst = np.asarray(atn2Lst).T.astype(float)
                delta = np.abs(atn1 - atn2Lst)
                delta = np.where(delta > 0.5 * float(self.boxSize[0]), delta - float(self.boxSize[0]), delta) # assume same box size for all dimension
                dist = np.sqrt((delta ** 2).sum(axis=-1))
                atomDf['dist'] = dist
        else:
            dist = []
            atn1 = np.array([atn.posX, atn.posY, atn.posZ]).astype(float)
            atn2Lst = [atomDf.posX.to_list(), atomDf.posY.to_list(), atomDf.posZ.to_list()]
            atn2Lst = np.asarray(atn2Lst).T.astype(float)
            for atn in atn2Lst:
                dist.append(self.calDist(atn1, atn))
            atomDf['dist'] = dist
        return atomDf

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
        if not df2.empty:
            # obtain the distance between the atoms and distance between atoms and the potential atoms within neighbor cell
            # df3 = df2.apply(lambda x: self.appendDist(x, atom, self.boxSize), axis=1)
            df3 = self.calAtnDistance(atom, df2)
            df3 = df3.iloc[1:] # remove atoms which connect to itself
            pd.to_numeric(df3.dist)
            df3.sort_values(by='dist', inplace=True)
            df3.to_csv('df3_sorted.csv')
            # print('cutoff: ', self.cutoff)
            df4 = df3.loc[(df3.dist < float(self.cutoff)) & (df3.molNum != atom.molNum)]
            df4.to_csv('df4.csv') # TODO: append files to check
            df_out1 = []
            if len(df4) == 0:
                return [[]]
            else:
                for index, row in df4.iterrows():
                    df_out1.append([atom.globalIdx, row.globalIdx, round(row.dist, 2), atom.rctNum, row.rctNum,
                                    atom.molCon, random.random(), row.rctPct, atom.molNum, row.molNum, atom.atomName, row.atomName])
                return df_out1
        else:
            return [[]]

    def getAllMolNames(self):
        atomNames = []
        rctFunc = self.rctInfo
        for i in rctFunc:
            if i[0][0] not in atomNames:
                atomNames.append(i[0][0])
            if i[1][0] not in atomNames:
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

    def getHname(self, atomDf):
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
        return filterIdx

    @countTime
    def checkHydrogen(self, atomDf):
        filterIdx = []
        p = Pool(processes=self.cpu)
        atomDf_split = np.array_split(atomDf, self.cpu)
        idxList = p.map(partial(self.getHname), atomDf_split)
        p.close(); p.join()
        for i in idxList:
            filterIdx += i

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

        # select all atoms
        if self.layerDir == 'x':
            df_tmp0 = atoms.loc[(atoms.rct == 'True') & (atoms.molName.isin(atomNames)) & (atoms.posX.astype(float) < self.dimLimit)]
        elif self.layerDir == 'y':
            df_tmp0 = atoms.loc[(atoms.rct == 'True') & (atoms.molName.isin(atomNames)) & (atoms.posY.astype(float) < self.dimLimit)]
        elif self.layerDir == 'z':
            df_tmp0 = atoms.loc[(atoms.rct == 'True') & (atoms.molName.isin(atomNames)) & (atoms.posZ.astype(float) < self.dimLimit)]
        else:
            raise ValueError(f'Wrong direction {self.layerDir}')


        # assign atoms to cell
        df_tmp0 = self.assignAtoms(df_tmp0)
        df_tmp = self.checkHydrogen(df_tmp0) # This check just to confirm selected atom can react
        # ##### START PARALLEL
        print('-----> start parallel searching!!')
        p = Pool(processes=self.cpu)
        dfSplit = np.array_split(df_tmp, self.cpu)
        # TODO: without parallel searching and check the results
        results = p.map(partial(self.parallel_getPairs, df_sum=df_tmp), dfSplit)
        p.close()
        p.join()
        tmpLst = []
        for l in results:
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

        names = ['acro', 'amon', 'dist', 'rctNum1', 'rctNum2', 'molCon', 'p', 'rctP', 'croMol', 'monMol', 'croAtomName', 'monAtomName']
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

    def cycleDetect(self):
        g = nx.Graph()
        g.add_nodes_from(self.tmpBonds.keys())
        for key, value in self.tmpBonds.items():
            g.add_edges_from(([(key, t) for t in value]))

        return nx.cycle_basis(g)

    def searchCycle(self, row):
        atoms = self.rctAtoms
        bonds = self.rctBonds
        tmpBonds = deepcopy(bonds)
        if row.acro in tmpBonds.keys():
            tmpBonds[row.acro].append(row.amon)
        else:
            tmpBonds[row.acro] = [row.amon]

        if row.amon in tmpBonds.keys():
            tmpBonds[row.amon].append(row.acro)
        else:
            tmpBonds[row.amon] = [row.acro]

        with open('rctBonds.txt', 'a') as f:
            f.write('\nOriginal con: ')
            f.write('\n\t{}'.format(self.rctBonds))
            f.write('\nRct pairs: ')
            f.write('\n\t{}\t{}'.format(row.acro, row.amon))
            f.write('\ntmpBonds: ')
            f.write('\n\t{}'.format(tmpBonds))

        self.tmpBonds = tmpBonds
        path = self.cycleDetect()
        if path == []:
            return True
        else:
            with open('cycle.txt', 'a') as f1:
                f1.write('cycle: \n')
                f1.write('\tpath: {}\n'.format(path))
            return False

    def checkCircuit(self, df_rctAtoms, row):
        cond = 0
        mol1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molNum'].to_list()
        mol2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molNum'].to_list()
        rctGrp1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['rctGroup'].to_list()[0]
        rctGrp2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['rctGroup'].to_list()[0]

        molCon1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['molCon'].to_list()[0].split(',')
        molCon2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['molCon'].to_list()[0].split(',')

        # 条件1： 两个分子不能形成环
        if mol1[0] in molCon2 or mol2[0] in molCon1 or mol1[0] == mol2[0]:
            # print('generate circuit, mol {}'.format(mol1[0]))
            keys = self.dumpPairs.keys()
            if row.acro in keys:
                self.dumpPairs[row.acro].append(row.amon)
            else:
                self.dumpPairs[row.acro] = [row.amon]

            if row.amon in keys:
                self.dumpPairs[row.amon].append(row.acro)
            else:
                self.dumpPairs[row.acro] = [row.acro]
            cond += 1
        else:
            pass

        # 条件2： 一条由单体形成的链不能首尾相连
        # chain1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['chain'].to_list()[0]
        # chain2 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.amon]['chain'].to_list()[0]
        # name1 = self.convChainsNum2Name(chain1, df_rctAtoms)
        # name2 = self.convChainsNum2Name(chain2, df_rctAtoms)
        # chain1 = chain1.split(',')
        # chain2 = chain2.split(',')
        # if len(name1) == 1 or len(name2) == 1:
        #     pass
        # else:
        #     croNames = []
        #     for n in self.croInfo:
        #         croNames.append(n[1])

        #     if chain1 == chain2 and any(i not in name1 for i in croNames):
        #         #TODO: add comment here, should be question here
        #         if mol1[0] in chain1 and mol2[0] in chain1:
        #             cond += 1
        #             with open('cond2.txt', 'a') as f:
        #                 f.write('Found STY chain: ')
        #                 f.write('croNames: {}'.format(croNames))
        #                 f.write('\n\tatom1: {}\tatom2: {}'.format(row.acro, row.amon))
        #                 f.write('mol1: {}\tmol2: {}\n'.format(mol1, mol2))
        #                 f.write('\n\tchain: {}'.format(name1))
        #         else:
        #             pass
        t1 = time.time()
        monNum = 0
        for nb in self.monInfo:
            monNum += int(nb[2])

        chain1 = df_rctAtoms[df_rctAtoms.loc[:, 'globalIdx'] == row.acro]['chain'].to_list()[0]
        chain1 = chain1.split(',')

        if len(chain1) == 1:
            pass
        else:
            tmp1 = np.asarray(chain1).astype(int)
            tmp2 = list(np.where(tmp1 > monNum)[0])
            if len(tmp2) > 1:
                pass
            else:
                tmpChain = []
                for i in range(len(chain1)):
                    if i < len(chain1) - 1:
                        tmpChain.append([chain1[i], chain1[i + 1]])

                tmpG = nx.Graph()
                tmpG.add_edges_from(tmpChain)
                newCon = [[mol1[0], mol2[0]]]
                tmpG.add_edges_from(newCon)
                if len(nx.cycle_basis(tmpG)) > 0:
                    cond += 1
                    with open('cond2.txt', 'a') as f:
                        f.write('Found STY chain: ')
                        f.write('\n\tatom1: {}\tatom2: {}'.format(row.acro, row.amon))
                        f.write('mol1: {}\tmol2: {}\n'.format(mol1, mol2))
                else:
                    pass

        t2 = time.time()
        # print('check cycle 2: ', t2 - t1)
        # 条件三：
        # 描述的是一个group中不能形成loop
        # T1-H1--R1--H1-T1, T2-H2--R2--H2-T2
        # 其中T1-H1-T2-H2是不允许形成的
        # 取消的原因：和下面的判断一样，所以可以取消这个参数
        # adj1 = df_rctAtoms[(df_rctAtoms.loc[:, 'molNum'] == mol1[0]) &
        #                    (df_rctAtoms.loc[:, 'rctGroup'] == rctGrp1) &
        #                    (df_rctAtoms.loc[:, 'globalIdx'] != row.acro)]['groupCon']
        # adj2 = df_rctAtoms[(df_rctAtoms.loc[:, 'molNum'] == mol2[0]) &
        #                    (df_rctAtoms.loc[:, 'rctGroup'] == rctGrp2) &
        #                    (df_rctAtoms.loc[:, 'globalIdx'] != row.amon)]['groupCon']

        # if adj1.to_list() == [] and adj2.to_list() == []:
        #     pass

        # else:
        #     grpsCon1 = adj1.to_list()[0].split('-')
        #     grpsCon2 = adj2.to_list()[0].split('-')
        #     info1 = '[{}, {}]'.format(mol1[0], rctGrp1)
        #     info2 = '[{}, {}]'.format(mol2[0], rctGrp2)
        #     if info1 in grpsCon1 or info2 in grpsCon2:
        #         cond += 1
        #         with open('cond3.txt', 'a') as f:
        #             f.write('mol1: {}\tmol2: {}\n'.format(mol1, mol2))
        #     else:
        #         pass

        # 条件4: 不能形成全部由可反应原子形成的环
        if self.searchCycle(row):
            pass
        else:
            cond += 1

        if cond == 0:
            return True
        else:
            return False

    def checkAtomsRepeat(self, inAtoms, lst):
        for a in lst:
            if a in inAtoms:
                return False
            else:
                continue
        return True

    def checkHT(self, at1Idx, at2Idx): # TODO: didn't work
        if self.HTProcess == 'False':
            return True
        else:
            atomsDf = self.gro.df_atoms
            at1Prop = atomsDf.loc[atomsDf.globalIdx == at1Idx].prop.values[0]
            at2Prop = atomsDf.loc[atomsDf.globalIdx == at2Idx].prop.values[0]
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

    def updateChains(self, mol1, mol2):
        tmpChains = {}
        nmol1 = 0; nmol2 = 0
        for i in range(len(self.chains)):
            c = self.chains[i]
            # 同一个分子不会出同时出现在两个链上
            if mol1 in c:
                tmpChains[mol1] = c
                nmol1 += 1
            if mol2 in c:
                tmpChains[mol2] = c
                nmol2 += 1

        if nmol1 > 0 and nmol2 > 0:
            if tmpChains[mol1] == tmpChains[mol2]:
                pass
            else:
                try:
                    self.chains.remove(tmpChains[mol1])
                    self.chains.remove(tmpChains[mol2])
                    new = tmpChains[mol1] + tmpChains[mol2]
                    self.chains.append(new)
                except:
                    print('Unexpected chains')

        elif nmol1 == 0 and nmol2 > 0:
            idx = self.chains.index(tmpChains[mol2])
            self.chains[idx].append(mol1)

        elif nmol2 == 0 and nmol1 > 0:
            idx = self.chains.index(tmpChains[mol1])
            self.chains[idx].append(mol2)

        elif nmol1 == 0 and nmol2 == 0:
            new = [mol1, mol2]
            self.chains.append(new)

    @countTime
    def updateAllChains(self):
        df_atoms = self.gro.df_atoms
        chains = self.chains
        for c in chains:
            tmpStr = ','.join(c)
            df_atoms.loc[df_atoms['molNum'].isin(c), 'chain'] = tmpStr

    def updateGroupCon(self, a1, a2, mol1, mol2, atomsDf):
        # 已经在找getPairs的一步判定了该连接方式是否征程
        group1 = atomsDf[atomsDf.globalIdx == str(a1)].rctGroup.values[0]
        group2 = atomsDf[atomsDf.globalIdx == str(a2)].rctGroup.values[0]

        grp1 = atomsDf.loc[atomsDf['globalIdx'] == str(a1), 'groupCon'].values[0]
        grp2 = atomsDf.loc[atomsDf['globalIdx'] == str(a2), 'groupCon'].values[0]
        if grp1 == '':
            atomsDf.loc[atomsDf['globalIdx'] == str(a1), 'groupCon'] = '[{}, {}]'.format(mol2, group2)
        else:
            tmpStr = grp1 + '-[{}, {}]'.format(mol2, group2)
            atomsDf.loc[atomsDf['globalIdx'] == str(a1), 'groupCon'] = tmpStr

        if grp2 == '':
            atomsDf.loc[atomsDf['globalIdx'] == str(a2), 'groupCon'] = '[{}, {}]'.format(mol1, group1)
        else:
            tmpStr = grp2 + '-[{}, {}]'.format(mol1, group1)
            atomsDf.loc[atomsDf['globalIdx'] == str(a2), 'groupCon'] = tmpStr

    def updateMolCon(self, atomsDf, a1, a2):
        mol1 = atomsDf[atomsDf.globalIdx == str(a1)].molNum.values[0]
        mol2 = atomsDf[atomsDf.globalIdx == str(a2)].molNum.values[0]

        con1 = list(atomsDf[atomsDf.molNum == mol1].molCon)[0]
        con2 = list(atomsDf[atomsDf.molNum == mol2].molCon)[0]

        tmp = '{},{}'.format(con1, mol2)
        atomsDf.loc[atomsDf['molNum'] == mol1, 'molCon'] = tmp

        tmp = '{},{}'.format(con2, mol1)
        atomsDf.loc[atomsDf['molNum'] == mol2, 'molCon'] = tmp

        self.updateChains(mol1, mol2)
        # print('Bonds formed between mol {} and mol {}'.format(mol1, mol2))
        self.updateGroupCon(a1, a2, mol1, mol2, atomsDf)

        # update tmp bond session
        if a1 in self.rctBonds.keys():
            self.rctBonds[a1].append(a2)
        else:
            self.rctBonds[a1] = [a2]

        if a2 in self.rctBonds.keys():
            self.rctBonds[a2].append(a1)
        else:
            self.rctBonds[a2] = [a1]

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

        print('----> {} bonds goes to check reactivity probability'.format(len(df_tmp0)))
        pcriteria = 0
        if len(df_tmp0) <= 2:
            pcriteria = 1

        cc = 0
        while(len(rowList1) == 0 and cc < 5):
            for index, row in df_tmp0.iterrows():
                if self.checkKineticRatio(pcriteria, row):
                    rowList1.append(row)

            if self.conv < 0.7:
                if len(rowList1) < 0.2 * self.desBonds:
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
        df_tmp2 = self.updateAllCon(df_tmp1)

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
        if self.cutoff > float(self.boxSize[0]) * 0.72:
            print('cutoff should smaller than the half of the box size, please modify it.')
            sys.exit()

        while (xdiv <= self.cutoff):
            if parts <= 1:
                print('Unusual occasion occurs! Please check the cutoff you set and rerun the script!')
                sys.exit()

            x = np.linspace(0, float(self.boxSize[0]), parts + 1)
            y = np.linspace(0, float(self.boxSize[1]), parts + 1)
            z = np.linspace(0, float(self.boxSize[2]), parts + 1)
            xdiv = x[1] - x[0]
            ydiv = y[1] - y[0]
            zdiv = z[1] - z[0]
            parts -= 1

        cell_id = []
        div_box = [xdiv, ydiv, zdiv]
        com = [0.5 * (x[1] - x[0]),
               0.5 * (y[1] - y[0]),
               0.5 * (z[1] - z[0])]
        xNum = 0
        yNum = 0
        zNum = 0
        xMax = 0
        yMax = 0
        zMax = 0

        for i in range(parts):  # x dir
            xCom = com[0] + xdiv * xNum
            yMax = yNum
            zMax = zNum
            xNum += 1
            yNum = 0
            zNum = 0
            for j in range(parts):  # y dir
                yCom = com[1] + ydiv * yNum
                yNum += 1
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
        # assign all unrcted atoms to different cell
        df1 = df_atoms.apply(lambda x: self.searchCell(x), axis=1)
        return df1

    def updateAllCon(self, df_tmp):
        atomsDf = self.gro.df_atoms
        rowList = []
        atomsList = []
        df_tmp.to_csv('tmp1.csv')
        for index, row in df_tmp.iterrows():
            if self.checkCircuit(atomsDf, row):
                if self.checkAtomsRepeat(atomsList, [row.acro, row.amon]):
                    atomsList.append(row.acro)
                    atomsList.append(row.amon)
                    rowList.append(row)
                    self.updateMolCon(atomsDf, row.acro, row.amon)
            else:
                continue

        df_out = pd.DataFrame(rowList)
        df_out.to_csv('tmp2.csv')
        return df_out
    
    @countTime
    def splitLayer(self):
        layerLst = []
        molNum = self.gro.df_atoms['molNum'].to_list()
        for nm in molNum:
            tmp_mol = self.gro.df_atoms.loc[(self.gro.df_atoms.molNum == nm)]
            if self.layerDir == 'x':
                tmp_layer_atns = tmp_mol.loc[(tmp_mol.posX.astype(float) < self.dimLimit)]
            elif self.layerDir == 'y':
                tmp_layer_atns = tmp_mol.loc[(tmp_mol.posY.astype(float) < self.dimLimit)]
            elif self.layerDir == 'z':
                tmp_layer_atns = tmp_mol.loc[(tmp_mol.posZ.astype(float) < self.dimLimit)]
            else:
                raise ValueError(f'Wrong direction is input, {self.layerDir}')
            if len(tmp_layer_atns) > int(0.8 * len(tmp_mol)):
                layerLst.append(nm)
        return set(layerLst)

    @countTime
    def collectBonds(self, count):
        pairs = []
        count = 0
        largest_dim = max(self.boxSize)
        # df_layer = self.gro.df_atoms.loc[(self.gro.df_atoms.posY.astype(float) < self.dimLimit)]
        # layerLst = list(set(df_layer['molNum'].to_list()))
        layerLst = self.splitLayer()
        print('layerLst: ', layerLst)
        with open('molLayer.txt', 'w') as f:
            for i in layerLst:
                f.write('{}\t'.format(i))

        while (len(pairs) == 0 and count < 10):
            while (len(pairs) == 0):
                print('-----> Current cutoff: {}'.format(round(self.cutoff, 2)))
                
                df_pairs = self.getRctDf()
                df_pairs.to_csv('all_bonds_within_cutoff.csv')
                if self.conv < 0.7:
                    if len(df_pairs) == 0 or len(df_pairs) < 0.2 * self.desBonds:
                        print('-----> Not enough pairs found in this step, increasing cutoff by 0.2nm')
                        self.cutoff += 0.2
                        if self.cutoff > 0.7 * float(largest_dim):
                            break
                        else:
                            continue
                    else:
                        break

                elif 0.7 <= self.conv <= 0.9:
                    if len(df_pairs) == 0 or len(df_pairs) < 0.2 * self.desBonds:
                        print('-----> Not enough pairs found in this step, increasing cutoff by 0.2nm')
                        self.cutoff += 0.2
                        if self.cutoff > 0.75 * float(largest_dim):
                            break
                        else:
                            continue
                    else:
                        break
                else:
                    if len(df_pairs) == 0:
                        print('-----> Not enough pairs found in this step, increasing cutoff by 0.1nm')
                        self.cutoff += 0.1
                        if self.cutoff > 0.6 * float(largest_dim):
                            break
                        else:
                            continue
                    else:
                        break

            df_pairs = self.finalRctPairs(df_pairs)
            pairs = df_pairs
            count += 1
            if len(pairs) == 0:
                self.cutoff += 0.2
            with open('count.txt', 'a') as f:
                f.write('\ncount: {}'.format(count))
        return pairs

    def updateBoxSize(self):
        self.boxSize = [float(self.gro.boxSize.split()[0]), float(self.gro.boxSize.split()[1]), float(self.gro.boxSize.split()[2])]
        print('self.boxSize: ', self.boxSize)
        # two occasions, 1. box size is same along all directions; 2. box size is different along one direction
        if self.layerDir == 'x':
            self.dimLimit = float(self.boxSize[0]) * self.boxLimit
        elif self.layerDir == 'y':
            self.dimLimit = float(self.boxSize[1]) * self.boxLimit
        elif self.layerDir == 'z':
            self.dimLimit = float(self.boxSize[2]) * self.boxLimit

        print('self.boxLimit: ', self.boxLimit)
        print('self.dimLimit: ', self.dimLimit)

    @countTime
    def sBonds(self):
        count = 5
        parts = 8
        self.updateBoxSize()
        self.genCell(parts)
        pairs = self.collectBonds(count)
        if len(pairs) == 0:
            return pairs, self.chains, self.mol, self.cutoff

        self.idx2Mol(pairs)

        # update all atoms belonged chain
        self.updateAllChains()
        print('-----> {} bonds are going to be generated'.format(len(pairs)))
        # print('Following bonds will be formed: \n')
        # for index, value in pairs.iterrows():
        #     print('\t', value.acro, '\t', value.amon, '\t',
        #           round(value.p, 2), '\t', value.rctP)
        return pairs, self.chains, self.mol, self.cutoff