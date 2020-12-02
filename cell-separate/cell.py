# separate the cell based on the gro coordinate
# 1. Read atoms from gro file
# 2. Build cells based on the system size
# 3. assign cell id to each atoms df

import pandas as pd
import numpy as np
import os
import sys

class cell(object):
    def __init__(self):
        self.id = ''
        self.center = ''
        self.length = ''

class readGRO(object):
    def __init__(self, name):
        self.groName = name
        self.atomsDF = ''

    def splitRow(self, row, idx=1):
        x = row
        molNum = x[0:5].strip(' ')
        molName = x[5:10].strip(' ')
        atomName = x[10:15].strip(' ')
        globalIdx = x[15:20].strip(' ')
        posX = x[21:28].strip(' ')
        posY = x[29:36].strip(' ')
        posZ = x[37:44].strip(' ')
        data = [molNum, molName, atomName, globalIdx, posX, posY, posZ]
        if idx == 0:
            return data[0]
        elif idx == 1:
            return data[1]
        elif idx == 2:
            return data[2]
        elif idx == 3:
            return data[3]
        elif idx == 4:
            return data[4]
        elif idx == 5:
            return data[5]
        elif idx == 6:
            return data[6]
        else:
            sys.exit()

    def readGRO(self):
        name = self.groName
        columns_name = [0, 'molNum', 'molName', 'atomName', 'globalIdx', 'x', 'y', 'z']
        df = pd.read_csv(name, names=columns_name, header=None, sep='\n', skip_blank_lines=False)
        sysName = df.loc[0][0];
        atNum = df.loc[1][0];
        boxSize = df.iloc[-1][0]
        df1 = df.iloc[2:-1]
        df1.loc[:, 'molNum'] = df1[0].apply(lambda x: self.splitRow(x, idx=0))
        df1.loc[:, 'molName'] = df1[0].apply(lambda x: self.splitRow(x, idx=1))
        df1.loc[:, 'atomName'] = df1[0].apply(lambda x: self.splitRow(x, idx=2))
        df1.loc[:, 'globalIdx'] = df1[0].apply(lambda x: self.splitRow(x, idx=3))
        df1.loc[:, 'posX'] = df1[0].apply(lambda x: self.splitRow(x, idx=4))
        df1.loc[:, 'posY'] = df1[0].apply(lambda x: self.splitRow(x, idx=5))
        df1.loc[:, 'posZ'] = df1[0].apply(lambda x: self.splitRow(x, idx=6))
        df_init = df1[['molNum', 'molName', 'atomName', 'globalIdx', 'posX', 'posY', 'posZ']]
        df_init.reset_index(drop=True, inplace=True)
        self.atomsDF = df_init
        return df_init, sysName, atNum, boxSize

def getId(x, maxId):
    outList = []
    tmpLst = [-1, 0, 1]
    for i in tmpLst:
        tmp = x + i
        if tmp < 0:
            tmp = maxId
        elif tmp > maxId:
            tmp = 0
        else:
            pass
        outList.append(tmp)
    return outList

def filterAtoms(atoms, atomsDf, maxCellId): # collect atoms based on cell id. itself and adjacent cell
    # maxCellId used for pdb condition
    cell0 = atoms.cellId
    tmpLst = [-1, 0, 1]
    xList = getId(cell0[0], maxCellId[0])
    yList = getId(cell0[1], maxCellId[1])
    zList = getId(cell0[2], maxCellId[2])

    cellSum = []
    for i in xList:
        for ii in yList:
            for iii in zList:
                cellSum.append([i, ii, iii])

    print('cellSum: ', cellSum)
    print('cellId： ', cell0)
    df_out = atomsDf.loc[atomsDf.cellId.to_list().isin(cellSum)]
    # return df_out

def searchAtoms(atomsDf, atomNames):
    df = atomsDf.loc[atomsDf.atomName.isin(atomNames)]
    return df

def genCell(boxSize):
    # boxSize = '3.000000'
    parts = 5
    x = np.linspace(0, float(boxSize[0]), parts + 1)
    y = np.linspace(0, float(boxSize[1]), parts + 1)
    z = np.linspace(0, float(boxSize[2]), parts + 1)
    xdiv = x[1] - x[0]
    ydiv = y[1] - y[0]
    zdiv = z[1] - z[0]
    cell_id = []; div_box = [xdiv, ydiv, zdiv]
    com = [0.5 * (x[1] - x[0]),
           0.5 * (y[1] - y[0]),
           0.5 * (z[1] - z[0])]
    xNum = 0; yNum = 0; zNum = 0
    xMax = 0; yMax = 0; zMax = 0

    for i in range(parts): # x dir
        xCom = com[0] + xdiv * xNum
        yMax = yNum; zMax = zNum
        xNum += 1; yNum = 0; zNum = 0
        for j in range(parts): # y dir
            yCom = com[1] + ydiv * yNum
            yNum += 1; zNum = 0
            for k in range(parts): # z dir
                zCom = com[2] + zdiv * zNum
                zNum += 1
                info = [[i, j, k], [xCom, yCom, zCom]] # [[id], [center coord]]
                cell_id.append(info)
        xMax = xNum
        maxCellId = [xMax - 1, yMax - 1, zMax - 1] # Since it starts from 0
    return cell_id, div_box, maxCellId


def searchCell(row, cell_id, box_div):
    xDiv = 0.5 * box_div[0]
    yDiv = 0.5 * box_div[1]
    zDiv = 0.5 * box_div[2]
    row['cellId'] = '0'
    row['comCoord'] = '0'
    coord = [float(row.posX), float(row.posY), float(row.posZ)]
    for c in cell_id:
        xlow, xhigh = [c[1][0] - xDiv, c[1][0] + xDiv]
        ylow, yhigh = [c[1][1] - yDiv, c[1][1] + yDiv]
        zlow, zhigh = [c[1][2] - zDiv, c[1][2] + zDiv]
        if coord[0] >= xlow and coord[0] <= xhigh and coord[1] >= ylow and coord[1] <= yhigh and coord[2] >=zlow and coord[2] <= zhigh:
            row['cellId'] = c[0]
            row['comCoord'] = c[1]
            return row
        else:
            continue

def assignAtoms(atomsDf, cell_id, box_div):
    df1 = atomsDf.apply(lambda x: searchCell(x, cell_id, box_div), axis=1)
    return df1

if __name__ == '__main__':
    name = os.path.join(os.getcwd(), 'cell-separate/init-1.gro')
    a = readGRO(name)
    df_init, sysName, atNum, boxSize = a.readGRO()
    cell_id, div_box, maxCellId = genCell(boxSize.split())
    df = assignAtoms(df_init, cell_id, div_box)
    print(df.head())

    atomNames = ['C1']
    df2 = searchAtoms(df, atomNames)
    df3 = filterAtoms(df2.iloc[1], df2, maxCellId)