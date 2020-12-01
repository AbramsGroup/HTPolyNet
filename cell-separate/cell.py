# separate the cell based on the gro coordinate
# 1. Read atoms from gro file
# 2. Build cells based on the system size
# 3. assign cell id to each atoms df

import pandas as pd
import numpy as np
import os

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

def genCell(boxSize):
    boxSize = '3.000000'
    parts = 2
    x = np.linspace(0, float(boxSize), parts + 1)
    y = np.linspace(0, float(boxSize), parts + 1)
    z = np.linspace(0, float(boxSize), parts + 1)
    xdiv = x[1] - x[0]
    ydiv = y[1] - y[0]
    zdiv = z[1] - z[0]
    cell_id = []
    com = [0.5 * (x[1] - x[0]),
           0.5 * (y[1] - y[0]),
           0.5 * (z[1] - z[0])]
    xNum = 0; yNum = 0; zNum = 0
    for i in range(parts): # x dir
        xCom = com[0] + xdiv * xNum
        xNum += 1; yNum = 0; zNum = 0
        for j in range(parts): # y dir
            yCom = com[1] + ydiv * yNum
            yNum += 1; zNum = 0
            for k in range(parts): # z dir
                zCom = com[2] + zdiv * zNum
                zNum += 1
                info = [i, j, k, [xCom, yCom, zCom]]
                cell_id.append(info)
    return cell_id

if __name__ == '__main__':
    name = os.path.join(os.getcwd(), 'cell-separate/init-1.gro')
    a = readGRO(name)
    df_init, sysName, atNum, boxSize = a.readGRO()
    a = genCell(boxSize.split())
