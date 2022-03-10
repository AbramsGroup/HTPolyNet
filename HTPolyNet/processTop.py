import pandas as pd
import topInfo
from shutil import copyfile
import os
class processTop(object):
    def __init__(self, name, repeat=False):
        self.name = name
        self.repeat = repeat
        if repeat:
            self.topName = '{}.top-bk'.format(name)
        else:
            self.topName = '{}.top'.format(name)
        self.itpName = '{}.itp'.format(name)

        self.aTypes = ''
        self.mTypes = ''
        self.bTypes = ''
        self.angTypes = ''
        self.dihTypes = ''
        self.impTypes = ''

        self.atoms = ''
        self.bonds = ''
        self.pairs = ''
        self.angles = ''
        self.dihs = ''
        self.imps = ''
        self.system = ''
        self.molecules = ''
        self.dupDihTypeKey = []
        self.topInfo = ''
        self.sumTop = []

        self.top = topInfo.top()
    def checkFile(self):
        # check if top has been processed
        if os.path.isfile(self.topName) and os.path.isfile(self.itpName):
            return True
        else:
            return False

    def getTopInfo(self):
        # read old top file, and move it to the top-bk file
        df1 = pd.read_csv(self.topName, names=['0'], comment=';', header=None, sep='\n', skip_blank_lines=True)
        dil_indx = list(df1.loc[df1['0'].str.startswith('[')].index)
        df_sep = []
        for i in range(len(dil_indx)):
            if i == 0:
                continue
            else:
                df_tmp = df1.iloc[dil_indx[i - 1] + 1:dil_indx[i], :]
                if '#' in df_tmp.to_string():
                    continue
                else:
                    df_sep.append(df_tmp)
        df_sep.append(df1.iloc[dil_indx[i] + 1:, :])
        copyfile(self.topName, '{}-bk'.format(self.topName))
        return df_sep

    def rmDihType(self, df,
                  dupKey=False):  # Since some dih type are mulitple defined, they don't need to add to the dih type section. Detail parameters are list in the dih section
        keys = [];
        dup_keys = []
        if dupKey:
            '''
            In this option, passed df is the df_dih
            '''
            for index, row in df.iterrows():
                if row.c0 != ' ':
                    a1Type, a2Type, a3Type, a4Type = self.subAtom2Atypes(row, self.atoms, length=4, search=True)
                    key = '{}-{}-{}-{}'.format(a1Type, a2Type, a3Type, a4Type)
                    dup_keys.append(key)
            self.dupDihTypeKey = dup_keys
        else:
            '''
            In this option, passed df is the df_dihTypes
            '''
            for index, row in df.iterrows():
                key = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
                if key in keys:
                    dup_keys.append(key)
                else:
                    keys.append(key)
            self.dupDihTypeKey = dup_keys

            for index, row in df.iterrows():
                key = '{}-{}-{}-{}'.format(row.ai, row.aj, row.ak, row.al)
                if key in dup_keys:
                    df.drop(index, inplace=True)
        return df

    def subAtom2Atypes(self, row, df_atoms, length=2, search=False):
        if length == 2:
            a1 = row.ai; a2 = row.aj
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']
            row.ai = a1Type; row.aj = a2Type
        elif length == 3:
            a1 = row.ai; a2 = row.aj; a3 = row.ak
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']
            a3Type = df_atoms.loc[int(a3)-1, 'type']
            row.ai = a1Type; row.aj = a2Type; row.ak = a3Type
        elif length == 4:
            a1 = row.ai; a2 = row.aj; a3 = row.ak; a4 = row.al
            a1Type = df_atoms.loc[int(a1)-1, 'type']
            a2Type = df_atoms.loc[int(a2)-1, 'type']
            a3Type = df_atoms.loc[int(a3)-1, 'type']
            a4Type = df_atoms.loc[int(a4)-1, 'type']
            if search:
                return [a1Type, a2Type, a3Type, a4Type]
            else:
                row.ai = a1Type; row.aj = a2Type; row.ak = a3Type; row.al = a4Type
        return row

    def extractType(self, df, df_atoms, keys='bonds'):
        if keys == 'bonds':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=2), axis=1)
        elif keys == 'angles':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=3), axis=1)
        elif keys == 'dih':
            df_new = df.copy(deep=True)
            df_new = df_new.apply(lambda x: self.subAtom2Atypes(x, df_atoms, length=4), axis=1)
        return df_new

    def sepData(self, row, length, idx=0):
        data = list(row.str.split())[0]
        if len(data) != length:
            try:
                a = data[idx]
                return a
            except:
                return ' '
            return row.to_string()
        else:
            return data[idx]

    def initSession(self, df_new, df_ori, cNames):
        newList = []
        for index, value in df_ori.iterrows():
            row = list(value.str.split())[0]
            newList.append(row)

        df_new = pd.DataFrame(newList, columns=cNames)
        # i = 0
        # for c in cNames:
        #     df_new.loc[:, c] = df_ori.apply(lambda x: self.sepData(x, len(cNames), idx=i), axis=1)
        #     i += 1
        return df_new

    def genTop(self, inLst):
        '''
        :param inLst: inLst is from the getTopInfo function.
        We assume top file is from the parmed.gromacs. Data is directly from the
        Antechamber. So the structure of the top file is fixed.
        The unprocessed top file sequence as follow:
        0. default
        1. atomstypes
        2. moleculetype
        3. atoms
        4. bonds
        5. pairs
        6. angles
        7. dihedrals
        8. impropers (not always appeared)
        9. system
        10.molecules
        Section 1, 10, 11 write in the top file
        Section 2 - 9 write in the new itp file
        Generate new sections of different types (bonds, angles, dih, imps) besides atomtype section

        Need to add a check on this sequence to promise this always work.
        Usually this won't be a matter unless parmed did some upgrade and change their writing structure
        '''
        #
        top0 = topInfo.top()
        default = inLst[0]
        system = inLst[-2]
        molecules = inLst[-1]
        

        atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
        moltypeNames = ['name', 'nrexcl']
        atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
        bNames = ['ai', 'aj', 'funct', 'c0', 'c1']
        pNames = ['ai', 'aj', 'funct']
        angNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
        dihNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        impNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
        df_atypes = pd.DataFrame(columns=atypeNames)
        df_mtypes = pd.DataFrame(columns=moltypeNames)
        df_atoms = pd.DataFrame(columns=atNames)
        df_bonds = pd.DataFrame(columns=bNames)
        df_pairs = pd.DataFrame(columns=pNames)
        df_angles= pd.DataFrame(columns=angNames)
        df_dih = pd.DataFrame(columns=dihNames)
        df_imp = pd.DataFrame(columns=impNames)

        df_lst = inLst[1:-2]
        df_atypes = self.initSession(df_atypes, df_lst[0], atypeNames).reset_index(drop=True)
        df_mtypes = self.initSession(df_mtypes, df_lst[1], moltypeNames).reset_index(drop=True)
        df_atoms = self.initSession(df_atoms, df_lst[2], atNames).reset_index(drop=True)
        df_bonds = self.initSession(df_bonds, df_lst[3], bNames).reset_index(drop=True)
        df_pairs = self.initSession(df_pairs, df_lst[4], pNames).reset_index(drop=True)
        df_angles = self.initSession(df_angles, df_lst[5], angNames).reset_index(drop=True)
        df_dih = self.initSession(df_dih, df_lst[6], dihNames).reset_index(drop=True)
        df_bTypes = self.extractType(df_bonds, df_atoms, keys='bonds').drop_duplicates().reset_index(drop=True)
        df_angTypes = self.extractType(df_angles, df_atoms, keys='angles').drop_duplicates().reset_index(drop=True)
        df_dihTypes = self.extractType(df_dih, df_atoms, keys='dih').drop_duplicates().reset_index(drop=True)
        df_impTypes = pd.DataFrame(impNames)

        df_dihTypes = self.rmDihType(df_dihTypes)
        if len(df_lst) == 8:
            df_imp = self.initSession(df_imp, df_lst[7], dihNames).reset_index(drop=True)
            df_impTypes = self.extractType(df_imp, df_atoms, keys='dih').drop_duplicates().reset_index(drop=True)

        else:
            df_imp = df_imp
            df_impTypes = pd.DataFrame(impNames)

        self.aTypes = df_atypes
        self.mTypes = df_mtypes
        self.bTypes = df_bTypes
        self.angTypes = df_angTypes
        self.dihTypes = df_dihTypes
        self.impTypes = df_impTypes

        self.atoms = df_atoms
        self.bonds = df_bonds
        self.pairs = df_pairs
        self.angles = df_angles
        self.dihs = df_dih
        self.imps = df_imp

        top0.default = pd.DataFrame(default)
        self.topInfo = [system, molecules]

        self.sumTop = [self.topInfo, self.aTypes, self.mTypes, self.bTypes, self.angTypes, self.dihTypes, self.impTypes,
                       self.atoms, self.bonds, self.pairs, self.angles, self.dihs, self.imps, self.dupDihTypeKey]

        top0.setInfo(self.sumTop)
        # if self.repeat:
        top0.outDf(self.name)

        self.top = top0

    def main(self):
        lst = self.getTopInfo()
        self.genTop(lst)

if __name__ == '__main__':
    a = processTop('STY')
    a1 = a.getTopInfo()
    a.genTop(a1)
