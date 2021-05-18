import re
import sys

class endCapping(object):
    def __init__(self, inGro, inTop, inCat):
        self.gro = inGro
        self.top = inTop
        self.cat = inCat
        if self.cat == 'VE-ST':
            pass

        else:
            pass

    def changeAtypes(self, pPairs):
        chargeCE = '0.014275'
        chargeC2 = '-0.090042'
        df_atoms = self.top.atoms
        for p in pPairs:
            if int(float(df_atoms.loc[(df_atoms.nr == p[0]), 'mass'])) == 12:
                hAtoms = self.findHydrogen(p[0])
                if len(hAtoms) > 0:
                    df_atoms.loc[(df_atoms.nr == p[0]), 'type'] = 'c2'
                    df_atoms.loc[(df_atoms.nr == p[0]), 'charge'] = chargeC2
                else:
                    df_atoms.loc[(df_atoms.nr == p[0]), 'type'] = 'ce'
                    df_atoms.loc[(df_atoms.nr == p[0]), 'charge'] = chargeCE

            else:
                print('Found wrong atoms during end capping: \n\t', df_atoms)

            if int(float(df_atoms.loc[(df_atoms.nr == p[1]), 'mass'])) == 12:
                hAtoms = self.findHydrogen(p[1])
                if len(hAtoms) > 0:
                    df_atoms.loc[(df_atoms.nr == p[1]), 'type'] = 'c2'
                    df_atoms.loc[(df_atoms.nr == p[1]), 'charge'] = chargeC2
                else:
                    df_atoms.loc[(df_atoms.nr == p[1]), 'type'] = 'ce'
                    df_atoms.loc[(df_atoms.nr == p[1]), 'charge'] = chargeCE
            else:
                print('Found wrong atoms during end capping: \n\t', df_atoms)


    def getPairs(self):
        pPairs = []
        pAtoms = self.gro.loc[self.gro.rct == 'True']
        for index, value in pAtoms.iterrows():
            a1GlobalIdx = pAtoms.loc[index, 'globalIdx']
            a1MolNum = pAtoms.loc[index, 'molNum']
            a1Group = pAtoms.loc[index, 'rctGroup']
            for i, v in pAtoms.iterrows():
                a2GlobalIdx = pAtoms.loc[i, 'globalIdx']
                a2MolNum = pAtoms.loc[i, 'molNum']
                a2Group = pAtoms.loc[i, 'rctGroup']
                if a1GlobalIdx == a2GlobalIdx:
                    continue
                else:
                    if a1MolNum == a2MolNum and a1Group == a2Group:
                        pPairs.append([a1GlobalIdx, a2GlobalIdx])
        return pPairs

    def searchCon(self, idx, df_bonds, df_new=[]):
        idx = str(idx)
        df_out = df_bonds[(df_bonds.ai == idx) | (df_bonds.aj == idx)] # df_bonds is the bond section in the topology
        con = []
        for index, row in df_out.iterrows():
            if row.ai == str(idx):
                con.append(row.aj)
            elif row.aj == str(idx):
                con.append(row.ai)
            else:
                print('bond connection met error, please check!')
                print('df_out: ', df_out)
                sys.exit()

        for i in df_new:
            if i[0] == str(idx):
                con.append(i[1])
            elif i[1] == str(idx):
                con.append(i[0])
            else:
                pass
        return con

    def findHydrogen(self, idx): # This idx is the pd index
        atomsDf = self.gro.df_atoms
        bonds = self.top.bonds
        con = []
        conSum = self.searchCon(idx, bonds)
        for i in conSum:
            if 'H' in atomsDf[atomsDf.loc[:, 'globalIdx'] == str(i)]['atomName'].values[0]:
                con.append(i)

        # sort the hydrogen based on the atom name
        con1 = []; atomName_ori = ''
        for idx in con:
            atName = atomsDf[atomsDf.loc[:, 'globalIdx'] == idx]['atomName'].values[0]
            if len(con1) == 0:
                con1 = [idx]
                atomName_ori = atName
            else:
                num0 = re.findall(r'\d+', atomName_ori)
                num1 = re.findall(r'\d+', atName)
                if len(num0) == 0:
                    break
                else:
                    if len(num1) == 0:
                        con1 = [idx]
                        atomName_ori = atName
                    elif int(num1[0]) < int(num0[0]):
                        con1 = [idx]
                        atomName_ori = atName
                    elif int(num1[0]) > int(num0[0]):
                        continue
        return con1

    def delHydrogen(self, pairs):
        '''
        1. Generate new bonds/pairs/angles/dihs
        2. Check new types if exist in the database
        3. Find and delete hydrogen
        '''
        atomsDf = self.gro.df_atoms
        inTop = self.top
        df_atoms = inTop.atoms  # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers
        hAtoms = []
        for p in pairs:
            hCon1 = self.findHydrogen(atomsDf, df_bonds, p[0])[0]
            hCon2 = self.findHydrogen(atomsDf, df_bonds, p[1])[0]
            hAtoms.append(hCon1);
            hAtoms.append(hCon2)

        for a in hAtoms:
            atomsDf.drop(atomsDf[atomsDf['globalIdx'] == str(a)].index, inplace=True)
            self.gro.df_atoms = atomsDf
            df_atoms.drop(df_atoms[df_atoms['nr'] == a].index, inplace=True)
            df_bonds.drop(df_bonds[df_bonds['ai'] == a].index, inplace=True)
            df_bonds.drop(df_bonds[df_bonds['aj'] == a].index, inplace=True)
            df_pairs.drop(df_pairs[df_pairs['ai'] == a].index, inplace=True)
            df_pairs.drop(df_pairs[df_pairs['aj'] == a].index, inplace=True)
            df_angs.drop(df_angs[df_angs['ai'] == a].index, inplace=True)
            df_angs.drop(df_angs[df_angs['aj'] == a].index, inplace=True)
            df_angs.drop(df_angs[df_angs['ak'] == a].index, inplace=True)
            df_dihs.drop(df_dihs[df_dihs['ai'] == a].index, inplace=True)
            df_dihs.drop(df_dihs[df_dihs['aj'] == a].index, inplace=True)
            df_dihs.drop(df_dihs[df_dihs['ak'] == a].index, inplace=True)
            df_dihs.drop(df_dihs[df_dihs['al'] == a].index, inplace=True)
            df_imps.drop(df_imps[df_imps['ai'] == a].index, inplace=True)
            df_imps.drop(df_imps[df_imps['aj'] == a].index, inplace=True)
            df_imps.drop(df_imps[df_imps['ak'] == a].index, inplace=True)
            df_imps.drop(df_imps[df_imps['al'] == a].index, inplace=True)

        self.top.atoms = df_atoms
        self.top.bonds = df_bonds
        self.top.pairs = df_pairs
        self.top.angs = df_angs
        self.top.dihedrals = df_dihs
        self.top.impropers = df_imps

    def updateTopIdx_cfa(self,x,myd,types='atoms'):
        if types == 'atoms':
            if x.nr not in myd.keys():
                print('Unknown atoms Id: ', x.nr)
                sys.exit()
            newidx1 = myd[x.nr];  x.new_idx = newidx1; x.nr = newidx1; x.cgnr = newidx1

        elif types == 'bonds':
            newidx1 = myd[x.ai]
            newidx2 = myd[x.aj]
            x.ai = str(newidx1)
            x.aj = str(newidx2)

        elif types == 'pairs':
            newidx1 = myd[x.ai]
            newidx2 = myd[x.aj]
            x.ai = str(newidx1)
            x.aj = str(newidx2)

        elif types == 'angles':
            newidx1 = myd[x.ai]
            newidx2 = myd[x.aj]
            newidx3 = myd[x.ak]
            x.ai = str(newidx1)
            x.aj = str(newidx2)
            x.ak = str(newidx3)

        elif types == 'dih':
            newidx1 = myd[x.ai]
            newidx2 = myd[x.aj]
            newidx3 = myd[x.ak]
            newidx4 = myd[x.al]
            x.ai = str(newidx1)
            x.aj = str(newidx2)
            x.ak = str(newidx3)
            x.al = str(newidx4)

        else:
            print('sth wrong')
        return x

    def mapUpdate(self, df_atoms, newIDx_from_oldIdx):
        # Using dataframe structure to build
        inTop = self.top
        # dataframes from the topology
        df_atoms = inTop.atoms  # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs  # 1--4 interactions
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers

        df_atoms_new = df_atoms.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='atoms'), axis=1).reset_index(
            drop=True)
        df_bonds_new = df_bonds.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='bonds'), axis=1)
        df_pairs_new = df_pairs.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='pairs'), axis=1)
        df_angs_new = df_angs.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='angles'), axis=1)
        df_dihs_new = df_dihs.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='dih'), axis=1)
        df_imps_new = df_imps.apply(lambda x: self.updateTopIdx_cfa(x, newIDx_from_oldIdx, types='dih'), axis=1)
        return df_atoms_new, df_bonds_new, df_pairs_new, df_angs_new, df_dihs_new, df_imps_new

    def updateIdx(self):
        # atomsDf is a from the gro file (coordinates)
        atomsDf = self.gro.df_atoms

        inTop = self.top
        # dataframes from the topology
        df_atoms = inTop.atoms  # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs  # 1--4 interactions
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers

        # globalIdx is old; store them in new column 'ori_idx'
        atomsDf['ori_idx'] = atomsDf['globalIdx']
        # creating a new index column that renumbers atoms
        atomsDf = atomsDf.reset_index(drop=True)
        # setting new globalIdx to be 1-initiated indices, as *strings*
        atomsDf['globalIdx'] = (atomsDf.index + 1).astype(str).to_list()
        # make old-to-new index dictionary (cfa)
        newIDx_from_oldIdx = {}
        for o, n in zip(atomsDf['ori_idx'], atomsDf['globalIdx']):
            newIDx_from_oldIdx[o] = n

        # create new dataframe that applies new globalIdx values to old values in original top df's
        df_atoms_new, df_bonds_new, df_pairs_new, df_angs_new, df_dihs_new, df_imps_new = \
            self.mapUpdate(df_atoms, newIDx_from_oldIdx)

        inTop.atoms = df_atoms_new
        inTop.bonds = df_bonds_new
        inTop.pairs = df_pairs_new
        inTop.angles = df_angs_new
        inTop.dihedrals = df_dihs_new
        inTop.impropers = df_imps_new
        self.gro.df_atoms = atomsDf;
        self.gro.atNum = len(atomsDf)

    def VECapping(self):
        # 1. Two atoms belong to the same molecule and same group --> unreact vinyl group
        # 2. change atom types
        # 3. remove h atoms
        # 4. update atom index

        pPairs = self.getPairs()
        self.changeAtypes(pPairs)
        self.delHydrogen(pPairs)
        self.updateIdx()
        self.top.avgCharge()
