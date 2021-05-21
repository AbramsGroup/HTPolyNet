import re
import sys

class endCapping(object):
    def __init__(self, inGro, inTop, inFFSum, inCappingBonds=[]):
        self.gro = inGro
        self.top = inTop
        self.cappingBonds = inCappingBonds
        self.topSum = inFFSum # [aTypes, bTypes, angTypes...], each is a dataframe
        self.VECapping()

    def changeAtypes(self, pPairs):
        df_atoms = self.top.atoms
        for p in pPairs:
            for a in p:
                atName = str(df_atoms.loc[(df_atoms.nr == a), 'atom'].to_list()[0])
                resName = str(df_atoms.loc[(df_atoms.nr == a), 'residue'].to_list()[0])
                tmpDf = self.topSum[-1]
                newAtomType = str(tmpDf.loc[(tmpDf.atom == atName) & (tmpDf.residue == resName), 'type'].to_list()[0])
                newAtomCharge = str(tmpDf.loc[(tmpDf.atom == atName) & (tmpDf.residue == resName), 'charge'].to_list()[0])

                print('newAtomType: ', newAtomType)
                print('newAtomCharge: ', newAtomCharge)

                df_atoms.loc[(df_atoms.nr == a), 'type'] = newAtomType
                df_atoms.loc[(df_atoms.nr == a), 'charge'] = newAtomCharge

    def getPairs(self):
        pPairs = []
        pAtoms = self.gro.df_atoms.loc[self.gro.df_atoms.rct == 'True']

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
                        if [a1GlobalIdx, a2GlobalIdx] in pPairs or [a2GlobalIdx, a1GlobalIdx] in pPairs:
                            continue
                        else:
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

        return con

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
            hCons1 = self.findHydrogen(p[0])
            hCons2 = self.findHydrogen(p[1])
            for a in hCons1:
                df_atoms.loc[(df_atoms.nr == a), 'type'] = 'ha'
            for a in hCons2:
                df_atoms.loc[(df_atoms.nr == a), 'type'] = 'ha'

            hCon1 = hCons1[0]
            hCon2 = hCons2[0]
            hAtoms.append(hCon1)
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

    def updateBasicType(self):
        self.top.atomtypes = self.top.atomtypes.append(self.topSum[0])
        self.top.bondtypes = self.top.bondtypes.append(self.topSum[1])
        self.top.angletypes = self.top.angletypes.append(self.topSum[2])
        self.top.dihtypes = self.top.dihtypes.append(self.topSum[3])
        self.top.imptypes = self.top.imptypes.append(self.topSum[4])
        # TODO: cannot remove duplicate like this
        self.top.atomtypes.drop_duplicates(inplace=True, ignore_index=True)
        self.top.bondtypes.drop_duplicates(inplace=True, ignore_index=True)
        self.top.angletypes.drop_duplicates(inplace=True, ignore_index=True)
        self.top.dihtypes.drop_duplicates(inplace=True, ignore_index=True)
        self.top.imptypes.drop_duplicates(inplace=True, ignore_index=True)

    def VECapping(self):
        # 1. Two atoms belong to the same molecule and same group --> unreact vinyl group
        # 2. change atom types
        # 3. remove h atoms
        # 4. update atom index

        pPairs = self.getPairs()
        print('Following atoms need to be capping: ')
        for p in pPairs:
            print('----> atoms {} and {}'.format(p[0], p[1]))
        self.updateBasicType()
        self.changeAtypes(pPairs)
        self.delHydrogen(pPairs)
        self.updateIdx()
        self.top.checkCharge()
