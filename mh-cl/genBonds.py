# -*- coding: utf-8 -*-
"""
Generate bonds
    - remove hydrogen atoms
    - add new connections
    - update atom index in each section
@author: huang
"""
import parameters
import sys

class genBonds(object):
    def __init__(self, gro, top, pairs, chargeMap, rctMols):
        self.gro = gro
        self.top = top
        self.pairs = pairs
        self.chargeMap = chargeMap # map
        self.rctMols = rctMols # list
        self.genPairs = []
        self.rctAtoms = []
        
        '''
        pairs data structure
        index   acro   amon   dist   rctNum1   rctNum2   molCon   p
        12      glbIdx glbIdx   
        '''
    
    def findHydrogen(self, atomsDf, bonds, idx): # This idx is the pd index
        con = []
        print('atom idx: ', idx)
        for index, row in bonds.iterrows():
            if str(row.ai) == str(idx):
                if 'H' in atomsDf[atomsDf.loc[:, 'globalIdx'] == str(row.aj)]['atomName'].values[0]:
                    con.append(row.aj)
            elif str(row.aj) == str(idx):
                if 'H' in atomsDf[atomsDf.loc[:, 'globalIdx'] == str(row.ai)]['atomName'].values[0]:
                    con.append(row.ai)
        print('con: ', con)
        return con
    
    def idx2Atypes(self, idx, df_atoms):
        aTypes = df_atoms.loc[int(idx)-1, 'type']
        return aTypes
    
    def checkNewTypes(self, pairs, sysTop, types='bonds'):
        sysTop = self.top
        df_atoms = sysTop.atoms # Top df
        if types == 'bonds':
            for p in pairs:
                bt = sysTop.bondtypes
                pairList = bt.loc[:, ['ai', 'aj']].values.tolist()
                a1Type = self.idx2Atypes(p[0], df_atoms)
                a2Type = self.idx2Atypes(p[1], df_atoms)
                if [a1Type, a2Type] in pairList or [a2Type, a1Type] in pairList:
                    continue
                else:
                    key = '{}-{}'.format(a1Type, a2Type)
                    print('{} pairs didn\'t show in the origin types. Searching the database...'.format(key))
                    param = parameters.dictBond[key]; lst_tmp = [a1Type, a2Type] + param
                    sysTop.addBondTypes(lst_tmp)
                    
        if types == 'angles':
            for p in pairs:
                angt = sysTop.angletypes
                pairList = angt.loc[:, ['ai', 'aj', 'ak']].values.tolist()
                a1Type = self.idx2Atypes(p[0], df_atoms)
                a2Type = self.idx2Atypes(p[1], df_atoms)
                a3Type = self.idx2Atypes(p[2], df_atoms)
                if [a1Type, a2Type, a3Type] in pairList or [a3Type, a2Type, a1Type] in pairList:
                    continue            
                else:
                    key = '{}-{}-{}'.format(a1Type, a2Type, a3Type)
                    print('{} pairs didn\'t show in the origin types. Searching the database...'.format(key))
                    param = parameters.dictAngle[key]; lst_tmp = [a1Type, a2Type, a3Type] + param
                    sysTop.addAngleTypes(lst_tmp)
        
        if types == 'dih':
            for p in pairs:
                diht = sysTop.dihtypes
                pairList = diht.loc[:, ['ai', 'aj', 'ak', 'al']].values.tolist()
                a1Type = self.idx2Atypes(p[0], df_atoms)
                a2Type = self.idx2Atypes(p[1], df_atoms)
                a3Type = self.idx2Atypes(p[2], df_atoms)
                a4Type = self.idx2Atypes(p[3], df_atoms)
                key1 =  '{}-{}-{}-{}'.format(a1Type, a2Type, a3Type, a4Type)
                key2 =  '{}-{}-{}-{}'.format(a4Type, a3Type, a2Type, a1Type)
                if [a1Type, a2Type, a3Type, a4Type] in pairList or [a4Type, a3Type, a2Type, a1Type] in pairList:
                    continue
                else:
                    if key1 in parameters.dictDihedral.keys():
                        param = parameters.dictDihedral[key1]
                        lst_tmp = [a1Type, a2Type, a3Type, a4Type] + param.split(',')
                        sysTop.addDihTypes(lst_tmp)
                    elif key2 in parameters.dictDihedral.keys():
                        param = parameters.dictDihedral[key2]
                        lst_tmp = [a4Type, a3Type, a2Type, a1Type] + param.split(',')
                        sysTop.addDihTypes(lst_tmp)
                    else:
                        sys.exit('Unknow dihedral type{}, need to find param for the pair'.format(key1))
                        
    def delHydrogen(self):
        '''
        1. Generate new bonds/pairs/angles/dihs
        2. Check new types if exist in the database
        3. Find and delete hydrogen
        '''
        
        pairs = []
        for index, row in self.pairs.iterrows():
            pairs.append([row.acro, row.amon])
#        pairs = self.pairs
        self.genPairs = pairs
        atomsDf = self.gro.df_atoms
        inTop = self.top
        df_atoms = inTop.atoms # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers
        hAtoms = []
        for p in pairs:
            print('pairs: ', p)
            hCon1 = self.findHydrogen(atomsDf, df_bonds, p[0])[0]
            hCon2 = self.findHydrogen(atomsDf, df_bonds, p[1])[0]
            hAtoms.append(hCon1); hAtoms.append(hCon2)
        
        print('Following atoms will be removed: ', hAtoms)
        for a in hAtoms:
#            gro_df = sysGro[0]
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
        
    def searchCon(self, idx, df_bonds):
        con = []
        for index, row in df_bonds.iterrows():
            if str(row.ai) == str(idx):
                con.append(row.aj)
            elif str(row.aj) == str(idx):
                con.append(row.ai)
        return con
    
    def genNewCon(self, pair, df_bonds):
        new_bonds = []; new_pairs = []; new_angles = []; new_dihedrals = []
        a1 = str(pair[0]); a2 = str(pair[1])
        new_bonds.append([a1, a2])
        con1 = self.searchCon(a1, df_bonds); con2 = self.searchCon(a2, df_bonds)
        for a in con1:
            a = str(a)
            if a != 1:#row.amon:
                #row.acro = a1; row.amon = a2
                new_angles.append([a, a1, a2])
                con3 = self.searchCon(a, df_bonds); con4 = self.searchCon(a2, df_bonds)
                for aa in con3:
                    aa = str(aa)
                    if aa != a1: # sth wrong with this
                        new_dihedrals.append([aa, a, a1, a2])
                        new_pairs.append([aa, a2])
                for aa in con4:
                    aa = str(aa)
                    if aa != a1:
                        new_dihedrals.append([a, a1, a2, aa])
                        new_pairs.append([a, aa])
        for a in con2:
            a = str(a)
            if a != 1:#row.acro:
                new_angles.append([a1, a2, a])
                con3 = self.searchCon(a1, df_bonds); con4 = self.searchCon(a, df_bonds)
                for aa in con3:
                    aa = str(aa)
                    if aa != a2:
                        if [aa, a1, a2, a] not in new_dihedrals:
                            new_dihedrals.append([aa, a1, a2, a])
                            new_pairs.append([aa, a])
                for aa in con4:
                    aa = str(aa)
                    if aa != a2:
                        if [a1, a2, a, aa] not in new_dihedrals:
                            new_dihedrals.append([a1, a2, a, aa])
                            new_pairs.append([a1, aa])
                            
        return new_bonds, new_pairs, new_angles, new_dihedrals
    
    def addNewCon(self):
        inTop = self.top
        df_atoms = inTop.atoms # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers
        new_bonds = []; new_pairs = []; new_angles = []; new_dihedrals = []
        
        pairs = self.genPairs
        for p in pairs:
            nBonds, nPairs, nAngles, nDihs = self.genNewCon(p, df_bonds)
            new_bonds += nBonds
            new_pairs += nPairs
            new_angles += nAngles
            new_dihedrals += nDihs
            
        #TODO: update charge
        # check and add new types to the corresponding type section
        print('checking and adding new types...')
        self.checkNewTypes(new_bonds, inTop, types='bonds')
        self.checkNewTypes(new_angles, inTop, types='angles')
        self.checkNewTypes(new_dihedrals, inTop, types='dih')

        inTop.atoms = df_atoms
        inTop.bonds = df_bonds
        inTop.pairs = df_pairs
        inTop.angles = df_angs
        inTop.dihedrals = df_dihs
        inTop.impropers = df_imps
        
        # add new pairs to the corresponding section
        print('generate new connection bonds/pairs/angles/dihs/imps...')
        inTop.addBonds(new_bonds)
        inTop.addPairs(new_pairs)
        inTop.addAngles(new_angles)
        inTop.addDih(new_dihedrals)
#        inTop.addImp(new_dihedrals) # Currently don't find improper, add it back when it needed
        
    def updateTopIdx(self, x, gro_df, types='atoms'):
        if types == 'atoms':
            idx1 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.nr)]['globalIdx'].values[0]
            x.nr = str(idx1); x.new_idx = str(idx1)
        elif types == 'bonds':
            idx1 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ai)]['globalIdx'].values[0]
            idx2 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.aj)]['globalIdx'].values[0]
            x.ai = str(idx1); x.aj = str(idx2)
        elif types == 'pairs':
            idx1 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ai)]['globalIdx'].values[0]
            idx2 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.aj)]['globalIdx'].values[0]
            x.ai = str(idx1); x.aj = str(idx2)
        elif types == 'angles':
            idx1 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ai)]['globalIdx'].values[0]
            idx2 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.aj)]['globalIdx'].values[0]
            idx3 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ak)]['globalIdx'].values[0]
            x.ai = str(idx1); x.aj = str(idx2); x.ak = str(idx3)
        elif types == 'dih':
            idx1 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ai)]['globalIdx'].values[0]
            idx2 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.aj)]['globalIdx'].values[0]
            idx3 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.ak)]['globalIdx'].values[0]
            idx4 = gro_df[gro_df.loc[:, 'ori_idx'] == str(x.al)]['globalIdx'].values[0]
            x.ai = str(idx1); x.aj = str(idx2); x.ak = str(idx3); x.al = str(idx4)
        else:
            print('sth wrong')
        return x
    
    def updateIdx(self):
        atomsDf = self.gro.df_atoms
        inTop = self.top
        df_atoms = inTop.atoms # Top df
        df_bonds = inTop.bonds
        df_pairs = inTop.pairs
        df_angs = inTop.angles
        df_dihs = inTop.dihedrals
        df_imps = inTop.impropers
        
        atomsDf['ori_idx'] = atomsDf['globalIdx']
        
        atomsDf = atomsDf.reset_index(drop=True)
        atomsDf['globalIdx'] = (atomsDf.index + 1).astype(str).to_list()
        df_atoms_new = df_atoms.apply(lambda x: self.updateTopIdx(x, atomsDf, types='atoms'), axis=1).reset_index(drop=True)
        df_bonds_new = df_bonds.apply(lambda x: self.updateTopIdx(x, atomsDf, types='bonds'), axis=1)
        df_pairs_new = df_pairs.apply(lambda x: self.updateTopIdx(x, atomsDf, types='pairs'), axis=1)
        df_angs_new = df_angs.apply(lambda x: self.updateTopIdx(x, atomsDf, types='angles'), axis=1)
        df_dihs_new = df_dihs.apply(lambda x: self.updateTopIdx(x, atomsDf, types='dih'), axis=1)
        df_imps_new = df_imps.apply(lambda x: self.updateTopIdx(x, atomsDf, types='dih'), axis=1)
        inTop.atoms = df_atoms_new
        inTop.bonds = df_bonds_new
        inTop.pairs = df_pairs_new
        inTop.angles = df_angs_new
        inTop.dihedrals = df_dihs_new
        inTop.impropers = df_imps_new
        self.gro.df_atoms = atomsDf; self.gro.atNum = len(atomsDf)
    
    def getSeq(self, df):
        seq = list(df.type)
        name = list(df.residue)[0]
        return seq, name
    
    def mapCharge(self, resname, seq): 
        print('seq: ', seq)
        charges = self.chargeMap
        for keys, value in charges.items():
            if resname in keys and seq == value[0]:
                print('find map!: ', keys)
                return value[1]
    
    def getMolIdx(self, molNum):
        df_atoms = self.gro.df_atoms
        rctIdx = list(df_atoms[df_atoms.molNum == molNum].globalIdx)
        df_atoms_top = self.top.atoms[(self.top.atoms.nr.isin(rctIdx))]
        seq, resname = self.getSeq(df_atoms_top)
        charges = self.mapCharge(resname, seq)
        
        self.top.atoms.loc[(self.top.atoms.nr.isin(rctIdx)), 'charge'] = charges
        
    def updateCharge(self): # Didn't test this function
        mols = self.rctMols
        m = mols[1]
        a = self.getMolIdx(m)
        return a
    
    def main(self):
        self.delHydrogen()
        self.addNewCon()
        self.updateIdx()
        
if __name__ == "__main__":
    import readGro
    import readTop
    import groInfo
    import topInfo
    atomsDf = groInfo.gro()
    topDf = topInfo.top()
    
    a2 = readGro.initGro()
    a3 = readTop.initTop()
    
    a2.setName('nvt-1')
    df_init, sysName, atNum, boxSize = a2.readGRO()
    atomsDf.setGroInfo(df_init, sysName, atNum, boxSize)
    a3.setName('init.top', 'init.itp')
    a3.genTopSession()
    topDf.setInfo(a3.sumTop)
    topDf.checkCharge()
    rctMols = ['11', '2', '6', '10']
    chargeMaps = b
    pairs = [['505', '145']]
    b0 = genBonds(atomsDf, topDf, pairs, chargeMaps, rctMols)
    aa0 = b0.updateCharge()
    aa1 = b0.top.atoms
    
#    b.main()
    
    