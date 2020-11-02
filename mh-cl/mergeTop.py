# -*- coding: utf-8 -*-
"""
Created on Tue Sep 29 14:45:45 2020

@author: huangming
"""

import topInfo
import pandas as pd
import time

def addIdx(row, atIdx, keys):
    if keys == 'atoms':
        row.nr = str(int(row.nr) + int(atIdx))
        row.cgnr = str(int(row.cgnr) + int(atIdx))
        
    elif keys == 'bonds':
        row.ai = str(int(row.ai) + int(atIdx))
        row.aj = str(int(row.aj) + int(atIdx))
    elif keys == 'pairs':
        row.ai = str(int(row.ai) + int(atIdx))
        row.aj = str(int(row.aj) + int(atIdx))
    elif keys == 'angles':
        row.ai = str(int(row.ai) + int(atIdx))
        row.aj = str(int(row.aj) + int(atIdx))
        row.ak = str(int(row.ak) + int(atIdx))
    elif keys == 'dih':
        row.ai = str(int(row.ai) + int(atIdx))
        row.aj = str(int(row.aj) + int(atIdx))
        row.ak = str(int(row.ak) + int(atIdx))
        row.al = str(int(row.al) + int(atIdx))
    return row

def updateIdx(df, atIdx, keys):
    df = df.apply(lambda x: addIdx(x, atIdx, keys), axis=1)
    return df

def mergeInfo(atomslist, bondslist, pairslist, anglelist, dihlist, implist, topSum):
    count = 0
    atomlst_new = []; bondlst_new = []; pairlst_new = []; anglelist_new = []; dihlist_new = []; implist_new = []
    t1 = time.time()
    for i in range(len(atomslist)):
        df_atoms = updateIdx(atomslist[i], count, 'atoms')
        df_bonds = updateIdx(bondslist[i], count, 'bonds')
        df_pairs = updateIdx(pairslist[i], count, 'pairs')
        df_angles = updateIdx(anglelist[i], count, 'angles')
        df_dih = updateIdx(dihlist[i], count, 'dih')
        df_imp = updateIdx(implist[i], count, 'dih')
        atomlst_new.append(df_atoms); bondlst_new.append(df_bonds); pairlst_new.append(df_pairs)
        anglelist_new.append(df_angles); dihlist_new.append(df_dih); implist_new.append(df_imp)
        count = int(df_atoms.nr.iloc[-1])
    t2 = time.time()
    print('tt_merge: ', t2 - t1)
    topSum.atoms = pd.concat(atomlst_new).reset_index(drop=True)
    topSum.bonds = pd.concat(bondlst_new).reset_index(drop=True)
    topSum.pairs = pd.concat(pairlst_new).reset_index(drop=True)
    topSum.angles = pd.concat(anglelist_new).reset_index(drop=True)
    topSum.dihedrals = pd.concat(dihlist_new).reset_index(drop=True)
    topSum.impropers = pd.concat(implist_new).reset_index(drop=True)
    
    return topSum
    
def mergeTypes(typelist):
    df = pd.concat(typelist).drop_duplicates().reset_index(drop=True)
    return df

def mergeTopList(topList):
    topInit = topList[0]
    topSum = topInfo.top()
    topSum.system = topInit.system
    topSum.molecules = topInit.molecules
    topSum.moleculetype = topInit.moleculetype
#    topSum.dupDihTypeKey = topInit.dupDihTypeKey
    dup = []
    for i in range(len(topList)):
        dup += topList[i].dupDihTypeKey
    dup = list(set(dup))
    topSum.dupDihTypeKey = dup
    atomtypes = []; bondtypes = []; angletypes = []; dihtypes = []; imptypes = []
    atomslist = {}; bondslist = {}; pairslist = {}; anglelist = {}; dihlist = {}; implist = {}
    
    for i in range(len(topList)):
        atomtypes.append(topList[i].atomtypes)
        bondtypes.append(topList[i].bondtypes)
        angletypes.append(topList[i].angletypes)
        dihtypes.append(topList[i].dihtypes)
        imptypes.append(topList[i].imptypes)
    
        topList[i].molNum = i + 1
        atomslist[i] = topList[i].atoms.copy()
        bondslist[i] = topList[i].bonds.copy()
        pairslist[i] = topList[i].pairs.copy()
        anglelist[i] = topList[i].angles.copy()
        dihlist[i] = topList[i].dihedrals.copy()
        implist[i] = topList[i].impropers.copy()
   
    topSum.atomtypes = mergeTypes(atomtypes)
    topSum.bondtypes = mergeTypes(bondtypes)
    topSum.angletypes = mergeTypes(angletypes)
    topSum.dihtypes = mergeTypes(dihtypes)
    topSum.imptypes = mergeTypes(imptypes)
    topSum = mergeInfo(atomslist, bondslist, pairslist, anglelist, dihlist, implist, topSum)
    
    return topSum