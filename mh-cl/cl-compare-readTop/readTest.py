# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 14:05:13 2020

@author: huang
"""
import pandas as pd

def readFiles(name):
    atypeNames = ['name', 'bond_type', 'mass', 'charge', 'ptype', 'sigma', 'epsilon']
    btypeNames = ['ai', 'aj', 'funct', 'c0', 'c1']
    angTypeNames = ['ai', 'aj', 'ak', 'funct', 'c0', 'c1']
    dihTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
    impTypeNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
    moltypeNames = ['name', 'nrexcl']
    
    atNames = ['nr', 'type', 'resnr', 'residue', 'atom', 'cgnr', 'charge', 'mass']
    bNames = ['ai', 'aj', 'funct']
    pNames = ['ai', 'aj', 'funct']
    angNames = ['ai', 'aj', 'ak', 'funct']
    dihNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
    impNames = ['ai', 'aj', 'ak', 'al', 'funct', 'c0', 'c1', 'c2']
    
    names = [atypeNames, btypeNames, angTypeNames, dihTypeNames, impTypeNames, moltypeNames, 
             atNames, bNames, pNames, angNames, dihNames, impNames]
    
    with open(name, 'r') as f:
        lines = f.read().split('\n')
        
        lst0 = []; lst_tmp = []
        for l in lines:
            if l.startswith('['):
                lst0.append(lst_tmp)
                lst_tmp = []
            else:
                lst_tmp.append(l.split())
    
    lst0 = lst0[1:] # since the first line is start with a '[', it will create an empty ele

    df_lst = []
    for i in range(len(lst0)):
        print('list length: ', len(lst0[i][0]))
        print('names: ', names[i])
        print('names length: ', len(names[i]))
        df_tmp = pd.DataFrame(lst0[i], columns=names[i])
        df_lst.append(df_tmp)
    
    return lst0, df_lst
    
name = 'init.itp'

a1, a2 = readFiles(name)