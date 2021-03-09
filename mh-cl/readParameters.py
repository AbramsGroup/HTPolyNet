# -*- coding: utf-8 -*-
"""

read parameters from basic/options.txt
@author: huang

"""

import pandas as pd

class parameters(object):
    def __init__(self):
        self.name = ''
        self.monInfo = ''
        self.croInfo = ''
        self.boxSize = ''
        self.monR_list = ''
        self.croR_list = ''
        self.cutoff = ''
        self.bondsRatio = ''
        self.maxBonds = ''
        self.HTProcess = ''
        self.rctInfo = ''
    
    def setName(self, name):
        self.name = name
        
    def readParam(self):
        monInfo = []
        croInfo = []
        
        monNum = ''
        monR_list = {}
        croNum = ''
        croR_list = {}
        df = pd.read_csv(self.name, header=None, sep='\n', skip_blank_lines=True)
        df = df[df[0].str.startswith('#') == False]
        baseList = df.iloc[:][0]
        
        # Get monomer and crosslinker info
        i = 1
        while i < 5:
            for l1 in baseList:
                if 'monName{}'.format(i) in l1:
                    monName = l1.split('=')[1].strip(' ')
                    for l2 in baseList:
                        key1 = 'monNum{}'.format(i)
                        if key1 in l2:
                            monNum = l2.split('=')[1].strip(' ')
                    
                    for l2 in baseList:
                        key2 = 'mon{}R_list'.format(i)
                        if key2 in l2:
                            monR_list_tmp = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
                    for l2 in baseList:
                        key3 = 'mon{}R_rNum'.format(i)
                        if key3 in l2:
                            monR_rNum = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
                    for l2 in baseList:
                        key4 = 'mon{}R_rct'.format(i)
                        if key4 in l2:
                            monR_rct = l2.split('=')[1].strip(' ').split('#')[0].split(',')

                    for l2 in baseList:
                        key4 = 'mon{}R_group'.format(i)
                        if key4 in l2:
                            monR_group = l2.split('=')[1].strip(' ').split('#')[0].split(',')

                    for idx in range(len(monR_list_tmp)):
                        monR_list_tmp[idx] = [monR_list_tmp[idx].strip(), monR_rNum[idx].strip(),
                                              monR_rct[idx].strip(), monR_group[idx].strip()]
                    
                    monInfo.append([i, monName, monNum, monR_list_tmp])
                    monR_list[monName] = monR_list_tmp
            i += 1
        
        i = 1
        while i < 5:
            for l1 in baseList:
                if 'croName{}'.format(i) in l1:
                    croName = l1.split('=')[1].strip(' ')
                    for l2 in baseList:
                        key1 = 'croNum{}'.format(i)
                        if key1 in l2:
                            croNum = l2.split('=')[1].strip(' ')
                    
                    for l2 in baseList:
                        key2 = 'cro{}R_list'.format(i)
                        if key2 in l2:
                            croR_list_tmp = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                                                
                    for l2 in baseList:
                        key3 = 'cro{}R_rNum'.format(i)
                        if key3 in l2:
                            croR_rNum = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
                    for l2 in baseList:
                        key4 = 'cro{}R_rct'.format(i)
                        if key4 in l2:
                            croR_rct = l2.split('=')[1].strip(' ').split('#')[0].split(',')

                    for l2 in baseList:
                        key4 = 'cro{}R_group'.format(i)
                        if key4 in l2:
                            croR_group = l2.split('=')[1].strip(' ').split('#')[0].split(',')

                    for idx in range(len(croR_list_tmp)):
                        croR_list_tmp[idx] = [croR_list_tmp[idx].strip(), croR_rNum[idx].strip(),
                                              croR_rct[idx].strip(), croR_group[idx].strip()]
                                                
                    croInfo.append([i, croName, croNum, croR_list_tmp])
                    croR_list[croName] = croR_list_tmp
            i += 1
                
        for line in baseList: # Basic Info
            if 'boxSize' in line:
                boxSize = line.split('=')[1].strip(' ')    
            if 'cutoff' in line:
                cutoff = float(line.split('=')[1].strip(' '))
            if 'bondsRatio' in line:
                bondsRatio = line.split('=')[1].strip(' ')
            if 'maxBonds' in line:
                maxBonds = line.split('=')[1].strip(' ')
            if 'HTProcess' in line:
                HTProcess = line.split('=')[1].strip(' ')

        rctInfo = []
        for line in baseList: # React Info
            if '+' in line:
                rct = [x.split() for x in line.split('+')]
                rctInfo.append(rct)
            
        self.monInfo = monInfo
        self.croInfo = croInfo
        self.boxSize = boxSize
        self.monR_list = monR_list
        self.croR_list = croR_list
        self.cutoff = cutoff
        self.bondsRatio = bondsRatio
        self.maxBonds = maxBonds
        self.rctInfo = rctInfo
        self.HTProcess = HTProcess