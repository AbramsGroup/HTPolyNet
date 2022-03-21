# -*- coding: utf-8 -*-
"""

read cfg file
@author: huang, abrams

"""

# import pandas as pd
# import os

def getOrDie(D,k,basetype=None,subtype=str,source='config'):
    '''
    Strict assignment provider from dictionary 'D' at key 'k'.
    If key k is in dict D, it is cast to type:subtype and returned
    otherwise, an error is raised
    '''
    s=D.get(k,f'Error: Keyword {k} not found in {source}')
    if not 'Error:' in s:
        if basetype==list:
            # print(k,s,type(s))
            if type(s)==list:
                r=list(map(subtype,s))
            else:
                raise Exception(f'Value at keyword {k} is not type {basetype}.')
        elif basetype!=list: # ignore subtype
            if basetype==None:
                r=s
            else:
                r=basetype(s)
    else:
        raise Exception(s)
    return r

class Configuration(object):
    def __init__(self):
        self.cfgFile = ''
        self.cappingMolPair = []
        self.cappingBonds = []
        self.unrctStruct = []
        self.monInfo = ''
        self.croInfo = ''
        self.boxSize = ''
        self.monR_list = ''
        self.croR_list = ''
        self.cutoff = ''
        self.bondsRatio = 0.0
        self.maxBonds = 0
        self.desBonds = 0
        self.desConv = 0.0
        self.HTProcess = ''
        self.CPU = ''
        self.GPU = ''
        self.trials = ''
        self.reProject = ''
        self.rctInfo = ''
        self.stepwise = ''
        self.cappingBonds = []
        self.boxLimit = 1
        self.layerConvLimit = 1
        self.layerDir       = ''

    @classmethod
    def read(cls,filename):
        inst=cls()
        inst.cfgFile=filename
        baseList=[]
        # read all lines, append each to baseList.
        # skip any lines beginning with '#'
        # and ignore anything after '#' on a line
        with open(filename,'r') as f:
            for l in f:
                if l[0]!='#':
                    try:
                        l=l[:l.index('#')].strip()
                    except:
                        pass
                    baseList.append(l)
        # parse each item in baseList to build baseDict
        baseDict={}
        keyCounts={}
        for l in baseList:
            if '=' in l:  # this item is a parameter assignment
                k,v=[a.strip() for a in l.split('=')]
                if not k in keyCounts:
                    keyCounts[k]=0
                keyCounts[k]+=1
                if ',' in v:  # is it a csv?
                    v=[a.strip() for a in v.split(',')]
                elif ' ' in v:  # does it have spaces?
                    v=[a.strip() for a in v.split(',')]
                elif '\t' in v: # does it have tabs?
                    v=[a.strip() for a in v.split('\t')]
                # repeated keys build lists
                if keyCounts[k]>1:
                    if keyCounts[k]==2:
                        baseDict[k]=[baseDict[k]]
                    baseDict[k].append(v)
                else:    
                    baseDict[k]=v
                # print(k,baseDict[k])
        inst.baseDict=baseDict

        rctInfo = []
        for line in baseList: # Reaction Info
            if '+' in line:
                rct = [x.strip() for x in line.split('+')]
                rctInfo.append(rct)
        inst.rctInfo=rctInfo
        inst.parseCfg()
        return inst

    def __str__(self):
        r=f'Configuration read in from {self.cfgFile}:\n'
        r+='    reProject? '+(self.reProject if self.reProject!='' else '<no>')+'\n'
        r+='    boxLimit '+str(self.boxLimit)+'\n'
        r+='    layerConvLimit '+str(self.layerConvLimit)+'\n'
        if type(self.boxSize)==list:
            r+='    boxSize '+', '.join(str(x) for x in self.boxSize)+'\n'
        else:
            r+='    boxSize '+', '.join(str(x) for x in [self.boxSize]*3)+'\n'
        r+='    cutoff '+str(self.cutoff)+'\n'
        r+='    bondsRatio '+str(self.bondsRatio)+'\n'
        r+='    HTProcess '+self.HTProcess+'\n'
        r+='    CPU '+str(self.CPU)+'\n'
        r+='    GPU '+str(self.GPU)+'\n'
        r+='    trials '+str(self.trials)+'\n'
        r+='    stepwise '+(', '.join(self.stepwise) if self.stepwise!='' else 'UNSET')+'\n'
        r+='    layerDir '+(self.layerDir if self.layerDir!='' else 'UNSET')+'\n'
        r+='Monomer info:\n'
        r+='     idx     name    count    reactive-atoms\n'
        for m in self.monInfo:
            r+=f'     {m[0]:<8d}{m[1]:<8s}{m[2]:<8d}\n'
            r+=f'                              atname  z       rct   grp\n'
            for rr in m[3]:
                r+=f'                              {rr[0]:<8s}{rr[1]:<8d}{rr[2]:<6s}{rr[3]:<6s}\n'
        r+='Crosslinker info:\n'
        r+='     idx     name    count    reactive-atoms\n'
        for c in self.croInfo:
            r+=f'     {c[0]:<8d}{c[1]:<8s}{c[2]:<8d}\n'
            r+=f'                              atname  z       rct   grp\n'
            for rr in c[3]:
                r+=f'                              {rr[0]:<8s}{rr[1]:<8d}{rr[2]:<6s}{rr[3]:<6s}\n'
        r+='Summary: Mol names are '+', '.join(self.molNames)+'\n'
        r+='Calculated conversion info:\n'
        r+=f'    Desired conversion ("bondsRatio"): {self.desConv:<.3f}\n'
        r+=f'    Maximum number of new bonds:       {self.maxBonds}\n'
        r+=f'       ->  Target number of new bonds: {self.desBonds}\n'

        if len(self.cappingMolPair)>0:
            r+='Capping info:\n'
            r+='    Mol    unreacted\n'
            for p in self.cappingMolPair:
                r+=f'    {p[0]:<6s} {p[1]:<6s}\n'
            r+='    Mol    Atom1  Atom2  BondOrder\n'
            for p in self.cappingBonds:
                r+=f'    {p[0]:<6s} {p[1]:<6s} {p[2]:<6s} {p[3]:<6s}\n'
            r+='    UnreactedNames\n'
            for p in self.unrctStruct:
                r+=f'    {p:>6s}\n'
        r+='React info:\n'
        for x in self.rctInfo:
            r+=f'     {x}\n'
        # print(r)
        return r

    # def setName(self, filename):
    #     # this method will be superseded
    #     self.name = filename
    def calMaxBonds(self):
        maxRct = 0
        for i in self.monInfo:
            molNum = i[2]
            tmp = 0
            for ii in i[3]:
                tmp += ii[1]
            maxRct += molNum * tmp
            
        for i in self.croInfo:
            molNum = i[2]
            tmp = 0
            for ii in i[3]:
                tmp += ii[1]
            maxRct += molNum * tmp
        
#        print('type of maxRct',type(maxRct))
        self.maxBonds = int(maxRct * 0.5)
#        print('type of maxBonds',type(self.maxBonds))
        self.desConv = self.bondsRatio
#        print('type of desConv',type(self.desConv))
        self.desBonds = int(self.desConv * self.maxBonds)

    def getMolNames(self):
        names = []
        for n in self.unrctStruct:
            names.append(n)
        for n in self.monInfo:
            names.append(n[1])
        for n in self.croInfo:
            names.append(n[1])
        self.molNames = names

    def mol2sNeeded(self):
        # is this right, or do we also need mol2 for "unreacted/inactive"
        L=self.monInfo+self.croInfo
        for n in L:
            yield n[1]

    def parseCfg(self):
        self.cappingMolPair=[]
        self.unrctStruct=[]
        i=1
        while f'mol{i}' in self.baseDict:
            tokens=self.baseDict[f'mol{i}']
            self.cappingMolPair.append(tokens)
            self.unrctStruct.append(tokens[1])
            i+=1
        if 'cappingBonds' in self.baseDict:
            for cp in self.baseDict['cappingBonds']:
                tokens=cp
                self.cappingBonds.append(tokens)

        monInfo=[]
        monR_list={}
        i=1
        while f'monName{i}' in self.baseDict:
            mname=self.baseDict[f'monName{i}']
            mnum=getOrDie(self.baseDict,f'monNum{i}',basetype=int,source=self.cfgFile)
            mrnames=self.baseDict.get(f'mon{i}R_list',f'Error: mon{i}R_list not found')
            mrnum=self.baseDict.get(f'mon{i}R_rNum',f'Error: mon{i}R_rNum not found')
            mrrct=self.baseDict.get(f'mon{i}R_rct',f'Error: mon{i}R_rct not found')
            mrgrp=self.baseDict.get(f'mon{i}R_group',f'Error: mon{i}R_group not found')
            mrlist=[[r,int(n),x,g] for r,n,x,g in zip(mrnames,mrnum,mrrct,mrgrp)]
            monInfo.append([i,mname,mnum,mrlist])
            monR_list[mname] = mrlist
            i+=1
        croInfo=[]
        croR_list={}
        i=1
        while f'croName{i}' in self.baseDict:
            cname=self.baseDict[f'croName{i}']
            cnum=getOrDie(self.baseDict,f'croNum{i}',basetype=int,source=self.cfgFile)
            crnames=self.baseDict.get(f'cro{i}R_list',f'Error: cro{i}R_list not found')
            crnum=self.baseDict.get(f'cro{i}R_rNum',f'Error: cro{i}R_rNum not found')
            crrct=self.baseDict.get(f'cro{i}R_rct',f'Error: cro{i}R_rct not found')
            crgrp=self.baseDict.get(f'cro{i}R_group',f'Error: cro{i}R_group not found')
            crlist=[[r,int(n),x,g] for r,n,x,g in zip(crnames,crnum,crrct,crgrp)]
            croInfo.append([i,cname,cnum,crlist])
            croR_list[cname] = crlist
            i+=1
        self.monInfo = monInfo
        self.croInfo = croInfo
        self.monR_list = monR_list
        self.croR_list = croR_list

        # basic parameters that need not be specified
        # along with their default values
        self.reProject=self.baseDict.get('reProject','')
        self.boxLimit=float(self.baseDict.get('boxLimit','1'))
        self.layerConvLimit=float(self.baseDict.get('layerConvLimit','1'))
        self.GPU=int(self.baseDict.get('GPU','0'))
        self.stepwise=self.baseDict.get('stepwise','')
        self.layerDir=self.baseDict.get('boxDir','')

        # basic parameters that must be specified in the config file
        # boxSize can be specified by a single float for a cubic box
        # or as three floats for an orthorhombic box
        try:
            self.boxSize=getOrDie(self.baseDict,'boxSize',basetype=list,subtype=float,source=self.cfgFile)
        except:
            self.boxSize=getOrDie(self.baseDict,'boxSize',basetype=float,source=self.cfgFile)
        self.cutoff=getOrDie(self.baseDict,'cutoff',basetype=float,source=self.cfgFile)
        self.bondsRatio=getOrDie(self.baseDict,'bondsRatio',basetype=float,source=self.cfgFile)
        self.HTProcess=getOrDie(self.baseDict,'HTProcess',basetype=str,source=self.cfgFile)
        self.CPU=getOrDie(self.baseDict,'CPU',basetype=int,source=self.cfgFile)
        self.trials=getOrDie(self.baseDict,'trials',basetype=int,source=self.cfgFile)

        # anything that can be calculated immediately from data in the config
        self.calMaxBonds()
        self.getMolNames()

    # def readCfg(self):
    #     # this method will be superseded
    #     monInfo = []
    #     croInfo = []
        
    #     monNum = ''
    #     monR_list = {}
    #     croNum = ''
    #     croR_list = {}
        
    #     df = pd.read_csv(self.name, header=None, sep='\n', skip_blank_lines=True)
    #     df = df[df[0].str.startswith('#') == False]
    #     baseList = df.iloc[:][0]

    #     # Get capping parameters
    #     i = 1
    #     while i < 10:
    #         for l1 in baseList:
    #             if 'mol{}'.format(i) in l1:
    #                 tmpMolPair = l1.split('=')[1].split(',')
    #                 self.cappingMolPair.append(tmpMolPair)
    #                 self.unrctStruct.append(tmpMolPair[1].strip())
    #         i += 1

    #     for l1 in baseList:
    #         if l1.startswith('cappingBonds'):
    #             self.cappingBonds.append(l1.split('=')[1].split(','))

    #     # Get monomer and crosslinker info
    #     i = 1
    #     while i < 5:
    #         for l1 in baseList:
    #             if 'monName{}'.format(i) in l1:
    #                 monName = l1.split('=')[1].strip(' ')
    #                 for l2 in baseList:
    #                     key1 = 'monNum{}'.format(i)
    #                     if key1 in l2:
    #                         monNum = l2.split('=')[1].strip(' ')
                    
    #                 for l2 in baseList:
    #                     key2 = 'mon{}R_list'.format(i)
    #                     if key2 in l2:
    #                         monR_list_tmp = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
    #                 for l2 in baseList:
    #                     key3 = 'mon{}R_rNum'.format(i)
    #                     if key3 in l2:
    #                         monR_rNum = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
    #                 for l2 in baseList:
    #                     key4 = 'mon{}R_rct'.format(i)
    #                     if key4 in l2:
    #                         monR_rct = l2.split('=')[1].strip(' ').split('#')[0].split(',')

    #                 for l2 in baseList:
    #                     key4 = 'mon{}R_group'.format(i)
    #                     if key4 in l2:
    #                         monR_group = l2.split('=')[1].strip(' ').split('#')[0].split(',')

    #                 for idx in range(len(monR_list_tmp)):
    #                     monR_list_tmp[idx] = [monR_list_tmp[idx].strip(), monR_rNum[idx].strip(),
    #                                           monR_rct[idx].strip(), monR_group[idx].strip()]
                    
    #                 monInfo.append([i, monName, monNum, monR_list_tmp])
    #                 monR_list[monName] = monR_list_tmp
    #         i += 1

    #     i = 1
    #     while i < 5:
    #         for l1 in baseList:
    #             if 'croName{}'.format(i) in l1:
    #                 croName = l1.split('=')[1].strip(' ')
    #                 for l2 in baseList:
    #                     key1 = 'croNum{}'.format(i)
    #                     if key1 in l2:
    #                         croNum = l2.split('=')[1].strip(' ')
                    
    #                 for l2 in baseList:
    #                     key2 = 'cro{}R_list'.format(i)
    #                     if key2 in l2:
    #                         croR_list_tmp = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                                                
    #                 for l2 in baseList:
    #                     key3 = 'cro{}R_rNum'.format(i)
    #                     if key3 in l2:
    #                         croR_rNum = l2.split('=')[1].strip(' ').split('#')[0].split(',')
                    
    #                 for l2 in baseList:
    #                     key4 = 'cro{}R_rct'.format(i)
    #                     if key4 in l2:
    #                         croR_rct = l2.split('=')[1].strip(' ').split('#')[0].split(',')

    #                 for l2 in baseList:
    #                     key4 = 'cro{}R_group'.format(i)
    #                     if key4 in l2:
    #                         croR_group = l2.split('=')[1].strip(' ').split('#')[0].split(',')

    #                 for idx in range(len(croR_list_tmp)):
    #                     croR_list_tmp[idx] = [croR_list_tmp[idx].strip(), croR_rNum[idx].strip(),
    #                                           croR_rct[idx].strip(), croR_group[idx].strip()]
                                                
    #                 croInfo.append([i, croName, croNum, croR_list_tmp])
    #                 croR_list[croName] = croR_list_tmp
    #         i += 1

    #     reProject = '' # para could be missing in the options file

    #     for line in baseList: # Basic Info
    #         if 'boxSize' in line:
    #             boxSize = line.split('=')[1].strip(' ').split()    
    #         if 'cutoff' in line:
    #             cutoff = float(line.split('=')[1].strip(' '))
    #         if 'bondsRatio' in line:
    #             bondsRatio = line.split('=')[1].strip(' ')
    #         if 'HTProcess' in line:
    #             HTProcess = line.split('=')[1].strip(' ')
    #         if 'CPU' in line:
    #             CPU = line.split('=')[1].strip(' ')
    #         if 'GPU' in line:
    #             GPU = line.split('=')[1].strip(' ')
    #         if 'trials' in line:
    #             trials = line.split('=')[1].strip(' ')
    #         if 'reProject' in line:
    #             reProject = line.split('=')[1].strip(' ')
    #         if 'stepwise' in line:
    #             tmpStr =  line.split('=')[1]
    #             stepwise = tmpStr.split(',')
    #         if 'boxLimit' in line:
    #             boxLimit = line.split('=')[1].strip(' ')
    #         if 'boxDir' in line:
    #             boxDir = line.split('=')[1].strip(' ')
    #         if 'layerConvLimit' in line:
    #             layerConvLimit = line.split('=')[1].strip(' ')
    #     rctInfo = []
    #     for line in baseList: # React Info
    #         if '+' in line:
    #             rct = [x.split() for x in line.split('+')]
    #             rctInfo.append(rct)
            
    #     self.monInfo = monInfo
    #     self.croInfo = croInfo
    #     self.boxSize = boxSize
    #     self.monR_list = monR_list
    #     self.croR_list = croR_list
    #     self.cutoff = cutoff
    #     self.bondsRatio = bondsRatio
    #     self.rctInfo = rctInfo
    #     try:
    #         self.CPU = CPU
    #     except:
    #         self.CPU = 0
    #     try:
    #         self.GPU = GPU
    #     except:
    #         self.GPU = 0
    #     self.trials = trials
    #     self.HTProcess = HTProcess
    #     self.reProject = reProject
    #     self.stepwise = stepwise
    #     self.boxLimit = boxLimit
    #     self.layerDir = boxDir
    #     self.layerConvLimit = layerConvLimit
        
# if __name__ == '__main__':
#     pass
