# -*- coding: utf-8 -*-
"""

read cfg file
@author: huang, abrams

"""

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
                raise KeyError(f'Value at keyword {k} is not type {basetype}.')
        elif basetype!=list: # ignore subtype
            if basetype==None:
                r=s
            else:
                r=basetype(s)
    else:
        raise KeyError(s)
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
        self.system='A Generic System Name'

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
        
        self.maxBonds = int(maxRct * 0.5)
        self.desConv = self.bondsRatio
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
        self.system=self.baseDict.get('system','A Generic System Name')
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
