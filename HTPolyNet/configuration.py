"""

read cfg file
@author: huang, abrams

"""
import json
import yaml
import os
import logging
import numpy as np
from copy import deepcopy
from itertools import product

from HTPolyNet.molecule import Molecule, Reaction

logger=logging.getLogger(__name__)

def _determine_sequence(m,moldict,atoms):
    if not moldict[m].generator:
        return [m],[atoms]
    thisseq=[]
    thisatm=[]
    for rid,mname in moldict[m].generator.reactants.items():
        atoms=[a for a in moldict[m].generator.atoms.values() if a['reactant']==rid]
        newseq,newatm=_determine_sequence(mname,moldict,atoms)
        thisatm.extend(newatm)
        thisseq.extend(newseq)
    return thisseq,thisatm

class Configuration:
    def __init__(self):
        self.cfgFile = ''
        self.Title = ''
        ''' List of (Molecule, count) '''
        self.initial_composition = []
        ''' Dictionary of name:Molecule '''
        self.molecules = {}
        ''' List of Reaction instances '''
        self.reactions = []
        ''' all other parameters in cfg file '''
        self.parameters = {}

    @classmethod
    def read(cls,filename):
        extension=filename.split('.')[-1]
        if extension=='json':
            return cls.read_json(filename)
        elif extension=='yaml' or extension=='yml':
            return cls.read_yaml(filename)
        else:
            raise Exception(f'Unknown config file extension {extension}')

    @classmethod
    def read_json(cls,filename):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=json.load(f)
        inst.parse()
        return inst

    @classmethod
    def read_yaml(cls,filename):
        inst=cls()
        inst.cfgFile=filename
        with open(filename,'r') as f:
            inst.basedict=yaml.safe_load(f)
        inst.parse()
        return inst
    
    def NewMolecule(self,mol_name):
        M=Molecule(mol_name)
        if 'symmetry_equivalent_atoms' in self.parameters:
            if mol_name in self.parameters['symmetry_equivalent_atoms']:
                M.symmetry_relateds=self.parameters['symmetry_equivalent_atoms'][mol_name]
        if 'stereocenters' in self.parameters:
            if mol_name in self.parameters['stereocenters']:
                M.stereocenters=self.parameters['stereocenters'][mol_name]
                extra_stereocenters=[]
                for stc in M.stereocenters:
                    for sc in M.symmetry_relateds:
                        if stc in sc:
                            sc_copy=sc.copy()
                            sc_copy.remove(stc)
                            extra_stereocenters.extend(sc_copy)
                M.stereocenters.extend(extra_stereocenters)
        return M

    def parse(self):
        """parse parse self.basedict to set Title, initial_composition, and lists of
           reactions and molecules.
           
        """
        self.Title=self.basedict.get('Title','No Title Provided')
        self.parameters=self.basedict
        if not 'ncpu' in self.parameters:
            self.parameters['ncpu']=os.cpu_count()
        '''
        add a molecule for each unique molecule specified in initial_composition
        '''
        self.initial_composition=self.basedict.get('initial_composition',[])
        for item in self.initial_composition:
            m=item['molecule']
            if m not in self.molecules:
                logger.debug(f'new {m}')
                self.molecules[m]=self.NewMolecule(m)
        '''
        add additional molecules that appear as either reactants or products of reactions
        '''
        rlist=self.basedict.get('reactions',[])
        self.reactions=[Reaction(r) for r in rlist]
        for R in self.reactions:
            '''
            add every reactant in this reaction to the list of molecules
            '''
            for rnum,rname in R.reactants.items():
                zrecs=[]
                logger.debug(f'{R.name} rname {rname}')
                for atnum,atrec in R.atoms.items():
                    if atrec['reactant']==rnum:
                        cprec=atrec.copy()
                        del cprec['reactant']
                        # this atom is in this reactant
                        zrecs.append(cprec)
                if not rname in self.molecules:
                    self.molecules[rname]=self.NewMolecule(rname)
                ''' provide molecule with records of atoms that have z values '''
                self.molecules[rname].update_zrecs(zrecs)
        for R in self.reactions:
            '''
            add product of this reaction or update generator if product is already in the list of molecules
            '''
            r=R.product
            if not r in self.molecules:
                self.molecules[r]=self.NewMolecule(r)
            self.molecules[r].generator=R

        ''' 
        
        '''
        # for mname,M in self.molecules.items():
        #     logger.debug(f'{M.name} {M.zrecs}')
        #     # M.sequence,M.reactive_atoms_seq=_determine_sequence(mname,self.molecules,[])
            # logger.debug(f'Sequence of {mname}: {M.sequence}')
            # logger.debug(f'Reactive atoms per sequence: {M.reactive_atoms_seq}')

    def symmetry_expand_reactions(self,molecules):
        extra_reactions=[]
        extra_molecules={}
        for R in self.reactions:
            atom_options=[]
            for atomKey,atomRec in R.atoms.items():
                reactantKey=atomRec['reactant']
                reactantName=R.reactants[reactantKey]
                resid=atomRec['resid']
                atomName=atomRec['atom']
                rmol=molecules[reactantName]
                asidx=rmol.TopoCoord.get_gro_attribute_by_attributes('sea-idx',{'resid':resid,'atomName':atomName})
                logger.debug(f'asidx {asidx}')
                if asidx>-1:
                    clu=rmol.atoms_w_same_attribute_as(find_dict={'atomName':atomName,'resid':resid},same_attribute='sea-idx',
                    return_attribute='atomName')
                    atom_options.append([[atomKey,c] for c in clu])
            logger.debug(f'{R.name} {list(R.reactants.values())} {R.product} {atom_options}')
            if len(R.reactants)>1:
                olist=list(product(*atom_options))
            else:
                olist=list(zip(*atom_options))
            idx=1
            for P in olist[1:]:
                logger.debug(f'{P}')
                newR=deepcopy(R)
                newR.name=R.name+f'-S{idx}'
                idx+=1
                for pp in P:
                    atomKey,atomName=pp
                    newR.atoms[atomKey]['atom']=atomName
                pname=''
                for b in newR.bonds:
                    i,j=b['atoms']
                    i_arec=newR.atoms[i]
                    j_arec=newR.atoms[j]
                    i_reKey=i_arec['reactant']
                    i_reNm=newR.reactants[i_reKey]
                    i_aNm=i_arec['atom']
                    j_reKey=j_arec['reactant']
                    j_reNm=newR.reactants[j_reKey]
                    j_aNm=j_arec['atom']
                    tbNm=f'{i_reNm}~{i_aNm}-{j_aNm}~{j_reNm}'
                    if len(R.reactants)>1:
                        if len(pname)==0:
                            pname=tbNm
                        else:
                            pname+='---'+tbNm
                    else:
                        pname=R.product+f'-{idx}'
                newR.product=pname
                newR.stage=R.stage
                logger.debug(f'new reaction {newR.name} product {newR.product}')
                extra_reactions.append(newR)
                newP=Molecule(name=newR.product,generator=newR)
                extra_molecules[newR.product]=newP
        self.reactions.extend(extra_reactions)
        return extra_molecules

    # def x_symmetry_expand_reactions(self,unique_molecules):
    #     # logger.debug('symmetry_expand_reactions')
    #     extra_reactions=[]
    #     current_molecules=unique_molecules
    #     extra_molecules={}
    #     trydict={}
    #     sclass={}
    #     for R in self.reactions:
    #         logger.debug(f'Symmetry-expanding reaction {R.name}')
    #         pro_seq,pro_ra=_determine_sequence(R.product,current_molecules,[])
    #         logger.debug(f'{R.product} {pro_seq} {pro_ra}')
    #         sseq=[]
    #         for mname,matrectlist in zip(pro_seq,pro_ra):
    #             # how many se versions of this reactant?
    #             # -> product of number of members of symmetry class of each reactive atom in this
    #             # reactant
    #             residue=current_molecules[mname]
    #             if not mname in sclass:
    #                 sclass[mname]={}
    #             sp=[]
    #             for atrec in matrectlist:
    #                 atomName=atrec['atom']
    #                 seaidx=residue.TopoCoord.get_gro_attribute_by_attributes('sea-idx',{'atomName':atomName})
    #                 if seaidx<=-1:
    #                     logger.debug(f'No atoms symmetric with {mname} {atomName}')
    #                     continue
    #                 clu=residue.atoms_w_same_attribute_as(find_dict={'atomName':atomName},    
    #                                             same_attribute='sea-idx',
    #                                             return_attribute='atomName')
    #                 # make sure atomName is the first element in clu
    #                 clu=list(clu)
    #                 assert atomName in clu
    #                 clu.remove(atomName)
    #                 clu.insert(0,atomName)
    #                 logger.debug(f'Atoms symmetry-equivalent to {mname} {atomName}: {clu}')
    #                 sclass[mname][atomName]=clu
    #                 logger.debug(f'sclass {sclass}')
    #                 sp.append(clu)

    #             # logger.debug(f'sp {sp}')
    #             if len(R.reactants)>1:
    #                 sseq.append(list(product(*sp)))
    #             else:
    #                 sseq.append(list(zip(*sp)))
    #         logger.debug(f'sclass {sclass}')
    #         logger.debug(f'Reaction reactant sealist: {sseq}')                      
    #         P=list(product(*[v for v in sseq]))
    #         logger.debug(f'P {P}')
    #         idx=1
    #         if len(P)==0:
    #             continue
    #         trydict[P[0]]=R.product
    #         for p in P[1:]: # skip the first permutation -- assume this is the explicit one!
    #             logger.debug(f'{R.name}-{idx} {p}({len(p)}) -> {R.product}-{idx}')
    #             trydict[p]=f'{R.product}-{idx}'
    #             newR=deepcopy(R)
    #             newR.name+=f'-{idx}'
    #             newR.product+=f'-{idx}'
    #             logger.debug(f'copy {R.name} with reactants {R.reactants} to {newR.name} with {newR.reactants}')
    #             ip=0
    #             pdiv={}
    #             if len(R.reactants)==1: # unimolecular reaction
    #                 pdiv[1]=p
    #             else:
    #                 # map patterns to previously created resnames
    #                 for nri,nR in R.reactants.items():
    #                     logger.debug(f'reactant {nri}: {nR}')
    #                     nres_nR=len(current_molecules[nR].sequence)
    #                     logger.debug(f'nres_nR {nres_nR}')
    #                     subp=p[ip:ip+nres_nR]
    #                     ip+=nres_nR
    #                     nresname=nR
    #                     if subp in trydict:
    #                         nresname=trydict[subp]
    #                     logger.debug(f'  {nR} {nresname} {subp}')
    #                     newR.reactants[nri]=nresname
    #                     pdiv[nri]=subp
    #             for aL,atomrec in R.atoms.items():
    #                 atomName=atomrec['atom']
    #                 z=atomrec['z']
    #                 resid=atomrec['resid']
    #                 reactant_idx=atomrec['reactant']
    #                 reactant=R.reactants[reactant_idx]
    #                 resname=current_molecules[reactant].sequence[resid-1]
    #                 logger.debug(f'need to alter atom: {atomName} resid: {resid} reactant: {reactant_idx} ({reactant})')
    #                 logger.debug(f'have {pdiv[reactant_idx]} {len(pdiv[reactant_idx])}')
    #                 new_atomName=atomName
    #                 for q in pdiv[reactant_idx]:
    #                     logger.debug(f'{q}')
    #                     for qq in q:
    #                         logger.debug(f'is it {qq}?')
    #                         isit=qq in sclass[resname][atomName]
    #                         logger.debug(f'{isit}')
    #                         if isit:
    #                             new_atomName=qq
    #                 newR.atoms[aL]={'atom':new_atomName,'resid':resid,'reactant':reactant_idx,'z':z}
    #                 # pdiv[] here has one element for each reactive atom in reactant, but
    #             extra_reactions.append(newR)
    #             newP=Molecule(name=newR.product,generator=newR)
    #             logger.debug(f'...generated {newR.name} as generator of {newP.name}')
    #             logger.debug(f'   {newR.name} uses reactants {newR.reactants}')
    #             extra_molecules[newR.product]=newP
    #             idx+=1
    #     # logger.debug(f'trydict {trydict}')
    #     # for R in self.reactions:
    #     #     logger.debug(R)
    #     # logger.debug(f'extra reactions:')
    #     # for R in extra_reactions:
    #     #     logger.debug(R)
    #     self.reactions.extend(extra_reactions)
    #     return extra_molecules

    def chain_expand_reactions(self,molecules):
        extra_reactions=[]
        extra_molecules={}
        monomers=[]
        dimer_lefts=[]
        dimer_rights=[]
        for mname,M in molecules.items():
            if len(M.sequence)==1 and len(M.TopoCoord.idx_lists['chain'])>0 and M.generator==None:
                monomers.append(M)
            elif len(M.sequence)==2:
                A=molecules[M.sequence[0]]
                if len(A.TopoCoord.idx_lists['chain'])>0:
                    dimer_lefts.append(M)
                A=molecules[M.sequence[1]]
                if len(A.TopoCoord.idx_lists['chain'])>0:
                    dimer_rights.append(M)
        for mon in monomers:
            cnms=[]
            for c in mon.TopoCoord.idx_lists['chain']:
                cnms.append([mon.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in c])
            logger.debug(f'Monomer {mon.name} has {len(mon.TopoCoord.idx_lists["chain"])} 2-chains: {mon.TopoCoord.idx_lists["chain"]} {cnms}')

        for dim in dimer_lefts:
            logger.debug(f'Dimer_left {dim.name} has sequence {dim.sequence}')
            logger.debug(f'-> chains: {dim.TopoCoord.idx_lists["chain"]}')
            for cl in dim.TopoCoord.idx_lists['chain']:
                nl=[dim.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl]
                logger.debug(f'  -> {nl}')
        for dim in dimer_rights:
            logger.debug(f'Dimer_right {dim.name} has sequence {dim.sequence}')
            logger.debug(f'-> chains: {dim.TopoCoord.idx_lists["chain"]}')
            for cl in dim.TopoCoord.idx_lists['chain']:
                nl=[dim.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':x}) for x in cl]
                logger.debug(f'  -> {nl}')

        # monomer head attacks dimer tail
        MD=product(monomers,dimer_lefts)
        for m,d in MD:
            for mb in m.TopoCoord.idx_lists['chain']:
                h_idx=mb[0]
                h_name=m.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
                # by definition, the dimer must have one chain of length 4
                D4=[]
                for dc in d.TopoCoord.idx_lists['chain']:
                    if len(dc)==4:
                        D4.append(dc)
                for DC in D4:
                    t_idx=DC[-1]
                    t_name=d.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
                    new_mname=f'{m.name}~{h_name}={t_name}~{d.name}'
                    '''construct reaction'''
                    R=Reaction()
                    R.reactants={1:m.name, 2:d.name}
                    R.atoms={'A':{'reactant':1,'resid':1,'atom':h_name,'z':1},
                            'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                    R.bonds=[{'atoms':['A','B'],'order':1}]
                    R.stage='template-only'
                    R.name=new_mname.lower()
                    R.product=new_mname
                    newP=Molecule(name=R.product,generator=R)
                    extra_molecules[R.product]=newP
                    logger.debug(f'monomer atom {m.name}_{h_name} will attack dimer atom {d.name}[{d.sequence[0]}1_{t_name}] -> {new_mname}:')
                    for ln in str(R).split('\n'):
                        logger.debug(ln)
                    extra_reactions.append(R)
        # dimer head attacks monomer tail
        MD=product(monomers,dimer_rights)
        for m,d in MD:
            for mb in m.TopoCoord.idx_lists['chain']:
                assert len(mb)==2,f'monomer {m.name} has a chain that is not length-2 -- this is IMPOSSIBLE'
                t_idx=mb[-1]
                t_name=m.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
                D4=[]
                for dc in d.TopoCoord.idx_lists['chain']:
                    if len(dc)==4:
                        D4.append(dc)
                for DC in D4:
                    h_idx=DC[0]
                    h_name=d.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
                    new_mname=f'{d.name}~{h_name}={t_name}~{m.name}'
                    '''construct reaction'''
                    R=Reaction()
                    R.reactants={1:d.name, 2:m.name}
                    R.atoms={'A':{'reactant':1,'resid':2,'atom':h_name,'z':1},
                            'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                    R.bonds=[{'atoms':['A','B'],'order':1}]
                    R.stage='template-only'
                    new_rxnname=new_mname.lower()
                    R.name=new_rxnname
                    R.product=new_mname
                    newP=Molecule(name=R.product,generator=R)
                    extra_molecules[R.product]=newP
                    logger.debug(f'dimer atom {d.name}[{d.sequence[1]}2_{h_name}] will attach monomer atom {m.name}_{t_name}-> {new_mname}:')
                    for ln in str(R).split('\n'):
                        logger.debug(ln)
                    extra_reactions.append(R)

        DD=product(dimer_rights,dimer_lefts)
        for dr,dl in DD:
            ''' head of dr attacks tail of dl '''
            for cr,cl in product(dr.TopoCoord.idx_lists['chain'],dl.TopoCoord.idx_lists['chain']):
                if len(cr)==4 and len(cl)==4:
                    h_idx=cr[0]
                    h_name=dr.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':h_idx})
                    t_idx=cl[-1]
                    t_name=dl.TopoCoord.get_gro_attribute_by_attributes('atomName',{'globalIdx':t_idx})
                    new_mname=f'{dr.name}~{h_name}={t_name}~{dl.name}'
                    '''construct reaction'''
                    R=Reaction()
                    R.reactants={1:dr.name, 2:dl.name}
                    R.atoms={'A':{'reactant':1,'resid':2,'atom':h_name,'z':1},
                             'B':{'reactant':2,'resid':1,'atom':t_name,'z':1}}
                    R.bonds=[{'atoms':['A','B'],'order':1}]
                    R.stage='template-only'
                    R.product=new_mname
                    R.name=R.product.lower()
                    newP=Molecule(name=R.product,generator=R)
                    extra_molecules[R.product]=newP
                    logger.debug(f'dimer atom {dr.name}-{dr.sequence[1]}2_{h_name} will attack dimer atom {dl.name}-{dl.sequence[0]}1_{t_name} -> {new_mname}:')
                    for ln in str(R).split('\n'):
                        logger.debug(ln)
                    extra_reactions.append(R)

        self.reactions.extend(extra_reactions)
        return extra_molecules

    def get_reaction(self,product_name):
        if '-' in product_name:
            base_product_name,sym=product_name.split('-')
        else:
            base_product_name,sym=product_name,'0'
        for r in self.reactions:
            if r.product==base_product_name:
                return r
        return None

    def calculate_maximum_conversion(self):
        N={}
        for item in self.initial_composition:
            N[item['molecule']]=item['count']
        Bonds=[]
        Atoms=[]
        for R in self.reactions:
            if R.stage=='cure':
                for b in R.bonds:
                    A,B=b['atoms']
                    a,b=R.atoms[A],R.atoms[B]
                    aan,ban=a['atom'],b['atom']
                    ari,bri=a['resid'],b['resid']
                    arnum,brnum=a['reactant'],b['reactant']
                    arn,brn=R.reactants[arnum],R.reactants[brnum]
                    if arnum!=brnum:
                        az,bz=a['z'],b['z']
                        ia=(aan,ari,arnum,arn,az)
                        ib=(ban,bri,brnum,brn,bz)
                        b=(ia,ib)
                        if ia not in Atoms and arn in N:
                            Atoms.append(ia)
                        if ib not in Atoms and brn in N:
                            Atoms.append(ib)
                        if b not in Bonds and arn in N and brn in N:
                            Bonds.append(b)
        # logger.debug(f'atomset: {Atoms}')
        Z=[]
        for a in Atoms:
            Z.append(a[4]*N[a[3]])
        # logger.debug(f'Z: {Z}')
        # logger.debug(f'bondset: {Bonds}')
        MaxB=[]
        for B in Bonds:
            a,b=B
            az=Z[Atoms.index(a)]
            bz=Z[Atoms.index(b)]
            MaxB.append(min(az,bz))
            Z[Atoms.index(a)]-=MaxB[-1]
            Z[Atoms.index(b)]-=MaxB[-1]
        # logger.debug(f'MaxB: {MaxB} {sum(MaxB)}')
        return sum(MaxB)
