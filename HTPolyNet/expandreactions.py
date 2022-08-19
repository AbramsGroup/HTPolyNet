from itertools import product
from copy import deepcopy
import logging
from HTPolyNet.molecule import Molecule, MoleculeList, MoleculeDict
from HTPolyNet.reaction import reaction_stage, Reaction, ReactionList, get_atom_options, generate_product_name, reactant_resid_to_presid

logger=logging.getLogger(__name__)

def symmetry_expand_reactions(reactions:ReactionList,symmetry_relateds:dict,MD:MoleculeDict):
    extra_reactions=[]
    extra_molecules={}
    logger.debug(f'begins: symmetry_relateds: {symmetry_relateds}')
    jdx=1
    for R in reactions:
        thisR_extra_reactions=[]
        thisR_extra_molecules={}
        atom_options=[]
        logger.debug(f'Symmetry expansion of {R.name} begins')
        # logger.debug(f'  Product {R.product} resname sequence {prod_seq_resn}')
        atom_options=get_atom_options(R,symmetry_relateds)#,reactions)
        logger.debug(f'  atom options: {atom_options}')
        if len(R.reactants)>1:
            olist=list(product(*atom_options))
        else:
            olist=list(zip(*atom_options))
        idx=1
        R.symmetry_versions=olist
        for P in olist[1:]:
            newR=deepcopy(R)
            newR.name=R.name+f'-S{idx}'
            logger.debug(f'Permutation {P}:')
            for pp in P:
                atomKey,atomName=pp
                newR.atoms[atomKey]['atom']=atomName
            pname=generate_product_name(newR)
            if len(pname)==0:
                pname=R.product+f'-{idx}'
            newR.product=pname 
            newR.stage=R.stage
            logger.debug(f'Primary:')
            for ln in str(newR).split('\n'):
                logger.debug(ln)
            thisR_extra_reactions.append(newR)
            thisR_extra_molecules[newR.product]=Molecule(name=newR.product,generator=newR)
            thisR_extra_molecules[newR.product].set_origin('symmetry_product')
            thisR_extra_molecules[newR.product].set_sequence_from_moldict(MD)
            for rR in [x for x in reactions if R.product in x.reactants.values()]:
                reactantKey=list(rR.reactants.keys())[list(rR.reactants.values()).index(R.product)]
                logger.debug(f'  product {newR.product} must replace reactantKey {reactantKey} in {rR.name}')
                nooR=deepcopy(rR)
                nooR.stage=rR.stage
                nooR.name=rR.name+f'-{reactantKey}:S{jdx}'
                nooR.reactants[reactantKey]=newR.product
                # update any atom names to reflect origin of this reactant
                for naK,naRec in {k:v for k,v in nooR.atoms.items() if v['reactant']==reactantKey}.items():
                    na_resid=naRec['resid'] # resid of reactant atom in target reactant
                    na_name=naRec['atom']
                    for p in P:
                        oaK,oa_name=p 
                        oaRec=R.atoms[oaK]
                        oa_reactatnName=R.reactants[oaRec['reactant']]
                        oa_resid=oaRec['resid']
                        oa_resid_in_o_product=reactant_resid_to_presid(R,oa_reactatnName,oa_resid,reactions)
                        # this atom is an atom in the permutation the resid in product matches
                        if na_resid == oa_resid_in_o_product:
                            nooR.atoms[naK]['resid']=oa_resid_in_o_product
                            nooR.atoms[naK]['atom']=oa_name
                noor_pname=generate_product_name(nooR)
                if len(noor_pname)==0:
                    noor_pname=rR.product+f'-{jdx}'
                nooR.product=noor_pname
                logger.debug(f'Secondary:')
                for ln in str(nooR).split('\n'):
                    logger.debug(ln)
                jdx+=1
                reactions.append(nooR)
                thisR_extra_molecules[nooR.product]=Molecule.New(nooR.product,nooR)
                thisR_extra_molecules[nooR.product].set_origin('symmetry_product')
                thisR_extra_molecules[newR.product].set_sequence_from_moldict(MD)
            idx+=1
        logger.debug(f'Symmetry expansion of reaction {R.name} ends')

        # done with this reaction
        extra_reactions.extend(thisR_extra_reactions)
        extra_molecules.update(thisR_extra_molecules)
    # done with all reactions

    return extra_reactions,extra_molecules

def chain_expand_reactions(molecules:MoleculeDict):
    ''' must be called after all grx attributes are set for all molecules '''
    extra_reactions:ReactionList=[]
    extra_molecules:MoleculeDict={}
    monomers:MoleculeList=[]
    dimer_lefts:MoleculeList=[]
    dimer_rights:MoleculeList=[]
    for mname,M in molecules.items():
        if len(M.sequence)==1 and len(M.TopoCoord.idx_lists['chain'])>0 and M.generator==None and M.parentname==M.name:
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
                R.stage=reaction_stage.param
                R.name=new_mname.lower()
                R.product=new_mname
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.set_origin('unparameterized')
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
                R.stage=reaction_stage.param
                new_rxnname=new_mname.lower()
                R.name=new_rxnname
                R.product=new_mname
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.set_origin('unparameterized')
                extra_molecules[R.product]=newP
                logger.debug(f'dimer atom {d.name}[{d.sequence[1]}2_{h_name}] will attack monomer atom {m.name}_{t_name}-> {new_mname}:')
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
                R.stage=reaction_stage.param
                R.product=new_mname
                R.name=R.product.lower()
                # newP=Molecule(name=R.product,generator=R)
                newP=Molecule.New(R.product,R).set_sequence_from_moldict(molecules)
                newP.set_origin('unparameterized')
                extra_molecules[R.product]=newP
                logger.debug(f'dimer atom {dr.name}-{dr.sequence[1]}2_{h_name} will attack dimer atom {dl.name}-{dl.sequence[0]}1_{t_name} -> {new_mname}:')
                for ln in str(R).split('\n'):
                    logger.debug(ln)
                extra_reactions.append(R)

    return extra_reactions, extra_molecules
