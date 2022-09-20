"""

.. module:: expandreactions
   :synopsis: handles expansion of reactions based on chains
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
from itertools import product
import logging
from HTPolyNet.molecule import Molecule, MoleculeList, MoleculeDict
from HTPolyNet.reaction import reaction_stage, Reaction, ReactionList

logger=logging.getLogger(__name__)

def chain_expand_reactions(molecules:MoleculeDict):
    """chain_expand_reactions handles generation of new reactions and molecular templates implied by any C-C chaining

    :param molecules: dictionary of molecular templates constructed from explicit declarations
    :type molecules: MoleculeDict
    :return: the list of new reactions and dictionary of new molecules
    :rtype: tuple(ReactionList,MoleculeDict)
    """
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
