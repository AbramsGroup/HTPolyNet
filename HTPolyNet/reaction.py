import logging
from copy import deepcopy
from enum import Enum
import networkx as nx
from HTPolyNet.plot import network_graph
logger=logging.getLogger(__name__)

class reaction_stage(Enum):
    """Enumerated reaction stage
    """
    build=0   # only used for building molecules that use parameterized oligomers
    param=1   # used to generate parameterized oligomers that are not cure reactions
    cure=2    # used to generate parameterized oligormers that ARE cure reactions
    cap=3     # used to generate parameterized capped monomers
    unset=99
    def __str__(self):
        return self.name

class Reaction:
    default_directives={
        'name':'',
        'atoms':{},
        'bonds':[],
        'reactants':{},
        'product':'',
        'stage':reaction_stage.unset,
        'procession':{},
        'probability':1.0
    }
    def __init__(self,jsondict={}):
        self.jsondict=jsondict
        self.name=jsondict.get('name',self.default_directives['name'])
        self.atoms=jsondict.get('atoms',self.default_directives['atoms'])
        self.bonds=jsondict.get('bonds',self.default_directives['bonds'])
        self.reactants=jsondict.get('reactants',self.default_directives['reactants'])
        self.product=jsondict.get('product',self.default_directives['product'])
        self.stage=reaction_stage[jsondict.get('stage',str(self.default_directives['stage']))]
        self.procession=jsondict.get('procession',self.default_directives['procession'])
        self.probability=jsondict.get('probability',self.default_directives['probability'])
        for label in jsondict:
            if not label in self.default_directives:
                logging.debug(f'Ignoring unknown reaction directive "{label}"')

        self.symmetry_versions=[]
    
    def __str__(self):
        retstr=f'Reaction "{self.name}" ({str(self.stage)})\n'
        for i,r in self.reactants.items():
            retstr+=f'reactant {i}: {r}\n'
        retstr+=f'product {self.product}\n'
        for i,a in self.atoms.items():
            retstr+=f'atom {i}: {a}\n'
        for b in self.bonds:
            retstr+=f'bond {b}\n'
        return retstr

ReactionList = list[Reaction]

def parse_reaction_list(baselist:ReactionList):
    ''' baselist is read directly in from cfg, just a list of unprocessed dicts
        generates any new reactions that are inferred by procession (polymerization)
    '''
    processives=[]
    for i,R in enumerate(baselist):
        if R.procession:
            nReactions=[]
            logger.debug(f'Reaction {R.name} is processive: {R.procession}')
            final_product_name=R.product
            R.product=f'{final_product_name}_I0'
            original_name=R.name
            R.name=f'{R.name}-i0'
            base_reactant_key=R.procession['increment_resid']
            for c in range(R.procession['count']):
                nR=deepcopy(R)
                nR.name=f'{original_name}-i{c+1}'
                nR.product=f'{final_product_name}_I{c+1}' if c<(R.procession['count']-1) else final_product_name
                nR.reactants[base_reactant_key]=f'{final_product_name}_I{c}'
                for ak,arec in nR.atoms.items():
                    if arec['reactant']==base_reactant_key: arec['resid']+=(c+1)
                nR.procession={}
                nReactions.append(nR)
            processives.append((i,nReactions))
    lasti=0
    ret_reactions=baselist.copy()
    for rec in processives:
        atidx,sublist=rec
        atidx+=lasti
        ret_reactions=ret_reactions[0:atidx+1]+sublist+ret_reactions[atidx+1:]
        lasti+=len(sublist)
    return ret_reactions

def extract_molecule_reactions(rlist:ReactionList):
    working_rlist=rlist.copy()
    G=nx.DiGraph()
    # Given an unsorted list of reactions (list of reactants +  one product), order a list of 
    # of molecules so that no product comes before any of its reactants
    if not rlist: return []
    molecule_react_order=[]
    reactants=set()
    products=set()
    # first molecules to generate are the input reactants
    for R in working_rlist:
        for r in R.reactants.values():
            reactants.add(r)
            G.add_edge(r,R.product)
        products.add(R.product)
    network_graph(G,'plots/reaction_network.png',arrows=True)
    input_reactants=reactants.intersection(reactants.symmetric_difference(products))
    logger.debug(f'Input reactants: {input_reactants}')
    for i in input_reactants:
        molecule_react_order.append((i,None))
    # once inputs are generated, any moleculed constructed only from inputs can be generated
    for R in rlist:
        if all([x in input_reactants for x in R.reactants.values()]):
            molecule_react_order.append((R.product,R))
            working_rlist.remove(R)
    
    while len(working_rlist)>0:
        molecules=[x[0] for x in molecule_react_order]
        for R in working_rlist:
            if all([x in molecules for x in R.reactants.values()]):
                molecule_react_order.append((R.product,R))
                molecules.append(R.product)
                working_rlist.remove(R)
            else:
                logger.debug(f'Cannot place {R.product} since not all of {R.reactants.values()} are in the list yet')

    return molecule_react_order

def get_r(mname,RL:ReactionList):
    if not RL: return None
    i=0
    while RL[i].product!=mname: 
        i+=1
        if i==len(RL):
            return None
    return RL[i]

def is_reactant(name:str,reaction_list:ReactionList,stage=reaction_stage.cure):
    reactants=[]
    for r in reaction_list:
        if r.stage==stage:
            for v in r.reactants.values():
                if not v in reactants:
                    reactants.append(v)
    return name in reactants
    
def product_sequence_resnames(R:Reaction,reactions:ReactionList):
    result=[]
    allProducts=[x.product for x in reactions]
    for rKey,rName in R.reactants.items():
        if rName in allProducts:
            result.extend(product_sequence_resnames(reactions[allProducts.index(rName)],reactions))
        else:
            result.extend([rName])
    return result

def molname_sequence_resnames(molname:str,reactions:ReactionList):
    prodnames=[x.product for x in reactions]
    if molname in prodnames:
        R=reactions[prodnames.index(molname)]
        return product_sequence_resnames(R,reactions)
    else:
        return [molname]  # not a product in the reaction list?  must be a monomer

def reactant_resid_to_presid(R:Reaction,reactantName:str,resid:int,reactions:ReactionList):
    reactants=[x for x in R.reactants.values()]
    if reactantName in reactants:
        ridx=reactants.index(reactantName)
        logger.debug(f'reactantName {reactantName} {ridx} in {reactants} of {R.name}')
        S=0
        for i in range(ridx):
            seq=molname_sequence_resnames(reactants[i],reactions)
            logger.debug(f'{i} {seq} {S}')
            S+=len(seq)
        return S+resid
    else:
        return -1

def generate_product_name(R:Reaction):
    pname=''
    for b in R.bonds:
        i,j=b['atoms']
        i_arec=R.atoms[i]
        j_arec=R.atoms[j]
        i_reKey=i_arec['reactant']
        i_reNm=R.reactants[i_reKey]
        i_aNm=i_arec['atom']
        j_reKey=j_arec['reactant']
        j_reNm=R.reactants[j_reKey]
        j_aNm=j_arec['atom']
        tbNm=f'{i_reNm}~{i_aNm}-{j_aNm}~{j_reNm}'
        if len(R.reactants)>1:
            if len(pname)==0:
                pname=tbNm
            else:
                pname+='---'+tbNm
    return pname

# def get_atom_options(R:Reaction,symmetry_relateds:dict): #,reactions:ReactionList):
#     # prod_seq_resn=product_sequence_resnames(R,reactions)
#     atom_options=[]
#     for atomKey,atomRec in R.atoms.items():
#         reactantKey=atomRec['reactant']
#         reactantName=R.reactants[reactantKey]
#         resid=atomRec['resid']
#         # resName=prod_seq_resn[resid-1]
#         atomName=atomRec['atom']
#         # logger.debug(f'resName {resName} atomName {atomName}')
#         symm_sets=symmetry_relateds.get(reactantName,[])
#         # logger.debug(f'symm_set {symm_sets}')
#         if symm_sets:
#             for symm_set in symm_sets:
#                 # logger.debug(f'symm_set {symm_set}')
#                 if atomName in symm_set:
#                     atom_options.append([[atomKey,c] for c in symm_set])
#     return atom_options

