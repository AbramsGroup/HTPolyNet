"""

.. module:: reaction
   :synopsis: handles Reaction objects
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
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
        """__init__ generates a Reaction object using directives in jsondict

        :param jsondict: directives for creating a new Reaction object, defaults to {}
        :type jsondict: dict, optional
        """
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
    """parse_reaction_list generates any new reactions inferred by procession (i.e., linear polymerization)

    :param baselist: list of unprocessed reaction dicts read directly from configuration file
    :type baselist: ReactionList
    :return: entire set of all reactions, including new ones added here
    :rtype: ReactionList
    """
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

def extract_molecule_reactions(rlist:ReactionList,plot=True):
    """extract_molecule_reactions establishes the correct order of reactions so that molecular templates are created in a sensible order (reactants before products for all reactions)

    :param rlist: list of all reactions
    :type rlist: ReactionList
    :param plot: flag to plot the reaction network, defaults to True
    :type plot: bool, optional
    :return: ordered list of tuples, reach of the from (name-of-product-molecule,Reaction)
    :rtype: list of tuples
    """
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
    if plot: network_graph(G,'plots/reaction_network.png',arrows=True,with_labels=True)
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
            # else:
            #     logger.debug(f'Cannot place {R.product} since not all of {R.reactants.values()} are in the list yet')

    return molecule_react_order

def get_r(mname,RL:ReactionList):
    """get_r returns the Reaction that generates molecule mname as a product, if it exists, otherwise None

    :param mname: name of molecule to search 
    :type mname: str
    :param RL: list of all Reactions
    :type RL: ReactionList
    :return: Reaction that generates mname, or None
    :rtype: Reaction or NoneType
    """
    if not RL: return None
    i=0
    while RL[i].product!=mname: 
        i+=1
        if i==len(RL):
            return None
    return RL[i]

def is_reactant(name:str,reaction_list:ReactionList,stage=reaction_stage.cure):
    """is_reactant tests to see if molecule 'name' appears as a reactant in any reaction in reaction_list

    :param name: name of molecule
    :type name: str
    :param reaction_list: list of all Reactions
    :type reaction_list: ReactionList
    :param stage: stage or reactions to search, defaults to reaction_stage.cure
    :type stage: reaction_stage(Enum), optional
    :return: True if molecule is a reactant, False otherwise
    :rtype: bool
    """
    reactants=[]
    for r in reaction_list:
        if r.stage==stage:
            for v in r.reactants.values():
                if not v in reactants:
                    reactants.append(v)
    return name in reactants
    
def product_sequence_resnames(R:Reaction,reactions:ReactionList):
    """product_sequence_resnames recursively generate sequence of product of R by traversing network of reactions; sequence is by definition an ordered list of monomer names

    :param R: a Reaction
    :type R: Reaction
    :param reactions: list of all Reactions
    :type reactions: ReactionList
    :return: the sequence of R's product molecule as a list of monomer names
    :rtype: list
    """
    result=[]
    allProducts=[x.product for x in reactions]
    for rKey,rName in R.reactants.items():
        if rName in allProducts:
            result.extend(product_sequence_resnames(reactions[allProducts.index(rName)],reactions))
        else:
            result.extend([rName])
    return result

def molname_sequence_resnames(molname:str,reactions:ReactionList):
    """molname_sequence_resnames determine monomer sequence of the molecule with name molname based on a traversal of the reaction network

    :param molname: name of molecule
    :type molname: str
    :param reactions: list of all Reactions
    :type reactions: ReactionList
    :return: monomer sequence of molecule as a list of monomer names
    :rtype: list
    """
    prodnames=[x.product for x in reactions]
    if molname in prodnames:
        R=reactions[prodnames.index(molname)]
        return product_sequence_resnames(R,reactions)
    else:
        return [molname]  # not a product in the reaction list?  must be a monomer

def reactant_resid_to_presid(R:Reaction,reactantName:str,resid:int,reactions:ReactionList):
    """reactant_resid_to_presid map the resid of a reactant monomer to its inferred resid in the associate product in reaction R, by traversing the reaction network

    :param R: a Reaction
    :type R: Reaction
    :param reactantName: the name of one of the reactants in R
    :type reactantName: str
    :param resid: the resid in that reactant we want mapped to the product
    :type resid: int
    :param reactions: list of all Reactions
    :type reactions: ReactionList
    :return: resid of this residue in the product of R
    :rtype: int
    """
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
    """generate_product_name automatically generate the name of the product based on the names of the reactants and the bonds in R

    :param R: a Reaction
    :type R: Reaction
    :return: suggested name of product
    :rtype: str
    """
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

