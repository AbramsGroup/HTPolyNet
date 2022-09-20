"""

.. module:: bondtemplate
   :synopsis: Manages bond templates (bonds defined by type) and reaction bonds (bonds defined by instances)
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
from copy import deepcopy

class BondTemplate:
    def __init__(self,names,resnames,intraresidue,order,bystander_resnames,bystander_atomnames,oneaway_resnames,oneaway_atomnames):
        """__init__ create a BondTemplate object

        :param names: names of atoms in the bond
        :type names: list-like container of two ints
        :param resnames: names of the two residues to which the atoms belong
        :type resnames: list-like container of two strs
        :param intraresidue: True if this is an intraresidue bond
        :type intraresidue: bool
        :param order: bond order (1=single, 2=double, ...)
        :type order: int
        :param bystander_resnames: lists of names of bystander residues (residues also bound to one of the atoms)
        :type bystander_resnames: list of two list-like containers of strs, one for each atom in names
        :param bystander_atomnames: lists of names of atoms in bystander residues
        :type bystander_atomnames: list of two list-like containers of strs (parallel to bystander_resnames)
        :param oneaway_resnames: names of one-away residues (residues bound one bond away from the new interresidue bond; only relevant for C=C free-radical polymerization)
        :type oneaway_resnames: list-like container of strs
        :param oneaway_atomnames: names of atoms in one-away residues
        :type oneaway_atomnames: list-like container of strs
        """
        self.names=names
        self.resnames=resnames
        self.intraresidue=intraresidue
        self.bystander_resnames=bystander_resnames
        self.bystander_atomnames=bystander_atomnames
        self.oneaway_resnames=oneaway_resnames
        self.oneaway_atomnames=oneaway_atomnames
        self.order=order
    def reverse(self):
        """reverse the order of all parallel lists in a BondTemplate object
        """
        self.names=self.names[::-1]
        self.resnames=self.resnames[::-1]
        self.bystander_resnames=self.bystander_resnames[::-1]
        self.bystander_atomnames=self.bystander_atomnames[::-1]
        self.oneaway_resnames=self.oneaway_resnames[::-1]
        self.oneaway_atomnames=self.oneaway_atomnames[::-1]
    def __str__(self):
        return f'BondTemplate {self.names} resnames {self.resnames} intraresidue? {self.intraresidue} order {self.order} bystander-resnames {self.bystander_resnames} bystander-atomnames {self.bystander_atomnames} oneaway-resnames {self.oneaway_resnames} oneaway-atomnames {self.oneaway_atomnames}'
    def __eq__(self,other):
        check=self.names==other.names
        check=check and self.intraresidue==other.intraresidue
        check=check and self.resnames==other.resnames
        check=check and self.bystander_resnames==other.bystander_resnames
        check=check and self.bystander_atomnames==other.bystander_atomnames
        check=check and self.oneaway_resnames==other.oneaway_resnames
        check=check and self.oneaway_atomnames==other.oneaway_atomnames
        return check
    def is_reverse_of(self,other):
        """is_reverse_of return True if self and other are reverse of each other

        :param other: another BondTemplate object
        :type other: BondTemplate
        :return: True if self and other are reverse copies of each other
        :rtype: bool
        """
        rb=deepcopy(other)
        rb.reverse()
        return self==rb

BondTemplateList=list[BondTemplate]

class ReactionBond:
    def __init__(self,idx,resids,order,bystanders,bystanders_atomidx,oneaways,oneaways_atomidx):
        """__init__ generate a new ReactionBond object

        :param idx: indices of two atoms that form the bond
        :type idx: list-like container of two ints
        :param resids: resids of the two atoms that form the bond
        :type resids: list-like container of two ints
        :param order: order of the bond (1=single, 2=double, ...)
        :type order: int
        :param bystanders: two lists of indices of interresidue resids already bound to each atom in idx 
        :type bystanders: list of two list-like containers of ints
        :param bystanders_atomidx: two lists of indices of interresidue atoms already bound to each atom in idx 
        :type bystanders_atomidx: list of two list-like containers of ints
        :param oneaways: list of one-away resids
        :type oneaways: list of ints
        :param oneaways_atomidx: list of one-away atom indices
        :type oneaways_atomidx: list of ints
        """
        self.idx=idx
        self.resids=resids
        self.bystander_resids=bystanders
        self.bystander_atomidx=bystanders_atomidx
        self.oneaway_resids=oneaways
        self.oneaway_atomidx=oneaways_atomidx
        self.order=order
    def reverse(self):
        self.idx=self.idx[::-1]
        self.resids=self.resids[::-1]
        self.bystander_resids=self.bystander_resids[::-1]
        self.bystander_atomidx=self.bystander_atomidx[::-1]
        self.oneaway_resids=self.oneaway_resids[::-1]
        self.oneaway_atomidx=self.oneaway_atomidx[::-1]
    def __str__(self):
        return f'ReactionBond {self.idx} resids {self.resids} order {self.order} bystander-resids {self.bystander_resids} oneaway-resids {self.oneaway_resids} oneaway-atomidx {self.oneaway_atomidx}'

ReactionBondList=list[ReactionBond]
