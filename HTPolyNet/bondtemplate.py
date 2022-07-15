import logging
from copy import deepcopy

class BondTemplate:
    def __init__(self,names,resnames,order,bystander_resnames,oneaway_resnames):
        self.names=names
        self.resnames=resnames
        self.bystander_resnames=bystander_resnames
        self.oneaway_resnames=oneaway_resnames
        self.order=order
    def reverse(self):
        self.names=self.names[::-1]
        self.resnames=self.resnames[::-1]
        self.bystander_resnames=self.bystander_resnames[::-1]
        self.oneaway_resnames=self.oneaway_resnames[::-1]
    def __str__(self):
        return f'BondTemplate {self.names} resnames {self.resnames} order {self.order} bystander-resnames {self.bystander_resnames} oneaway-resnames {self.oneaway_resnames}'
    def __eq__(self,other):
        check=self.names==other.names
        check=check and self.resnames==other.resnames
        check=check and self.bystander_resnames==other.bystander_resnames
        check=check and self.oneaway_resnames==other.oneaway_resnames
        return check
    def is_reverse_of(self,other):
        rb=deepcopy(other)
        rb.reverse()
        return self==rb

BondTemplateList=list[BondTemplate]

class ReactionBond:
    def __init__(self,idx,resids,order,bystanders,oneaways):
        self.idx=idx
        self.resids=resids
        self.bystander_resids=bystanders
        self.oneaway_resids=oneaways
        self.order=order
    def reverse(self):
        self.idx=self.idx[::-1]
        self.resids=self.resids[::-1]
        self.bystander_resids=self.bystander_resids[::-1]
        self.oneaway_resids=self.oneaway_resids[::-1]
    def __str__(self):
        return f'ReactionBond {self.idx} resids {self.resids} order {self.order} bystander-resids {self.bystander_resids} oneaway-resids {self.oneaway_resids}'

ReactionBondList=list[ReactionBond]