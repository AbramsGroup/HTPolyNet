import logging
from copy import deepcopy

class BondTemplate:
    def __init__(self,names,resnames,order,bystander_resnames,bystander_atomnames,oneaway_resnames,oneaway_atomnames):
        self.names=names
        self.resnames=resnames
        self.bystander_resnames=bystander_resnames
        self.bystander_atomnames=bystander_atomnames
        self.oneaway_resnames=oneaway_resnames
        self.oneaway_atomnames=oneaway_atomnames
        self.order=order
    def reverse(self):
        self.names=self.names[::-1]
        self.resnames=self.resnames[::-1]
        self.bystander_resnames=self.bystander_resnames[::-1]
        self.bystander_atomnames=self.bystander_atomnames[::-1]
        self.oneaway_resnames=self.oneaway_resnames[::-1]
        self.oneaway_atomnames=self.oneaway_atomnames[::-1]
    def __str__(self):
        return f'BondTemplate {self.names} resnames {self.resnames} order {self.order} bystander-resnames {self.bystander_resnames} bystander-atomnames {self.bystander_atomnames} oneaway-resnames {self.oneaway_resnames} oneaway-atomnames {self.oneaway_atomnames}'
    def __eq__(self,other):
        check=self.names==other.names
        check=check and self.resnames==other.resnames
        check=check and self.bystander_resnames==other.bystander_resnames
        check=check and self.bystander_atomnames==other.bystander_atomnames
        check=check and self.oneaway_resnames==other.oneaway_resnames
        check=check and self.oneaway_atomnames==other.oneaway_atomnames
        return check
    def is_reverse_of(self,other):
        rb=deepcopy(other)
        rb.reverse()
        return self==rb

BondTemplateList=list[BondTemplate]

class ReactionBond:
    def __init__(self,idx,resids,order,bystanders,bystanders_atomidx,oneaways,oneaways_atomidx):
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

class PassBond:
    def __init__(self,bondtuple,reactantname,distance,probability=1.0,order=1):
        assert len(bondtuple)==2
        assert type(bondtuple)==tuple
        self.bond=bondtuple
        self.reactantname=reactantname
        self.distance=distance
        self.probability=probability
        self.order=order
    def members(self):
        return self.bond[0],self.bond[1]