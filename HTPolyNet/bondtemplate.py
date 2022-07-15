class BondTemplate:
    def __init__(self,names,resnames,order,bystander_resnames,oneaway_resnames):
        self.names=names
        self.resnames=resnames
        self.bystander_resnames=bystander_resnames
        self.oneaway_resnames=oneaway_resnames
        self.order=order
    def __str__(self):
        return f'BondTemplate {self.idx} resnames {self.resnames} order {self.order} bystander-resnames {self.bystander_resnames} oneaway-resnames {self.oneaway_resnames}'
    def __eq__(self,other):
        check=(self.names==other.names or self.names==other.names[::-1])
        check=check and (self.resnames==other.resnames or self.resnames==other.resnames[::-1])
        check=check and (self.bystander_resnames==other.bystander_resnames or self.bystander_resnames==other.bystander_resnames[::-1])
        check=check and (self.oneaway_resnames==other.oneaway_resnames or self.oneaway_resnames==other.oneaway_resnames[::-1])
        return check

BondTemplateList=list[BondTemplate]