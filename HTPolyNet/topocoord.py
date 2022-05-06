# class for methods that need to work with both Topology and Coordinates
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topology import Topology

class TopoCoord:
    def __init__(self):
        self.Coordinates=Coordinates()
        self.Topology=Topology()

    def make_bonds(self,pairs):
        idx_to_ignore=self.Coordinates.find_sacrificial_H(pairs,self.Topology)
        self.Topology.add_bonds(pairs,ignores=idx_to_ignore)
        idx_to_delete=self.Coordinates.find_sacrificial_H(pairs,self.Topology,rename=True)
        assert type(idx_to_delete)==list
        return idx_to_delete

    def delete_atoms(self,atomlist):
        self.Coordinates.delete_atoms(atomlist)
        idx_mapper=self.Topology.delete_atoms(atomlist)
        assert type(idx_mapper)==dict
        return idx_mapper
