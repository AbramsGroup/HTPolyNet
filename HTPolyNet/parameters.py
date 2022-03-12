# extra GAFF parameters handling

from HTPolyNet.libraries import *
import json
import os
class ExtraGAFFParams:
    def __init__(self):
        lp=IdentifyLibraryResourcePaths()
        extra_gaff=lp['gaff_extra']
        with open(os.path.join(extra_gaff,'bonds.json'),'r') as f:
            self.extra_bonds=json.load(f)
        with open(os.path.join(extra_gaff,'angles.json'),'r') as f:
            self.extra_angles=json.load(f)
        with open(os.path.join(extra_gaff,'dihedrals.json'),'r') as f:
            self.extra_dihedrals=json.load(f)
