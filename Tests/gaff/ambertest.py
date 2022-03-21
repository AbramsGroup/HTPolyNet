import HTPolyNet.ambertools as ha
import HTPolyNet.libraries as li
import os

lp=li.IdentifyLibraryResourcePaths(['mol2'])
os.system(f'cp {os.path.join(lp["mol2"],"mST.mol2")} .')
p=ha.Parameterization()
p.GAFFParameterize('mST.mol2','testmsT')
os.system(f'cp {os.path.join(lp["mol2"],"DGEBA.mol2")} .')
p.GAFFParameterize('DGEBA.mol2','testDBEGA',extra_antechamber_params='-eq 1 -pl 10',parmed_save_inline=False)


