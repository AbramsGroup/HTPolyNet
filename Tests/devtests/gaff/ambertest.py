import HTPolyNet.ambertools as ha
import HTPolyNet.libraries as li
import parmed as pm
import os

lp=li.IdentifyLibraryResourcePaths(['mol2'])
os.system(f'cp {os.path.join(lp["mol2"],"mST.mol2")} .')
p=ha.Parameterization()
p.GAFFParameterize('mST.mol2','testmsT')
os.system(f'cp {os.path.join(lp["mol2"],"DGEBA.mol2")} .')
p.GAFFParameterize('DGEBA.mol2','testDGEBA',extra_antechamber_params='-eq 1 -pl 10')
s1=pm.load_file('testmsT.top')
s2=pm.load_file('testDGEBA.top')
print(s1)
print(s2)
s3=s1*10+s2*10
s3.save('tens.top')



