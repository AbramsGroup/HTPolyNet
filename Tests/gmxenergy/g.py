from HTPolyNet.plot import trace
# pref='3-relax-stage-{n}-{ens}'
# l=[]
# for n in range(10):
#     for ens in ['min','nvt','npt']:
#         l.append(pref.format(n=n,ens=ens))
l='test'
trace(['Potential'],[l],outfile='test.png')
