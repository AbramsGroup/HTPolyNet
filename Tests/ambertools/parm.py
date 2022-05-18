import parmed
import argparse as ap
parser=ap.ArgumentParser()
parser.add_argument('-acrd',default='')
parser.add_argument('-atop',default='')
parser.add_argument('-gro',default='')
parser.add_argument('-top',default='')
a=parser.parse_args()
file=parmed.load_file(f'{a.atop}', xyz=f'{a.acrd}')
file.save(a.gro,overwrite=True)
file.save(f'complete-{a.top}',overwrite=True)
itp='.'.join(a.top.split('.')[0:-1])+'.itp'
file.save(a.top,parameters=itp,overwrite=True)
