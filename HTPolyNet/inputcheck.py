from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.configuration import Configuration
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.command import Command
import os

def input_check(args):
    lib='./lib/molecules'
    C=Configuration.read(args.config,plot_reaction_network=False)
    icdict={x['molecule']:x['count'] for x in C.initial_composition}
    natoms=0
    tmass=0.0
    mmass={}
    for mname,M in C.molecules.items():
        matoms=0
        mmass[mname]=0.0
        if os.path.exists(os.path.join(lib,'parameterized',f'{mname}.top')):
            M.TopoCoord=TopoCoord(grofilename=os.path.join(lib,'parameterized',f'{mname}.gro'),topfilename=os.path.join(lib,'parameterized',f'{mname}.top'))
            matoms=M.TopoCoord.Coordinates.A.shape[0]
            mmass[mname]=M.TopoCoord.Topology.total_mass()
        elif os.path.exists(os.path.join(lib,'inputs',f'{mname}.mol2')):
            c=Coordinates.read_mol2(os.path.join(lib,'inputs',f'{mname}.mol2'))
            matoms=c.A.shape[0]
        elif os.path.exists(os.path.join(lib,'inputs',f'{mname}.pdb')):
            # print(os.path.join(lib,f"{mname}.{fmt}"))
            out,err=Command(f'grep -c ^ATOM {os.path.join(lib,"inputs",f"{mname}.pdb")}').run(ignore_codes=[1])
            matoms=int(out)
            out,err=Command(f'grep -c ^HETATM {os.path.join(lib,"inputs",f"{mname}.pdb")}').run(ignore_codes=[1])
            matoms+=int(out)
        if matoms>0 and mname in icdict:
            natoms+=icdict[mname]*matoms
            tmass+=icdict[mname]*mmass[mname]
            print(f'Molecule {mname}: {matoms} atoms, {icdict[mname]} molecules')
    print(f'{args.config}: {natoms} atoms in initial system.')
    if tmass>0.0:
        for mname in icdict:
            wtpct=100.0*icdict[mname]*mmass[mname]/tmass
            print(f'Molecule {mname}: {wtpct:0.2f} wt-%')
