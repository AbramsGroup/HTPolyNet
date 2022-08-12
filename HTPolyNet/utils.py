
from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.gromacs import gmx_energy_trace
import pandas as pd
from scipy.constants import physical_constants
import os
import numpy as np

def density_from_gro(gro,mollib='./lib/molecules/parameterized',units='SI'):
    C=Coordinates.read_gro(gro)
    resnames=list(set(C.A['resName']))
    templates={x:TopoCoord.from_top_gro(os.path.join(mollib,f'{x}.top'),os.path.join(mollib,f'{x}.gro')) for x in resnames}
    mass=0.0
    counts={}
    for nm,T in templates.items():
        if not nm in counts:
            counts[nm]={}
            for an in T.Topology.D['atoms']['atom']:
                counts[nm][an]=0
        these=C.A[C.A['resName']==nm]['atomName']
        for n in these:
            counts[nm][n]+=1
    for nm,adict in counts.items():
        T=templates[nm].Topology.D['atoms']
        for aname,c in adict.items():
            mass+=c*T[T['atom']==aname]['mass'].values[0]
    # for i,r in C.A.iterrows():
    #     resname=r['resName']
    #     tadf=templates[resname].Topology.D['atoms']
    #     atomName=r['atomName']
    #     m=tadf[tadf['atom']==atomName]['mass'].values[0]
    #     # print(f'{atomName} {m}')
    #     mass+=m
    volume=np.prod(C.box.diagonal())
    fac=1.0
    if units=='SI':
        mfac=physical_constants['atomic mass constant'][0]
        lfac=1.e-9
        vfac=lfac**3
        fac=mfac/vfac
    return fac*mass/volume

_system_dirs=['densification','precure',r'iter-{iter:d}','capping','postcure']
_md_ensembles={'min':[],'nvt':['Temperature'],'npt':['Temperature','Density']}
_indir_pfx={}
_indir_pfx['densification']=[r'densified-{ens:s}']
_indir_pfx['precure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
_indir_pfx['iter-n']=[r'1-cure_drag-stage-{stage:d}-{ens:s}',r'3-cure_relax-stage-{stage:d}-{ens:s}',r'4-cure_equilibrate-{ens:s}']
_indir_pfx['capping']=[r'7-cap_relax-stage-{stage:d}-{ens:s}',r'8-cap_equilibrate-{ens:s}']
_indir_pfx['postcure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
def _concat_from_edr(df,edr,names,add_if_missing=[('Density',0.0)]):
    xshift=0.0
    if not df.empty: xshift=df.iloc[-1]['time (ps)']
    if len(names)==0: return df,xshift
    # print(f'{add_if_missing}')
    this_df=gmx_energy_trace(edr,names,xshift=xshift)
    for m in add_if_missing:
        nm,val=m
        if not nm in this_df:
            this_df[nm]=np.ones(this_df.shape[0],dtype=float)*val
    df=pd.concat((df,this_df),ignore_index=True)
    return df,xshift

def density_evolution(proj_dir):
    if not os.path.exists(proj_dir): return
    sysd=os.path.join(proj_dir,'systems')
    df=pd.DataFrame()
    xshift=0.0
    transition_times=[xshift]
    interval_labels=[]
    markers=[]
    for subd in _system_dirs:
        # print(f'subd {subd}')
        if subd==r'iter-{iter:d}' or subd=='capping':
            # loooking for iterations
            markers.append(df.iloc[-1]['time (ps)'])
            # print(iter_mark0)
            iter=1
            iter_subd=os.path.join(sysd,subd if subd=='capping' else subd.format(iter=iter))
            while os.path.exists(iter_subd):
                dirkey='iter-n' if subd==r'iter-{iter:d}' else 'capping'
                # print(f'iter_subd {iter_subd} dirkey {dirkey}')
                for pfx in _indir_pfx[dirkey]:
                    # print(f'pfx {pfx}')
                    if r'stage' in pfx:
                        stg=1
                        stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,ens=x))+r'.edr') for x in _md_ensembles])
                        # print(f'stg_present {stg_present}')
                        while stg_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(iter_subd,pfx.format(stage=stg,ens=ens))
                                gro=edr_pfx+r'.gro'
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                                if os.path.exists(edr_pfx+r'.edr'):
                                    density=0.0
                                    if ens!='npt':
                                        if os.path.exists(gro):
                                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                                            # print(f'density from gro {density}')
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys,add_if_missing=[('Density',density)])
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                            stg+=1
                            stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,ens=x))) for x in _md_ensembles])
                    if r'equil' in pfx:
                        eq_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(ens=x))+r'.edr') for x in _md_ensembles])
                        if eq_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                                if os.path.exists(edr_pfx+r'.edr'):
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys)
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                iter+=1        
                iter_subd=os.path.join(sysd,subd.format(iter=iter)) if r'iter' in subd else 'I_BET_THIS_FILE_DNE'
            markers.append(df.iloc[-1]['time (ps)'])
        else:
            this_subd=os.path.join(sysd,subd)
            for pfx in _indir_pfx[subd]:
                # print(f'this_subd {this_subd} pfx {pfx}')
                if 'ens' in pfx:
                    for ens,qtys in _md_ensembles.items():
                        edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                        gro=edr_pfx+r'.gro'
                        # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                        if os.path.exists(edr_pfx+r'.edr'):
                            density=0.0
                            if ens=='nvt' and os.path.exists(gro):
                                density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                            df,xshift=_concat_from_edr(df,edr_pfx,qtys,add_if_missing=[('Density',density)])
                            transition_times.append(xshift)
                            interval_labels.append([edr_pfx])
                else:
                    edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                    gro=edr_pfx+r'.gro'
                    # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                    if os.path.exists(edr_pfx+r'.edr'):
                        density=0.0
                        if os.path.exists(gro):
                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                        df,xshift=_concat_from_edr(df,edr_pfx,['Temperature'],add_if_missing=[('Density',density)])
                        transition_times.append(xshift)
                        interval_labels.append([edr_pfx])
    return df,transition_times,markers,interval_labels
