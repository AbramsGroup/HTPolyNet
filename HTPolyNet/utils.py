from HTPolyNet.coordinates import Coordinates
from HTPolyNet.topocoord import TopoCoord
from HTPolyNet.gromacs import gmx_energy_trace
from scipy.constants import physical_constants
import pandas as pd
import os
import numpy as np
import networkx as nx
import logging

logger=logging.getLogger(__name__)

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
    volume=np.prod(C.box.diagonal())
    fac=1.0
    if units=='SI':
        mfac=physical_constants['atomic mass constant'][0]
        lfac=1.e-9
        vfac=lfac**3
        fac=mfac/vfac
    return fac*mass/volume

_system_dirs=['densification','precure',r'iter-{iter:d}','capping','postcure']
_md_ensembles={'nvt':['Temperature','Potential'],'npt':['Temperature','Potential','Density'],'default':['Temperature','Potential']}
_indir_pfx={}
_indir_pfx['densification']=[r'densified-{ens:s}',r'densified-repeat-{repeat:d}-{ens:s}']
_indir_pfx['precure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
_indir_pfx['iter-n']=[r'1-cure_drag-stage-{stage:d}-{ens:s}',r'3-cure_relax-stage-{stage:d}-{ens:s}',r'4-cure_equilibrate-{ens:s}']
_indir_pfx['capping']=[r'7-cap_relax-stage-{stage:d}-{ens:s}',r'8-cap_equilibrate-{ens:s}']
_indir_pfx['postcure']=[r'preequilibration-{ens:s}',r'annealed',r'postequilibration-{ens:s}']
def _concat_from_edr(df,edr,names,add=[],add_if_missing=[('Density',0.0)]):
    xshift=0.0
    if not df.empty: xshift=df.iloc[-1]['time(ps)']
    if len(names)==0: return df,xshift
    # print(f'{add_if_missing}')
    this_df=gmx_energy_trace(edr,names,xshift=xshift)
    for m in add:
        nm,val=m
        this_df[nm]=np.ones(this_df.shape[0],dtype=type(val))*val
    for m in add_if_missing:
        nm,val=m
        if not nm in this_df:
            this_df[nm]=np.ones(this_df.shape[0],dtype=type(val))*val
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
    nbonds=0
    for subd in _system_dirs:
        # print(f'subd {subd}')
        if subd==r'iter-{iter:d}' or subd=='capping' or subd=='densification':
            # loooking for iterations
            # print(iter_mark0)
            iter=1
            iter_subd_rel=subd if subd in ['capping','densification'] else subd.format(iter=iter)
            iter_subd=os.path.join(sysd,iter_subd_rel)
            while os.path.exists(iter_subd):
                logger.info(f'Reading edrs in {iter_subd_rel}...')
                if r'iter' in subd:
                    markers.append(df.iloc[-1]['time(ps)'])
                dirkey='iter-n' if subd==r'iter-{iter:d}' else subd
                if subd==r'iter-{iter:d}':
                    if os.path.exists(os.path.join(iter_subd,'2-cure_update-bonds.csv')):
                        bdf=pd.read_csv(os.path.join(iter_subd,'2-cure_update-bonds.csv'),sep='\s+',header=0,index_col=None)
                        nbonds+=bdf.shape[0]
                    else:
                        if os.path.exists(os.path.join(iter_subd,'6-cap_update-bonds.csv')):
                            bdf=pd.read_csv(os.path.join(iter_subd,'6-cap_update-bonds.csv'),sep='\s+',header=0,index_col=None)
                            nbonds+=bdf.shape[0]
                # print(f'iter_subd {iter_subd} dirkey {dirkey}')
                for pfx in _indir_pfx[dirkey]:
                    # print(f'pfx {pfx}')
                    if r'stage' in pfx or r'repeat' in pfx:
                        stg=1
                        stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=x))+r'.edr') for x in _md_ensembles])
                        # print(f'{iter_subd} {stg} {pfx} stg_present {stg_present}')
                        while stg_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=ens))
                                gro=edr_pfx+r'.gro'
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                                if os.path.exists(edr_pfx+r'.edr'):
                                    density=0.0
                                    if ens!='npt':
                                        if os.path.exists(gro):
                                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                                            # print(f'density from gro {density}')
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                            stg+=1
                            stg_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(stage=stg,repeat=stg,ens=x))+r'.edr') for x in _md_ensembles])
                            # print(f'->{iter_subd} {stg} {pfx} stg_present {stg_present}')
                    elif r'equil' in pfx or r'densif' in pfx:
                        eq_present=any([os.path.exists(os.path.join(iter_subd,pfx.format(ens=x))+r'.edr') for x in _md_ensembles])
                        # print(f'{iter_subd} eq_present? {eq_present}')
                        if eq_present:
                            for ens,qtys in _md_ensembles.items():
                                edr_pfx=os.path.join(iter_subd,pfx.format(ens=ens))
                                gro=edr_pfx+r'.gro'
                                # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'),ens,qtys)
                                if os.path.exists(edr_pfx+r'.edr'):
                                    density=0.0
                                    if ens!='npt':
                                        if os.path.exists(gro):
                                            density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                                    df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                                    transition_times.append(xshift)
                                    interval_labels.append([edr_pfx])
                iter+=1       
                iter_subd_rel=subd.format(iter=iter)
                iter_subd=os.path.join(sysd,iter_subd_rel) if r'iter' in subd else 'I_BET_THIS_FILE_DNE'
            # if r'iter' in subd: markers.append(df.iloc[-1]['time (ps)'])
        else:
            this_subd=os.path.join(sysd,subd)
            logger.info(f'Reading edrs in {subd}...')

            for pfx in _indir_pfx[subd]:
                # print(f'this_subd {this_subd} pfx {pfx}')
                if 'ens' in pfx:
                    for ens,qtys in _md_ensembles.items():
                        edr_pfx=os.path.join(this_subd,pfx.format(ens=ens))
                        gro=edr_pfx+r'.gro'
                        # print(edr_pfx,os.path.exists(edr_pfx+r'.edr'))
                        # print(gro,os.path.exists(gro))
                        if os.path.exists(edr_pfx+r'.edr'):
                            density=0.0
                            if ens=='nvt' and os.path.exists(gro):
                                density=density_from_gro(gro,os.path.join(proj_dir,'molecules/parameterized'))
                            df,xshift=_concat_from_edr(df,edr_pfx,qtys,add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
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
                        df,xshift=_concat_from_edr(df,edr_pfx,_md_ensembles['default'],add=[('nbonds',nbonds)],add_if_missing=[('Density',density)])
                        transition_times.append(xshift)
                        interval_labels.append([edr_pfx])
    return df,transition_times,markers,interval_labels

def graph_from_bondsfile(bondsfile):
    G=nx.Graph()
    df=pd.read_csv(bondsfile,header=0,index_col=None,sep='\s+')
    for i,r in df.iterrows():
        G.add_edge(r['mi'],r['mj'])
    return G

def mwbxl(G:nx.Graph,crosslinker='GMA',monomer='STY'):
    xG=G.copy()
    # traverse edges and label hetero/homo participants
    for n in xG.nodes():
        xG.nodes[n]['neighbor_count']=0
    for u,v in xG.edges():
        n,m=xG.nodes[u],xG.nodes[v]
        n['neighbor_count']+=1
        m['neighbor_count']+=1
        if n['molecule_name']!=m['molecule_name']:
            n['monomer_type']='hetero'
            m['monomer_type']='hetero'
        else:
            n['monomer_type']='homo'
            m['monomer_type']='homo'
    wipe_us=[]
    for n in xG.nodes():
        if xG.nodes[n]['molecule_name']==crosslinker:
            wipe_us.append(n)
    xG.remove_nodes_from(wipe_us)
    chaintype=[]
    ch=nx.connected_components(xG)
    mwhist={}
    nb=0
    for c in ch:
        nhetero=0
        for n in c:
            node=xG.nodes[n]
            if node.get('monomer_type','none')=='hetero':
                nhetero+=1
        if nhetero==0:
            chaintype.append('isolated')
        elif nhetero==1:
            if len(c)==1:
                node=list(xG.nodes.values())[0]
                if node['neighbor_count']==1:
                    chaintype.append('dangling')
                elif node['neighbor_count']==2:
                    chaintype.append('bridging')
                    if len(c) not in mwhist:
                        mwhist[len(c)]=0
                    mwhist[len(c)]+=1
            else:
                chaintype.append('dangling')
        elif nhetero==2:
            chaintype.append('bridging')
            if len(c) not in mwhist:
                mwhist[len(c)]=0
            mwhist[len(c)]+=1
            nb+=1
    maxlen=max(list(mwhist.keys()))
    mw=[]
    counts=[]
    for n in range(maxlen):
        cnt=mwhist.get(n,0)
        mw.append(n)
        counts.append(cnt)
    df=pd.DataFrame({'n':mw,'counts':counts})
    return df

def clusters(G:nx.Graph):
    counts={}
    for c in sorted(nx.connected_components(G),key=len,reverse=True):
        size=len(c)
        if not size in counts:
            counts[size]=0
        counts[size]+=1
    sizes=[]
    numbers=[]
    for s,c in counts.items():
        sizes.append(s)
        numbers.append(c)
    df=pd.DataFrame({'sizes':sizes,'counts':numbers})
    df.sort_values(by='sizes',ascending=False,inplace=True)
    return df
