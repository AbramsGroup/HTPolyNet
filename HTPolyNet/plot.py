from HTPolyNet.gromacs import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import logging
import networkx as nx
from datetime import datetime
from glob import glob
# import argparse as ap

logger=logging.getLogger(__name__)

# prevents "RuntimeError: main thread is not in main loop" tk bug
plt.switch_backend('agg')
#matplotlib.use('agg')

def trace(qty,edrs,outfile='plot.png',**kwargs):
    # disable debug-level logging and above since matplotlib has a lot of debug statements
    logging.disable(logging.DEBUG)
    df=pd.DataFrame()
    cmapname=kwargs.get('colormap','plasma')
    size=kwargs.get('size',(8,6))
    yunits=kwargs.get('yunits',None)
    avgafter=kwargs.get('avgafter',0)
    cmap=cm.get_cmap(cmapname)
    xshift=0.0
    chkpt=[]
    for edr in edrs:
        if not df.empty:
            xshift=df.tail(1).iloc[0,0]
        data=gmx_energy_trace(edr,[qty],xshift=xshift)
        lastchkpt=0
        if len(chkpt)>0:
            lastchkpt=chkpt[-1]
        chkpt.append(data.shape[0]+lastchkpt)
        df=pd.concat((df,data),ignore_index=True)
    fig,ax=plt.subplots(1,1,figsize=size)
    nseg=len(chkpt)
    beg=0
    for c in df.columns[1:]:
        for seg in range(nseg):
            ax.plot(df.iloc[beg:chkpt[seg],0],df[c].iloc[beg:chkpt[seg]],label=(c if seg==0 else None),color=cmap(seg/nseg))
            beg=chkpt[seg]
        if 'avgafter' in kwargs:
            if avgafter>0:
                pass
            else:
                avgafter=df['time (ps)'].iloc[-1]/2
            sdf=df[df['time (ps)']>avgafter]
            avg=sdf[c].mean()
            ax.plot(df.iloc[:,0],[avg]*df.shape[0],'k-',alpha=0.3,label=f'{avg:0.2f}')
        else:
            avg=df[c].mean()
    if not yunits:
        plt.ylabel(qty)
    else:
        plt.ylabel(f'{qty} ({yunits})')
    plt.xlabel('time (ps)')
    plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)
    return avg

def network_graph(G,filename,**kwargs):
    logging.disable(logging.DEBUG)
    fig,ax=plt.subplots(1,1,figsize=(8,8))
    nx.draw_networkx(G,ax=ax,arrows=False,node_size=200)
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)
    
def cure_graph(logfiles,filename,**kwargs):
    xmax=kwargs.get('xmax',-1)
    logging.disable(logging.DEBUG)
    df:dict[pd.DataFrame]={}
    for logfile in logfiles:
        with open(logfile,'r') as f:
            lines=f.read().split('\n')
        print(f'read {len(lines)} lines from {logfile}')
        # extracting lines of these formats:
        # 2022-07-28 13:28:52,379 HTPolyNet.HTPolyNet.CURE INFO> ********** Connect-Update-Relax-Equilibrate (CURE) begins
        # 2022-07-28 13:13:37,328 HTPolyNet.HTPolyNet.CURE INFO> Current conversion: 0.26 (26/100)
        data={}
        data['time']=[]
        data['conv']=[]
        for l in lines:
            tok=l.split()
            ntoks=len(tok)
            if ntoks<8:
                continue
            if tok[2]!='HTPolyNet.runtime.CURE' and tok[2]!='HTPolyNet.HTPolyNet.CURE':
                continue
            if tok[3]!='INFO>':
                continue
            if tok[4]=='**********' and tok[5]=='Connect-Update-Relax-Equilibrate' and tok[7]=='begins':
                data['time'].append(datetime.strptime(' '.join(tok[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                data['conv'].append(0.0)
            if tok[4]=='Current' and tok[5]=='conversion:':
                data['time'].append(datetime.strptime(' '.join(tok[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                data['conv'].append(float(tok[6]))
        df[logfile]=pd.DataFrame(data)
        df[logfile]['elapsed']=(df[logfile]['time']-df[logfile].loc[0]['time']).astype(int)/1.e9/3600.0
    fig,ax=plt.subplots(1,2,sharex=True,figsize=(10,6))
    ax[0].set_ylim([0,1])
    ax[0].set_xlabel('runtime (h)')
    ax[0].set_ylabel('conversion')
    ax[1].set_xlabel('runtime (h)')
    ax[1].set_ylabel('iteration')
    if xmax>-1:
        ax[0].set_xlim([0,xmax])
    for logfile in logfiles:
        ax[0].plot(df[logfile]['elapsed'],df[logfile]['conv'])
        ax[1].plot(df[logfile]['elapsed'],df[logfile].index+1)
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)

def density_evolution():
    D=glob('proj-[0-9]') + glob('proj-[1-9][0-9]')
    for e in D:
        d=os.path.join(e,'systems')
        edrs=[os.path.join(d,'init/npt-1')]
        iter=1
        while os.path.exists(os.path.join(d,f'iter-{iter}')):
            print(iter)
            this_stgs=[]
            dstg=1
            while os.path.exists(os.path.join(d,f'iter-{iter}/1-drag-stage-{dstg}-npt.edr')):
                this_stgs.append(os.path.join(d,f'iter-{iter}/1-drag-stage-{dstg}-npt'))
                edrs.append(os.path.join(d,f'iter-{iter}/1-drag-stage-{dstg}-npt'))
                dstg+=1
            rstg=1
            while os.path.exists(os.path.join(d,f'iter-{iter}/3-relax-stage-{rstg}-npt.edr')):
                this_stgs.append(os.path.join(d,f'iter-{iter}/3-relax-stage-{rstg}-npt'))
                edrs.append(os.path.join(d,f'iter-{iter}/3-relax-stage-{rstg}-npt'))
                rstg+=1
            this_stgs.append(os.path.join(d,f'iter-{iter}/4-equilibrate-post'))
            edrs.append(os.path.join(d,f'iter-{iter}/4-equilibrate-post'))
            trace('Density',this_stgs,outfile=os.path.join(d,f'iter-{iter}/density-trace.png'),size=(16,4),yunits='kg/m3')
            iter+=1

        rstg=1
        this_stgs=[]
        while os.path.exists(os.path.join(d,f'postcure/5-relax-stage-{rstg}-npt.edr')):
            this_stgs.append(os.path.join(d,f'postcure/5-relax-stage-{rstg}-npt'))
            edrs.append(os.path.join(d,f'postcure/5-relax-stage-{rstg}-npt'))
            rstg+=1
        this_stgs.append(os.path.join(d,f'postcure/6-equilibrate-post'))
        edrs.append(os.path.join(d,f'postcure/6-equilibrate-post'))
        print(this_stgs)
        trace('Density',this_stgs,outfile=os.path.join(d,f'postcure/density-trace.png'),size=(16,4),yunits='kg/m3')
        #print(edrs)    
        trace('Density',edrs,outfile=f'{e}-density.png',size=(16,4),yunits='kg/m3')
    
