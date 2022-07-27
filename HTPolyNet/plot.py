from HTPolyNet.gromacs import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import logging
import networkx as nx
from datetime import datetime
logger=logging.getLogger(__name__)

# prevents "RuntimeError: main thread is not in main loop" tk bug
plt.switch_backend('agg')

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
#            print(df.tail(1).to_string())
            xshift=df.tail(1).iloc[0,0]
#        print(f'xshift {xshift}')
        data=gmx_energy_trace(edr,[qty],xshift=xshift)
        # print(data.columns)
        lastchkpt=0
        if len(chkpt)>0:
            lastchkpt=chkpt[-1]
        chkpt.append(data.shape[0]+lastchkpt)
        df=pd.concat((df,data),ignore_index=True)
    # print(df.columns)
    fig,ax=plt.subplots(1,1,figsize=size)
    nseg=len(chkpt)
    beg=0
    for c in df.columns[1:]:
        for seg in range(nseg):
            ax.plot(df.iloc[beg:chkpt[seg],0],df[c].iloc[beg:chkpt[seg]],label=(c if seg==0 else None),color=cmap(seg/nseg))
            beg=chkpt[seg]
    if avgafter>0:
        pass
    else:
        avgafter=df['time (ps)']/2
    sdf=df[df['time (ps)']>avgafter]
    avg=sdf[c].mean()
    ax.plot(df.iloc[:,0],[avg]*df.shape[0],'k-',alpha=0.3,label=f'{avg:0.2f}')
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
    df={}
    for logfile in logfiles:
        with open(logfile,'r') as f:
            lines=f.read().split('\n')
        # extracting lines of this format:
        # 2022-06-06 10:49:38,673 Iter 52 current conversion: 0.927 (927/1000)
        data={}
        data['time']=[]
        data['conv']=[]
        data['iter']=[]
        for l in lines:
            # CURE iteration 1/150 begins in
            tok=l.split()
            if len(tok)>2:
                if len(tok)>6:
                    if tok[2]=='CURE' and tok[5]=='begins' and tok[4][0:2]=='1/':
                        data['time'].append(datetime.strptime(' '.join(tok[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                        data['conv'].append(0.0)
                        data['iter'].append(0)
                if tok[2]=='Iter' and tok[4]=='current':
                    data['time'].append(datetime.strptime(' '.join(tok[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                    data['conv'].append(float(tok[6]))
                    data['iter'].append(int(tok[3]))
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
        ax[1].plot(df[logfile]['elapsed'],df[logfile]['iter'])
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)
