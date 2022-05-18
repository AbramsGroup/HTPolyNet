from HTPolyNet.gromacs import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd

def trace(qtys,edrs,outfile='plot.png',**kwargs):
    df=pd.DataFrame()
    cmapname=kwargs.get('colormap','plasma')
    cmap=cm.get_cmap(cmapname)
    xshift=0.0
    chkpt=[]
    for edr in edrs:
        if not df.empty:
#            print(df.tail(1).to_string())
            xshift=df.tail(1).iloc[0,0]
#        print(f'xshift {xshift}')
        data=gmx_energy_trace(edr,qtys,xshift=xshift)
        lastchkpt=0
        if len(chkpt)>0:
            lastchkpt=chkpt[-1]
        chkpt.append(data.shape[0]+lastchkpt)
        df=pd.concat((df,data),ignore_index=True)
    fig,ax=plt.subplots(1,1,figsize=(8,6))
    nseg=len(chkpt)
    beg=0
    for c in df.columns[1:]:
        for seg in range(nseg):
            ax.plot(df.iloc[beg:chkpt[seg],0],df[c].iloc[beg:chkpt[seg]],label=c,color=cmap(seg/nseg))
            beg=chkpt[seg]
    plt.legend()
    plt.savefig(outfile)
