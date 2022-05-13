from HTPolyNet.gromacs import *
import matplotlib.pyplot as plt
import pandas as pd

def trace(qtys,edrs,outfile='plot.png',**kwargs):
    df=pd.DataFrame()
    xshift=0.0
    for edr in edrs:
        if not df.empty:
            print(df.tail(1).to_string())
            xshift=df.tail(1).iloc[0,0]
        print(f'xshift {xshift}')
        data=gmx_energy_trace(edr,qtys,xshift=xshift)
        df=pd.concat((df,data),ignore_index=True)
    fig,ax=plt.subplots(1,1,figsize=(8,6))
    for c in df.columns[1:]:
        ax.plot(df.iloc[:,0],df[c],label=c)
    plt.savefig(outfile)
