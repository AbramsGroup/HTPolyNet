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

def scatter(df,xcolumn,columns=[],outfile='plot.png',**kwargs):
    logging.disable(logging.DEBUG)
    cmapname=kwargs.get('colormap','plasma')
    size=kwargs.get('size',(8,6))
    yunits=kwargs.get('yunits',None)
    cmap=cm.get_cmap(cmapname)
    fig,ax=plt.subplots(1,1,figsize=size)
    ax.set_xlabel(xcolumn)
    for n in columns:
        ax.scatter(df[xcolumn],df[n],label=n)
    plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    logging.disable(logging.NOTSET)


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
                avgafter=df['time(ps)'].iloc[-1]/2
            sdf=df[df['time(ps)']>avgafter]
            avg=sdf[c].mean()
            ax.plot(df.iloc[:,0],[avg]*df.shape[0],'k-',alpha=0.3,label=f'{avg:0.2f}')
        else:
            avg=df[c].mean()
    if not yunits:
        plt.ylabel(qty)
    else:
        plt.ylabel(f'{qty} ({yunits})')
    plt.xlabel('time(ps)')
    plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)
    return avg

def global_trace(df,names,outfile='plot.png',transition_times=[],markers=[],interval_labels=[],y2names=[],**kwargs):
    # disable debug-level logging and above since matplotlib has a lot of debug statements
    default_units={'Temperature':'K','Pressure':'bar','Density':'kg/m^3','Potential':'kJ/mol'}
    units=kwargs.get('units',default_units)
    logging.disable(logging.DEBUG)
    size=kwargs.get('size',(16,4*len(names)))
    legend=kwargs.get('legend',False)
    fig,ax=plt.subplots(len(names),1,figsize=size)
    plt.xlabel('time(ps)')
    cmapname=kwargs.get('colormap','plasma')
    # yunits=kwargs.get('yunits',None)
    cmap=cm.get_cmap(cmapname)
    # print(f'in global_trace:\n{df.head().to_string()}')

    interval_times=[]
    if interval_labels:
        for i in range(1,len(transition_times)):
            interval_times.append((transition_times[i]+transition_times[i-1])/2)
    for l,t in zip(interval_labels,interval_times):
        logger.info(f'{t} {l}')
    assert len(interval_labels)==len(interval_times)
    L,R=-1,-1
    if len(markers)>1:
        L,R=markers[0],markers[-1]
        in_tt=[x for x in transition_times if L<x<R]
        marked_df=df[(df['time(ps)']>L)&(df['time(ps)']<R)]
        fig,ax=plt.subplots(len(names)*2,1,figsize=size)
        for i,colname in enumerate(names):
            out_ax=ax[0] if len(names)==1 else ax[i*2]
            in_ax=ax[1] if len(names)==1 else ax[i*2+1]
            out_ax.plot(df['time(ps)'],df[colname],label=colname)
            ylabel=colname
            if ylabel in units: ylabel+=f' ({units[ylabel]})'
            out_ax.set_ylabel(ylabel)
            out_ax.set_xlabel('time(ps)')
            if len(y2names)>i:
                out_ax2=out_ax.twinx()
                out_ax2.plot(df['time(ps)'],df[y2names[i]],label=y2names[i],color='black')
                ylabel=y2names[i]
                if ylabel in units: ylabel+=f' ({units[ylabel]})'
                out_ax2.set_ylabel(ylabel)
                out_ax2.set_xlabel('time(ps)')
            if len(transition_times)>0:
                colors=[cmap(i/len(transition_times)) for i in range(len(transition_times))]
                ylim=out_ax.get_ylim()
                out_ax.vlines(transition_times,ylim[0],ylim[1],color=colors,linewidth=0.75,alpha=0.5)
                # for x,l in zip(interval_times,interval_labels):
                #     out_ax.text(x,0.9*ylim[1],l,fontsize=8)
            in_ax.plot(marked_df['time(ps)'],marked_df[colname],label=colname)
            ylabel=colname
            if ylabel in units: ylabel+=f' ({units[ylabel]})'
            in_ax.set_xlabel('time(ps)')
            in_ax.set_yabel(ylabel)
            if len(y2names)>i:
                in_ax2=in_ax.twinx()
                in_ax2.plot(marked_df['time(ps)'],marked_df[y2names[i]],label=y2names[i],color='black')
                ylabel=y2names[i]
                if ylabel in units: ylabel+=f' ({units[ylabel]})'
                in_ax2.set_ylabel(ylabel)
                in_ax2.set_xlabel('time(ps)')
            if len(transition_times)>0:
                colors=[cmap(i/len(transition_times)) for i in range(len(transition_times))]
                ylim=in_ax.get_ylim()
                in_ax.vlines(in_tt,ylim[0],ylim[1],color=colors,linewidth=0.75,alpha=0.5)
                # for x,l in zip(interval_times,interval_labels):
                #     if L<x<R:
                #         out_ax.text(x,0.9*ylim[1],l,fontsize=8)
    else:
        fig,ax=plt.subplots(len(names),1,figsize=size)
        for i,colname in enumerate(names):
            the_ax=ax if len(names)==1 else ax[i]
            the_ax.plot(df['time(ps)'],df[colname],label=colname)
            ylabel=colname
            if ylabel in units: ylabel+=f' ({units[ylabel]})'
            the_ax.set_ylabel(ylabel)
            the_ax.set_xlabel('time(ps)')
            if len(y2names)>i:
                the_ax2=the_ax.twinx()
                the_ax2.plot(df['time(ps)'],df[y2names[i]],label=y2names[i],color='black')
                ylabel=y2names[i]
                if ylabel in units: ylabel+=f' ({units[ylabel]})'
                the_ax2.set_ylabel(ylabel)
            if len(transition_times)>0:
                colors=[cmap(i/len(transition_times)) for i in range(len(transition_times))]
                ylim=the_ax.get_ylim()
                the_ax.vlines(transition_times,ylim[0],ylim[1],color=colors,linewidth=0.5,alpha=0.5)

    # plt.xlabel('time(ps)')
    if legend:
        plt.legend()
    plt.savefig(outfile)
    plt.close(fig)
    # re-establish previous logging level
    logging.disable(logging.NOTSET)

def network_graph(G,filename,**kwargs):
    logging.disable(logging.DEBUG)
    arrows=kwargs.get('arrows',False)
    figsize=kwargs.get('figsize',(32,32))
    node_size=kwargs.get('node_size',200)
    cmap=cm.get_cmap('plasma')
    fig,ax=plt.subplots(1,1,figsize=figsize)
    ax.axis('off')
    molnames=list(set([n['molecule_name'] for k,n in G.nodes.items()]))
    nmolname=len(molnames)
    cx=[]
    for n in G.nodes.values():
        idx=molnames.index(n['molecule_name'])
        cx.append((float(idx)+1)/(nmolname+2))
    nx.draw_networkx(G,ax=ax,arrows=arrows,node_size=node_size,node_color=cx,cmap=cmap,with_labels=False)
    
    plt.savefig(filename)
    plt.close(fig)
    logging.disable(logging.NOTSET)

_template_1='2022-08-11 17:40:36,969 HTPolyNet.runtime.my_logger INFO> ********* Connect-Update-Relax-Equilibrate (CURE) begins **********'
_template_1_token_idx=[2,3,5,7]
_template_2='2022-08-11 17:43:35,512 HTPolyNet.runtime.do_cure INFO> Iteration 1 current conversion 0.210 or 63 bonds'
_template_2_token_idx=[2,3,6,7]
_template_2_data_idx={'iter':(int,5),'conv':(float,8),'nbonds':(int,10)}
def _token_match(l,template,pat_idx):
    if len(l.split())!=len(template.split()): return
    return all([l.split()[t]==template.split()[t] for t in pat_idx])
def _parse_data(dat,l,idx_dict):
    tok=l.split()
    for k,v in idx_dict.items():
        conv,s=v
        dat[k].append(conv(tok[s]))

def diagnostics_graphs(logfiles,filename,**kwargs):
    xmax=kwargs.get('xmax',-1)
    figsize=kwargs.get('figsize',(12,6))
    logging.disable(logging.DEBUG)
    df:dict[pd.DataFrame]={}
    for logfile in logfiles:
        bn,ex=os.path.splitext(logfile)
        with open(logfile,'r') as f:
            lines=f.read().split('\n')
        logger.info(f'read {len(lines)} lines from {logfile}')
        data={}
        data['time']=[]
        data['iter']=[]
        data['conv']=[]
        data['nbonds']=[]
        counter=0
        for l in lines:
            if _token_match(l,_template_1,_template_1_token_idx):
                # print('you should only see this once')
                counter+=1
                assert not counter>1
                data['time'].append(datetime.strptime(' '.join(l.split()[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                data['iter'].append(0)
                data['conv'].append(0.0)
                data['nbonds'].append(0)
            elif _token_match(l,_template_2,_template_2_token_idx):
                data['time'].append(datetime.strptime(' '.join(l.split()[0:2]),'%Y-%m-%d %H:%M:%S,%f'))
                # print('data tok',f'{l.split()}')
                _parse_data(data,l,_template_2_data_idx)
        # print('data',f'{data}')
        df[logfile]=pd.DataFrame(data)
        time_idx=list(df[logfile].columns).index('time')
        df[logfile]['elapsed']=(df[logfile]['time']-df[logfile].iloc[0,time_idx]).astype(int)/1.e9/3600.0
        df[logfile].to_csv(f'{bn}.csv',index=False,sep=' ',header=True)
    fig,ax=plt.subplots(1,2,sharex=True,figsize=figsize)
    ax[0].set_ylim([0,1])
    ax[0].set_xlabel('runtime (h)')
    ax[0].set_ylabel('conversion')
    ax[1].set_xlabel('runtime (h)')
    ax[1].set_ylabel('iteration')
    if xmax>-1:
        ax[0].set_xlim([0,xmax])
    for logfile in logfiles:
        ax[0].plot(df[logfile]['elapsed'],df[logfile]['conv'],label=logfile)
        ax[1].plot(df[logfile]['elapsed'],df[logfile].index+1,label=logfile)
    plt.legend()
    plt.savefig(filename)
    plt.close(fig)

    logging.disable(logging.NOTSET)

