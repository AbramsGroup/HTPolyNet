"""

.. module:: gromacs
   :synopsis: methods for handling the gmx suite of executables
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import pandas as pd
import numpy as np
from itertools import product
from collections import namedtuple
from HTPolyNet.command import Command
import HTPolyNet.software as sw
logger=logging.getLogger(__name__)


def insert_molecules(composition,boxSize,outName,inputs_dir='.',**kwargs):
    """insert_molecules launcher for gmx insert-molecules

    :param composition: dictionary of molecule_name:count
    :type composition: dict
    :param boxSize: size of box as 3 floats or 1 float (if cubic box)
    :type boxSize: list(3,float) or float
    :param outName: basename of output files
    :type outName: str
    :param inputs_dir: directory to search for input structures, defaults to '.'
    :type inputs_dir: str, optional
    :raises Exception: if gmx-insert molecules fails to insert the requested number of molecules
    """
    if type(boxSize)==int:
        boxSize=float(boxSize)
    if type(boxSize)==float:
        boxSize=[boxSize]*3
    box=' '.join([f'{x:.8f}' for x in boxSize])
    scale=kwargs.get('scale',0.4) # our default vdw radius scaling
    for name,num in composition.items():  # composition determines order
        ci=os.path.join(inputs_dir,f'{name}.gro')
        if os.path.isfile(f'{outName}.gro'):
            logger.debug(f'gmx insert-molecules inserts into existing {outName}.gro')
            ''' final gro file exists; we must insert into it '''
            c=Command(f'{sw.gmx} {sw.gmx_options} insert-molecules',f=f'{outName}.gro',ci=ci,nmol=num,o=outName,box=box,scale=scale)
        else:
            ''' no final gro file yet; make it '''
            c=Command(f'{sw.gmx} {sw.gmx_options} insert-molecules',ci=ci,nmol=num,o=outName,box=box,scale=scale)
        out,err=c.run()
        out+=err
        # logger.debug(f'Output of "{sw.gmx} insert-molecules"')
        # for ln in out.split('\n'):
        #     logger.debug(ln)
        if 'Added' in out:
            outlines=out.split('\n')
            for ol in outlines:
                if 'Added' in ol:
                    tokens=ol.split()
                    numadded=int(tokens[1])
                    break
            if numadded!=num:
                msg=f'{sw.gmx} insert-molecules did not add enough {name}.\nOnly {numadded} out of {num} were placed.\nIncrease your initial boxsize or lower your initial density.'
                for l in msg.split('\n'):
                    logger.debug(l)
                raise Exception('need bigger box')

def grompp_and_mdrun(gro='',top='',out='',mdp='',boxSize=[],single_molecule=False,**kwargs):
    """grompp_and_mdrun launcher for grompp and mdrun

    :param gro: input gro file, defaults to ''
    :type gro: str, optional
    :param top: input top file, defaults to ''
    :type top: str, optional
    :param out: output file basename, defaults to ''
    :type out: str, optional
    :param mdp: input mdp file, defaults to ''
    :type mdp: str, optional
    :param boxSize: explicit box size, defaults to []
    :type boxSize: list, optional
    :param single_molecule: if true, a single-molecule system is simulated, defaults to False
    :type single_molecule: bool, optional
    """
    logger.debug(kwargs)
    quiet=kwargs.get('quiet',True)
    ignore_codes=kwargs.get('ignore_codes',[-11])
    maxwarn=kwargs.get('maxwarn',4)
    mdrun_options=kwargs.get('mdrun_options',{})
    for option in ['rdd','dds','dlb','npme','nt','ntpmi','ntomp','ntomp_pme','nb','tunepme','pme','pmefft','bonded','update']:
        if option in kwargs and not option in mdrun_options:
            mdrun_options[option]=kwargs[option]

    #tunepme=kwargs.get('tunepme','yes')
    if gro=='' or top=='' or out=='' or mdp=='':
        raise Exception('grompp_and_mdrun requires gro, top, out, and mdp filename prefixes.')
    infiles=[f'{gro}.gro',f'{top}.top',f'{mdp}.mdp']
    assert all([os.path.exists(x) for x in infiles])
    if len(boxSize)>0:
        logger.debug(f'Resizing to {boxSize}')
        c=Command(f'{sw.gmx} {sw.gmx_options} editconf',f=f'{gro}.gro',o=gro,
                     box=' '.join([f'{x:.8f}' for x in boxSize]))
        c.run(quiet=quiet)
    # nsteps=kwargs.get('nsteps',-2)
    c=Command(f'{sw.gmx} {sw.gmx_options} grompp',f=f'{mdp}.mdp',c=f'{gro}.gro',p=f'{top}.top',o=f'{out}.tpr',maxwarn=maxwarn)
    c.run(quiet=quiet)
    if single_molecule:
        c=Command(f'{sw.mdrun_single_molecule}',deffnm=out,**mdrun_options)
    else:
        c=Command(f'{sw.mdrun}',deffnm=out,**mdrun_options)
    c.run(quiet=quiet,ignore_codes=ignore_codes)
    if os.path.exists(f'{out}.gro'):
        pass
        # logger.info(f'grompp_and_mdrun completed.  Check {gro}.gro.')
    else:
        logger.error(f'{sw.mdrun} ended prematurely; {out}.gro not found.')
        raise Exception(f'{sw.mdrun} ended prematurely; {out}.gro not found.')

def get_energy_menu(edr,**kwargs):
    """get_energy_menu gets the menu provided by 'gmx energy' when applied to a particular edr file

    :param edr: name of edr file
    :type edr: str
    :return: menu dictionary
    :rtype: dict
    """
    assert os.path.exists(edr+'.edr'),f'Error: {edr} not found'
    with open('_menugetter_','w') as f:
        f.write('\n')
    c=Command(f'{sw.gmx} {sw.gmx_options} energy -f {edr}.edr -o {edr}-out.xvg -xvg none < _menugetter_ > _menu_ 2>&1')
    c.run(ignore_codes=[1,2])
    with open('_menu_','r') as f:
        lines=f.read().split('\n')
    topline=lines.index('End your selection with an empty line or a zero.')
    botline=topline
    while len(lines[botline])>0:
        botline+=1
    menulines=lines[topline+2:botline]
    m=''.join(menulines)
    n=' '.join(m.split()).split()
    num,label=n[::2],n[1::2]
    menu={l:n for l,n in zip(label,num)}
    os.remove('_menu_')
    os.remove('_menugetter_')
    return menu

def gmx_energy_trace(edr,names=[],report_averages=False,keep_files=False,**kwargs):
    """Generate traces of data in edr file

    :param edr: name of edr file
    :type edr: str
    :param names: list of data names, defaults to []
    :type names: list
    :param report_averages:flag to indicate if averages of all data are to be computed here, default False
    :type report_averages: bool, optional
    :param keep_files: flag indicating caller would like to keep the raw input and output files for gmx energy, default False
    :type keep_files: bool
    :return: dataframe of traces
    :rtype: pandas DataFrame
    """
    assert os.path.exists(edr+'.edr'),f'Error: {edr}.edr not found'
    assert len(names)>0,f'Nothing to plot'
    xshift=kwargs.get('xshift',0)
    menu=get_energy_menu(edr)
    if not any([i in menu for i in names]):
        logger.debug(f'None of {names} in menu {menu}')
        return pd.DataFrame()
    logger.debug(f'report_averages? {report_averages} {names}')
    namvals=[]
    with open(f'{edr}-gmx.in','w') as f:
        for i in names:
            if i in menu:
                f.write(f'{menu[i]}\n')
                namvals.append((i,int(menu[i])))
        f.write('\n')
    namvals.sort(key=lambda x: x[1])
    c=Command(f'{sw.gmx} {sw.gmx_options} energy -f {edr}.edr -o {edr}-out.xvg -xvg none < {edr}-gmx.in')
    c.run()
    colnames=[x[0] for x in namvals]
    colnames.insert(0,'time(ps)')
    data=pd.read_csv(f'{edr}-out.xvg',sep='\s+',header=None,names=colnames)
    data.iloc[:,0]+=xshift
    ndata=data.shape[0]
    if report_averages:
        for i in names:
            data[f'Running-average-{i}']=data[i].expanding(1).mean()
            data[f'Rolling-10-average-{i}']=data[i].rolling(window=ndata//10).mean()
            for ln in data.iloc[-1][[i,f'Running-average-{i}',f'Rolling-10-average-{i}']].to_string(float_format='{:.2f}'.format).split('\n'):
                logger.info(f'{ln}')
    if not keep_files:
        os.remove(f'{edr}-out.xvg')
        os.remove(f'{edr}-gmx.in')
    return data       

# make a bunch of 3-character filename prefixes so parallel invocations don't collide
_abc='abcdefghijklmnopqrstuwxyz'
_fnames=[''.join(i) for i in product(_abc,_abc,_abc)]
def gromacs_distance(idf,gro,new_column_name='r',pfx='tmp',force_recalculate=False,keep_files=False):
    """Use 'gmx distance' to measure interatomic distances

    :param idf: dataframe of atom indexes in pairs ['ai','aj']
    :type idf: pandas DataFrame
    :param gro: name of gromacs input file for 'gmx distance' to use
    :type gro: str
    :param new_column_name: name of column in idf where distances are stored, default 'r'
    :type new_column_name: str, optional
    :param force_recalculate: flag to force calculation of distances even if a distance column exists in idf, default False
    :type force_recalculate: boolean, optional
    :param keep_files: flag indicating caller would like to keep the raw input and output files for gmx energy default False
    :type keep_files: bool, optional
    :return: list of distances parallel to idf columns
    :rtype: numpy.ndarray
    """
    if type(idf)==tuple: # this is being called in parallel
        i,idf=idf # unpack index and actual data frame
        pfx=_fnames[i]
        logger.debug(f'packet {i} using fname {pfx}; dataframe size {idf.shape[0]}')
    npair=idf.shape[0]
    # logger.debug(f'idf dtype {idf["ai"].dtype}')
    if npair==0 or ('r' in idf and not force_recalculate):
        return idf
    ''' create the index file '''
    with open(f'{pfx}.ndx','w') as f:
        f.write('[ bonds-1 ]\n')
        idf[['ai','aj']].to_csv(f,sep=' ',header=False,index=False)
        # logger.debug('wrote tmp.ndx')
    ''' create the user-input file '''
    with open (f'{pfx}.in','w') as f:
        f.write('0\n\n')
    c=Command(f'{sw.gmx} {sw.gmx_options} distance -f {gro} -n {pfx}.ndx -oall {pfx}.xvg -xvg none < {pfx}.in')
    c.run()
    ''' read the xvg file that 'gmx distance' creates '''
    with open (f'{pfx}.xvg','r') as f:
        datastr=f.read().split()
    datastr=datastr[1:]
    # logger.debug(f'0 {datastr[0]} - {len(datastr)-1} {datastr[-1]}')
    data=np.array([float(x) for x in datastr])
    nd=len(data)
    ''' clean up '''
    os.remove(f'{pfx}.ndx')
    os.remove(f'{pfx}.in')
    os.remove(f'{pfx}.xvg')
    ''' add the lengths back to the input dataframe '''
    idf[new_column_name]=data
    assert npair==nd
    return idf

def mdp_to_dict(mdp_filename):
    with open(mdp_filename,'r') as f:
        lines=f.read().split('\n')
    all_dict={}
    for l in lines:
        if len(l)>0:
            if '=' in l:
                k,v=l.split('=')
                k=k.strip()
                v=v.strip()
                all_dict[k]=v
            else:
                logger.debug(f'line {l} in {mdp_filename} skipped')
    return all_dict

def mdp_get(mdp_filename,key):
    all_dict=mdp_to_dict(mdp_filename)
    return all_dict.get(key,'NOT FOUND!')

def mdp_modify(mdp_filename,opt_dict,new_filename=None,add_if_missing=True):
    """Modify a gromacs mdp file

    :param mdp_filename: name of mdp file to modify; overwritten if new_filename==None
    :type mdp_filename: str
    :param opt_dict: keyword:value dictionary of mdp options
    :type opt_dict: dict
    :param new_filename: name of outputile, defaults to None
    :type new_filename: str, optional
    :param add_if_missing: Flag indicating whether to insert key:value into mdp file if not already there, defaults to True
    :type add_if_missing: bool, optional
    """
    all_dict=mdp_to_dict(mdp_filename)
    # logger.debug(f'mdp_modify: all_dict: {all_dict}')
    for k,v in opt_dict.items():
        if not k in all_dict:
            if add_if_missing:
                logger.debug(f'adding {k} = {v} to {mdp_filename}')
                all_dict[k]=v
        else:
            all_dict[k]=v
    if new_filename:
        with open(new_filename,'w') as f:
            for k,v in all_dict.items():
                f.write(f'{k} = {v}\n')
        # logger.debug(f'wrote {new_filename}.')
    else:
        with open(mdp_filename,'w') as f:
            for k,v in all_dict.items():
                f.write(f'{k} = {v}\n')
        # logger.debug(f'wrote {mdp_filename}.')

def gmx_traj_info(trr):
    Result=namedtuple('gmx_check','nframes time')
    c=Command(f'{sw.gmx} {sw.gmx_options} check -f {trr}')
    out,err=c.run()
    out+=err
    for l in out.split('\n'):
        tok=l.split()
        if len(tok)==0:
            continue
        if tok[0]=='Step':
            # print(tok)
            nframes=int(tok[1])
            ts_ps=float(tok[2])
            break
    result=Result(nframes,(nframes-1)*ts_ps)
    return result
    
def gro_from_trr(pfx,nzero=2,b=0,outpfx=''):
    if not outpfx:
        outpfx=pfx
    with open('tmp.in','w') as f:
        f.write('0\n')
    c=Command(f'{sw.gmx} {sw.gmx_options} trjconv -f {pfx}.trr -s {pfx}.tpr -o {outpfx}.gro -sep -nzero {nzero} -b {b} < tmp.in')
    out,err=c.run()
    os.remove('tmp.in')
    out+=err