"""

.. module:: gromacs
   :synopsis: methods for handling the gmx suite of executables
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import pandas as pd
import numpy as np
from collections import namedtuple

from pytrr import GroTrrReader
from HTPolyNet.command import Command
import HTPolyNet.software as sw
logger=logging.getLogger(__name__)
# mdp_library = {'liquid-densify':'liquid-densify-npt',
#                'sea':'sea-nvt',
#                'minimize-single-molecule':'minimize-single-molecule',
#                'minimize':'minimize',
#                'equilibrate':'equilibrate-npt',
#                'drag-minimize':'drag-minimize',
#                'drag-nvt':'drag-nvt',
#                'drag-npt':'drag-npt',
#                'relax-minimize':'relax-minimize',
#                'relax-nvt':'relax-nvt',
#                'relax-npt':'relax-npt'}

def insert_molecules(composition,boxSize,outName,**kwargs):
    ''' launcher for `gmx insert-molecules`
        monomers:  dictionary of Molecule instances keyed on molecule name
        composition: dictionary keyed on molecule name with value molecule count
        boxSize:  3-element list of floats OR a single float (cubic box)
        outName:  output filename basename.  If {outName}.gro exists,
                  insertions are made into it.
    '''
    if type(boxSize)==int:
        boxSize=float(boxSize)
    if type(boxSize)==float:
        boxSize=[boxSize]*3
    box=' '.join([f'{x:.8f}' for x in boxSize])
    scale=kwargs.get('scale',0.4) # our default vdw radius scaling
    for name,num in composition.items():  # composition determines order
        if os.path.isfile(f'{outName}.gro'):
            logger.debug(f'gmx insert-molecules inserts into existing {outName}.gro')
            ''' final gro file exists; we must insert into it '''
            c=Command(f'{sw.gmx} {sw.gmx_options} insert-molecules',f=f'{outName}.gro',ci=f'{name}.gro',nmol=num,o=outName,box=box,scale=scale)
        else:
            ''' no final gro file yet; make it '''
            c=Command(f'{sw.gmx} {sw.gmx_options} insert-molecules',ci=f'{name}.gro',nmol=num,o=outName,box=box,scale=scale)
        out,err=c.run()
        out+=err
        logger.debug(f'Output of "{sw.gmx} insert-molecules"')
        for ln in out.split('\n'):
            logger.debug(ln)
        if 'Added' in out:
            outlines=out.split('\n')
            for ol in outlines:
                if 'Added' in ol:
                    tokens=ol.split()
                    numadded=int(tokens[1])
                    break
            if numadded!=num:
                logger.error(f'{sw.gmx} insert molecules did not add enough {name}; only {numadded} out of {num} were placed.  Increase your initial boxsize or lower your initial density.')
                raise Exception('need bigger box')

def grompp_and_mdrun(gro='',top='',out='',mdp='',boxSize=[],**kwargs):
    ''' launcher for grompp and mdrun
        gro: prefix for input coordinate file
        top: prefix for input topology
        out: prefix for desired output files
        boxsize: (optional) desired box size; triggers editconf before grompp
    '''
    quiet=kwargs.get('quiet',True)
    ignore_codes=kwargs.get('ignore_codes',[-11])
    maxwarn=kwargs.get('maxwarn',4)
    rdd=kwargs.get('rdd',0)
    dds=kwargs.get('dds',0.8)
    dlb=kwargs.get('dlb','yes')
    tunepme=kwargs.get('tunepme','yes')
    if gro=='' or top=='' or out=='' or mdp=='':
        raise Exception('grompp_and_mdrun requires gro, top, out, and mdp filename prefixes.')
    infiles=[f'{gro}.gro',f'{top}.top',f'{mdp}.mdp']
    assert all([os.path.exists(x) for x in infiles])
    if len(boxSize)>0:
        c=Command(f'{sw.gmx} {sw.gmx_options} editconf',f=f'{gro}.gro',o=gro,
                     box=' '.join([f'{x:.8f}' for x in boxSize]))
        c.run(quiet=quiet)
    # nsteps=kwargs.get('nsteps',-2)
    c=Command(f'{sw.gmx} {sw.gmx_options} grompp',f=f'{mdp}.mdp',c=f'{gro}.gro',p=f'{top}.top',o=f'{out}.tpr',maxwarn=maxwarn)
    c.run(quiet=quiet)
    c=Command(f'{sw.mdrun}',deffnm=out,rdd=rdd,dds=dds,dlb=dlb,tunepme=tunepme)
    c.run(quiet=quiet,ignore_codes=ignore_codes)
    if os.path.exists(f'{out}.gro'):
        pass
        # logger.info(f'grompp_and_mdrun completed.  Check {gro}.gro.')
    else:
        logger.error(f'{sw.mdrun} ended prematurely; {gro}.gro not found.')
        raise Exception(f'{sw.mdrun} ended prematurely; {gro}.gro not found.')

def get_energy_menu(edr,**kwargs):
    """Get then menu provided by 'gmx energy' when applied to a particular edr file

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

def gmx_energy_trace(edr,names=[],**kwargs):
    """Generate traces of data in edr file

    :param edr: name of edr file
    :type edr: str
    :param names: list of data names, defaults to []
    :type names: list, optional
    :return: dataframe of traces
    :rtype: pandas DataFrame
    """
    assert os.path.exists(edr+'.edr'),f'Error: {edr}.edr not found'
    assert len(names)>0,f'Nothing to plot'
    xshift=kwargs.get('xshift',0)
    report_averages=kwargs.get('report_averages',False)
    menu=get_energy_menu(edr)
    with open('gmx.in','w') as f:
        for i in names:
            if i in menu:
                f.write(f'{menu[i]}\n')
            f.write('\n')
    if any([i in menu for i in names]):
        c=Command(f'{sw.gmx} {sw.gmx_options} energy -f {edr}.edr -o {edr}-out.xvg -xvg none < gmx.in')
        c.run()
        os.remove('gmx.in')
        colnames=names
        colnames.insert(0,'time (ps)')
        data=pd.read_csv(f'{edr}-out.xvg',sep='\s+',header=None,names=colnames)
        data.iloc[:,0]+=xshift
        ndata=data.shape[0]
        if report_averages:
            for i in names[1:]:
                data[f'Running-average-{i}']=data[i].expanding(1).mean()
                data[f'Rolling-10-average-{i}']=data[i].rolling(window=ndata//10).mean()
                for ln in data.iloc[-1][[i,f'Running-average-{i}',f'Rolling-10-average-{i}']].to_string().split('\n'):
                    logger.info(f'{ln}')
        return data
    else:
        return pd.DataFrame()
        
# def density_trace(edr,**kwargs):
#     """Report density trace from edr file

#     :param edr: edr file name
#     :type edr: str
#     """
#     msg=kwargs.get('msg','Density:')
#     with open('gmx.in','w') as f:
#         f.write('22\n\n')
#     c=Command(f'{sw.gmx} {sw.gmx_options} energy -f {edr}.edr -o {edr}-density.xvg -xvg none < gmx.in')
#     c.run()
#     os.remove('gmx.in')
#     density=pd.read_csv(f'{edr}-density.xvg',sep='\s+',names=['time(ps)','density(kg/m^3)'])
#     density['Running-average-density']=density['density(kg/m^3)'].expanding(1).mean()
#     density['Rolling-average-10']=density['density(kg/m^3)'].rolling(window=10).mean()
#     for ln in density.iloc[-1].to_string().split('\n'):
#         logger.info(f'{ln}')

def gromacs_distance(idf,gro,new_column_name='r',force_recalculate=False):
    """Use 'gmx distance' to measure interatomic distances

    :param idf: dataframe of atom indexes in pairs ['ai','aj']
    :type idf: pandas DataFrame
    :param gro: name of gromacs input file for 'gmx distance' to use
    :type gro: str
    :param new_column_name: name of column in idf where distances are stored, default 'r'
    :type new_column_name: str, optional
    :param force_recalculate: flag to force calculation of distances even if a distance column exists in idf, default False
    :type force_recalculate: boolean, optional
    :return: list of distances parallel to idf columns
    :rtype: numpy.ndarray
    """
    npair=idf.shape[0]
    # logger.debug(f'idf dtype {idf["ai"].dtype}')
    if 'r' in idf and not force_recalculate:
        return None
    ''' create the index file '''
    with open('tmp.ndx','w') as f:
        f.write('[ bonds-1 ]\n')
        idf[['ai','aj']].to_csv(f,sep=' ',header=False,index=False)
        # logger.debug('wrote tmp.ndx')
    ''' create the user-input file '''
    with open ('gmx.in','w') as f:
        f.write('0\n\n')
    c=Command(f'{sw.gmx} {sw.gmx_options} distance -f {gro} -n tmp.ndx -oall tmp.xvg -xvg none < gmx.in')
    c.run()
    ''' read the xvg file that 'gmx distance' creates '''
    with open ('tmp.xvg','r') as f:
        datastr=f.read().split()
    datastr=datastr[1:]
    # logger.debug(f'0 {datastr[0]} - {len(datastr)-1} {datastr[-1]}')
    data=np.array([float(x) for x in datastr])
    nd=len(data)
    ''' clean up '''
    os.remove('tmp.ndx')
    os.remove('gmx.in')
    os.remove('tmp.xvg')
    ''' add the lengths back to the input dataframe '''
    idf[new_column_name]=data
    assert npair==nd
    return data

def encluster(i,j,c):
    # NOT USED
    if c[i]==c[j]:
        return True
    if c[j]>c[i]:
        for k in range(len(c)):
            if c[k]==c[j]:
                c[k]=c[i]
    elif c[j]<c[i]:
        for k in range(len(c)):
            if c[k]==c[i]:
                c[k]=c[j]
    return False        

def symm(d,thresh=0.1,outfile=None):
    ''' Builds and returns an atom-idx-ordered list of sea-cluster indexes.
        Any two atoms with the same sea-cluster-index are considered
        symmetry equivalent.
        
        NOT USED

        Parameters:
           - d: interatomic distance matrix
           - thresh:  threshhold below which two sorted columns of d are considered
                      the "same" when compared to the magnitude of their straight
                      difference
    '''
    na=d.shape[0]
    # initialize per-atom cluster id's to atom indexes (0-)
    cluster_ids=np.arange(na)
    # sort every column of distance matrix
    d=np.sort(d,axis=0)
    # make a list of columns
    c=np.array([d[:,i] for i in range(d.shape[1])])
    l=[]
    # compute magnitude of distance between columns for 
    # all unique pairs of columns, store as a list
    # of 3-tuples (i,j,distance)
    for i,ri in enumerate(c):
        for j,rj in enumerate(c):
            if i<j: 
                ij=ri-rj
                l.append((i,j,np.sqrt(ij.dot(ij))))

    # sort the list by distances (not strictly necessary)
    L=sorted(l,key=lambda x: x[2])
    if outfile:
        with open(outfile,'w') as f:
            for i,j,r in L:
                f.write(f'{i} {j} {r:0.8f}\n')

    # clusterize the i,j entries of each element of L
    all_done=False
    cpass=0
    while not all_done:
        all_done=True
        for i,j,r in L:
            if r<thresh:
                all_done=encluster(i,j,cluster_ids)
        cpass+=1

    # reindex cluster id 
    cpop={}
    for i,c in enumerate(cluster_ids):
        if not c in cpop:
            cpop[c]=[]
        cpop[c].append(i)
    rc=set(cpop.keys())
    for i,r in enumerate(rc):
        for j in cpop[r]:
            cluster_ids[j]=i
    return cluster_ids

def analyze_sea(deffnm,thresh=0.1):
    ''' Builds and returns an atom-idx-ordered list of sea-cluster indexes.
        Any two atoms with the same sea-cluster-index are considered
        symmetry equivalent.
        
        NOT USED

        Parameters:
           - deffnm: ouput filename prefix for the Gromacs mdrun that generated the 
                     trajectory file we are analyzing here
           - thresh: threshhold below which two sorted columns of d are considered
                     the "same" when compared to the magnitude of their straight
                     difference

        The main job of this method is to compute the time-averaged interatomic
        distance matrix.  This matrix, if computed from a "hot" md simulation,
        should reveal atoms that are topologically symmetric, since the set of 
        average interatomic distances from atom A to all other atoms and the set
        of average interatomic distances from atom B to all other atoms are the 
        same if A and B are symmetry-equivalent.
    '''
    if not os.path.exists(f'{deffnm}.trr'):
        logger.error(f'{deffnm}.trr not found.')
        return []
    logger.debug(f'SEA analysis from {deffnm}.trr')
    with GroTrrReader(f'{deffnm}.trr') as trrfile:
        d=np.array((0,))
        nframes=0
        for mobydick in trrfile:
            frame_natoms=mobydick['natoms']
            if d.shape==(1,):
                d=np.zeros((frame_natoms,frame_natoms))
            data=trrfile.get_data()
            # tally all interatomic distances
            for i,ri in enumerate(data['x']):
                for j,rj in enumerate(data['x']):
                    rij=ri-rj
                    dist=np.sqrt(rij.dot(rij))
                    d[i][j]+=dist
            nframes+=1
        # averages over frames
        d/=nframes
        logger.debug(f'{deffnm}.trr: {nframes} frames')
        # send the distance matrix to be processed, return
        # the atom-ordered list of sea-cluster-idx's
        return symm(d,thresh=thresh,outfile=f'{deffnm}-symmanalysis.dat')

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
    