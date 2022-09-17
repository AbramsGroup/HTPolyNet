from http.cookiejar import CookiePolicy
import logging
import os
import argparse as ap
import textwrap
import shutil
import networkx as nx
import numpy as np
from HTPolyNet.banner import banner, banner_message
from HTPolyNet.runtime import Runtime,logrotate
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.plot import diagnostics_graphs,global_trace,network_graph
from HTPolyNet.stringthings import my_logger
from HTPolyNet.utils import density_evolution, graph_from_bondsfile, init_molecule_graph,mwbxl,clusters
from HTPolyNet.inputcheck import input_check
from HTPolyNet.postprocess import postprocess, tladder

logger=logging.getLogger(__name__)

def info(args):
    print('This is some information on your installed version of HTPolyNet')
    l=pfs.lib_setup()
    software.sw_setup()
    print(l.info())
    print(software.to_string())

def run(args):

    logrotate(args.diag)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.diag,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    if not args.no_banner: banner(logger.info)
    my_logger('HTPolyNet runtime begins',logger.info)
    userlib=args.lib if os.path.exists(args.lib) else None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.do_workflow(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolyNet runtime ends',logger.info)

def parameterize(args):

    logrotate(args.diag)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.diag,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    if not args.no_banner: banner(logger.info)
    # banner(logger.info)
    my_logger('HTPolyNet parameterization begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.generate_molecules(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolynet parameterization ends',logger.info)

def htpolynet_cure_plots(args):
    logs=args.logs
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(format='%(levelname)s> %(message)s',level=loglevel_numeric)

    if not args.no_banner: banner(logger.info)
    # banner(print)
    if len(logs)>0:
        diagnostics_graphs(logs,args.plotfile)
    if args.proj:
        if args.t:
            df,transition_times,cure_markers,interval_labels=density_evolution(args.proj)
            global_trace(df,['Temperature','Density','Potential'],args.t,transition_times=transition_times,markers=[],interval_labels=interval_labels,y2names=['nbonds','nbonds'],legend=True)
            if args.o:
                logger.info(f'All data to {args.o}')
                with open(args.o,'w') as f:
                    f.write(df.to_string(index=False,float_format='{:.3f}'.format)+'\n')
        if args.g or args.mwbxl or args.clusters:
            G=init_molecule_graph(args.proj)
            n=1
            while os.path.exists(os.path.join(args.proj,f'systems/iter-{n}/2-cure_update-bonds.csv')):
                logger.info(f'iter-{n}/2-cure_update-bonds.csv')
                g=graph_from_bondsfile(os.path.join(args.proj,f'systems/iter-{n}/2-cure_update-bonds.csv'))
                G=nx.compose(G,g)
                if args.g and args.byiter: 
                    network_graph(G,os.path.join(args.proj,f'plots/iter-{n}-{args.g}'))
                n+=1
            if args.g: network_graph(G,args.g)
            # print('Covalent cluster sizes (in monomers):')
            if args.clusters:
                clu=clusters(G)
                clu.to_csv(args.clusters,sep=' ',header=True,index=False)
                print(f'{args.clusters} created.')
            if args.mwbxl:
                am=mwbxl(G)
                print(f'Avg homo-N between xlinks: {np.average(am["n"],weights=am["counts"]):.2f}')
                am.to_csv(args.mwbxl,sep=' ',index=False,header=True)
                print(f'{args.mwbxl} created.')

def fetch_example(args):
    l=pfs.system()
    nm=args.n
    kp=args.k
    possibles=l.get_example_names()
    if nm.isdecimal():
        fullname=[x for x in possibles if x.startswith(nm)][0]
    else:
        if nm=='all':
            for fullname in possibles:
                shutil.copy(os.path.join(l.get_example_depot_location(),f'{fullname}.tgz'),'.')
                os.system(f'tar zxf {fullname}.tgz')
                if not kp:
                    os.remove(f'{fullname}.tgz')
            return
        fullname=nm
    shutil.copy(os.path.join(l.get_example_depot_location(),f'{fullname}.tgz'),'.')
    os.system(f'tar zxf {fullname}.tgz')
    if not kp:
        os.remove(f'{fullname}.tgz')


def cli():
    """cli Command-line interface
    """

    l=pfs.lib_setup()
    example_names=l.get_example_names()
    example_ids=[x.split('-')[0] for x in example_names]
    commands={}
    commands['run']=run
    commands['parameterize']=parameterize
    commands['info']=info
    commands['plots']=htpolynet_cure_plots
    commands['fetch-example']=fetch_example
    commands['input-check']=input_check
    commands['postprocess']=postprocess

    helps={}
    helps['run']='build a system using instructions in the config file and any required molecular structure inputs'
    helps['parameterize']='parameterize monomers and oligomer templates using instructions in the config file'
    helps['info']='print some information to the console'
    helps['plots']='generate some plots that summarize aspects of the current completed build'
    helps['fetch-example']='fetch and unpack example(s) from the HTPolyNet.Library: '+', '.join([f'"{x}"' for x in l.get_example_names()])
    helps['input-check']='reports number of atoms that would be in initial system based on config'
    helps['postprocess']='perform specified postprocessing calculations on final results in one or more project directory'

    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,help=helps[k])
        command_parsers[k].set_defaults(func=commands[k])

    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['run'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory.  Anything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['run'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['run'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['run'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['run'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['run'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['run'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['parameterize'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['parameterize'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['parameterize'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory.  Anything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['parameterize'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['parameterize'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['parameterize'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['parameterize'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['parameterize'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['parameterize'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['plots'].add_argument('-logs',type=str,default='',nargs='+',help='names of diagnostic log files (1 or more)')
    command_parsers['plots'].add_argument('-proj',type=str,default='',help='name of project directory')
    command_parsers['plots'].add_argument('-t',type=str,default='',help='Plot density and temperature traces for entire build in specified project directory to this file')
    command_parsers['plots'].add_argument('-o',type=str,default='',help='dump density/temperature trace data to this file')
    command_parsers['plots'].add_argument('-g',type=str,default='',help='Plot graph network of resids and save to this file name')
    command_parsers['plots'].add_argument('-byiter',default=False,action='store_true',help='Plot graph network of resids for each iter separately')
    command_parsers['plots'].add_argument('-mwbxl',default='',help='Compute home-N between crosslinks, save to this file')
    command_parsers['plots'].add_argument('-clusters',default='',help='Save cluster size distribution to this file')
    command_parsers['plots'].add_argument('--plotfile',type=str,default='cure-info.png',help='name of plot file to generate')
    command_parsers['plots'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['plots'].add_argument('--loglevel',type=str,default='info',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['fetch-example'].add_argument('-n',type=str,choices=example_ids+['all'],help='number of example tarball to unpack from '+', '.join(example_names))
    command_parsers['fetch-example'].add_argument('-k',default=False,action='store_true',help='keep tarballs')

    command_parsers['input-check'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['input-check'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')

    command_parsers['postprocess'].add_argument('-proj',type=str,default='',nargs='+',help='name of project directory')
    command_parsers['postprocess'].add_argument('-Tladder',type=tladder, default='',help='run a temperature-ladder simulation to measure density=f(T); format (T0,T1,Ntemps,ps_per_run,ps_per_rise,ps_warmup')
    command_parsers['postprocess'].add_argument('-cfg',type=str,default='',help='config file used to generate this project')
    command_parsers['postprocess'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['postprocess'].add_argument('--loglevel',type=str,default='info',help='Log level for messages written to diagnostic log (debug|info)')

    args=parser.parse_args()
    args.func(args)
