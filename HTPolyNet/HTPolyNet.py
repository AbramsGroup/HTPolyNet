import logging
import os
import argparse as ap
import textwrap
import shutil

from HTPolyNet.banner import banner, banner_message
from HTPolyNet.runtime import Runtime,logrotate
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.plot import diagnostics_graphs,global_trace
from HTPolyNet.stringthings import my_logger
from HTPolyNet.utils import density_evolution

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

    banner(logger.info)
    my_logger('HTPolyNet runtime begins',logger.info)
    userlib=args.lib if os.path.exists(args.lib) else None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.build(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
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

    banner(logger.info)
    my_logger('HTPolyNet parameterization begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.generate_molecules(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolynet parameterization ends',logger.info)

def htpolynet_cure_plots(args):
    logs=args.logs
    banner(print)
    if len(logs)>0:
        diagnostics_graphs(logs,args.plotfile)
    if args.proj:
        df,transition_times,markers,interval_labels=density_evolution(args.proj)
        global_trace(df,['Temperature','Density'],'global_traces.png',transition_times=transition_times,markers=markers,interval_labels=interval_labels,y2names=['nbonds','nbonds'],legend=True)
        if args.o:
            print(f'All data to {args.o}')
            with open(args.o,'w') as f:
                f.write(df.to_string(header=['time(ps)','Temperature','nbonds','Density'],index=False,float_format='{:.3f}'.format)+'\n')

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

    helps={}
    helps['run']='build a system using instructions in the config file and any required molecular structure inputs'
    helps['parameterize']='parameterize monomers and oligomer templates using instructinos in the config file'
    helps['info']='print some information to the console'
    helps['plots']='generate some plots that summarize aspects of the current completed build'
    helps['fetch-example']='fetch and unpack example(s) from the HTPolyNet.Library: '+', '.join([f'"{x}"' for x in l.get_example_names()])

    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    subparsers=parser.add_subparsers()
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,help=helps[k])
        command_parsers[k].set_defaults(func=commands[k])

    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['run'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory\nAnything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['run'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['run'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['run'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['run'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['run'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['parameterize'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['parameterize'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['parameterize'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['parameterize'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['parameterize'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['parameterize'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['parameterize'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')

    command_parsers['plots'].add_argument('-logs',type=str,default='',nargs='+',help='names of diagnostic log files (1 or more)')
    command_parsers['plots'].add_argument('-proj',type=str,default='',help='name of project directory')
    command_parsers['plots'].add_argument('-o',type=str,default='',help='name of global trace output data file')
    command_parsers['plots'].add_argument('--plotfile',type=str,default='cure-info.png',help='name of plot file to generate')

    command_parsers['fetch-example'].add_argument('-n',type=str,choices=example_ids+['all'],help='number of example tarball to unpack from '+', '.join(example_names))
    command_parsers['fetch-example'].add_argument('-k',default=False,action='store_true',help='keep tarballs')

    args=parser.parse_args()
    args.func(args)
