import logging
import os
import sys
import argparse as ap

from HTPolyNet.banner import banner
from HTPolyNet.runtime import Runtime,logrotate
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from HTPolyNet.plot import cure_graph,density_evolution
from HTPolyNet.stringthings import my_logger

logger=logging.getLogger(__name__)
parser=ap.ArgumentParser()

def info():
    print('This is some information on your installed version of HTPolyNet')
    l=pfs.lib_setup()
    print(l.info())
    software.info()

def run():
    parser.add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    parser.add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    parser.add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory\nAnything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    parser.add_argument('-log',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    parser.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    parser.add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    parser.add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    parser.add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
    args=parser.parse_args(sys.argv[2:])

    logrotate(args.log)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
    console=logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter=logging.Formatter('%(levelname)s> %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)

    banner(logger.info)
    my_logger('HTPolyNet runtime begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.build(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolyNet runtime ends',logger.info)

def parameterize():
    parser.add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    parser.add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    parser.add_argument('-log',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    parser.add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    parser.add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    parser.add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    parser.add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
    args=parser.parse_args(sys.argv[2:])

    logrotate(args.log)

    ''' set up logger with all debug+-level messages going to the diagnostic log file and info to console '''
    loglevel_numeric=getattr(logging, args.loglevel.upper())
    logging.basicConfig(filename=args.log,filemode='w',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
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

def htpolynet_cure_plots():
    parser.add_argument('logs',type=str,default='',nargs='+',help='names of diagnostic log files')
    parser.add_argument('--plotfile',type=str,default='cure-info.png',help='name of plot file to generate')
    args=parser.parse_args(sys.argv[2:])
    logs=args.logs
    banner(print)
    cure_graph(logs,args.plotfile)
    density_evolution()

def cli():
    """cli Command-line interface
    """
    commands={}
    commands['run']=run
    commands['parameterize']=parameterize
    commands['info']=info
    commands['plots']=htpolynet_cure_plots

    parser=ap.ArgumentParser()
    parser.add_argument('command',type=str,default='',help='command ('+', '.join(list(commands.keys()))+')')
    args=parser.parse_args(sys.argv[1:2])

    # TODO: use subparsers!

    if args.command in commands:
        commands[args.command]()
    else:
        print(f'HTPolyNet command {args.command} not recognized')

