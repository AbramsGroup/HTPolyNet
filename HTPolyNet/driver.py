"""

.. module:: driver
   :synopsis: manages the HTPolyNet application, provides the command-line interface entry point
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import logging
import os
import argparse as ap
import textwrap
import shutil
import glob
from .banner import banner, banner_message
from .runtime import Runtime,logrotate
from .command import Command
import HTPolyNet.projectfilesystem as pfs
import HTPolyNet.software as software
from .plot import plots
from .stringthings import my_logger
from .inputcheck import input_check
from .postsim import postsim
from .analyze import analyze

logger=logging.getLogger(__name__)

def info(args):
    """info handles the info subcommmand

    :param args: parsed arguments
    :type args: argparse.Namespace
    """
    print('This is some information on your installed version of HTPolyNet')
    l=pfs.lib_setup()
    software.sw_setup()
    print(l.info())
    print(software.to_string())
    possibles=l.get_example_names()
    print('Available examples using htpolynet fetch-example')
    for p in possibles:
        if p[0].isdecimal():
            print(f'   {p}')

def run(args):
    """run handles the run subcommand

    :param args: parsed arguments
    :type args: argparse.Namespace
    """
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
    """parameterize handles the parameterize subcommand

    :param args: parsed arguments
    :type args: argparse.Namespace
    """
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
    my_logger('HTPolyNet parameterization begins',logger.info)
    userlib=args.lib
    if not os.path.exists(args.lib):
        userlib=None
    software.sw_setup()
    pfs.pfs_setup(root=os.getcwd(),topdirs=['molecules','systems','plots'],verbose=True,projdir=args.proj,reProject=args.restart,userlibrary=userlib)
    a=Runtime(cfgfile=args.config,restart=args.restart)
    a.generate_molecules(force_checkin=args.force_checkin,force_parameterization=args.force_parameterization)
    my_logger('HTPolynet parameterization ends',logger.info)

def fetch_example(args):
    """fetch_example handles the fetch-example subcommand

    :param args: parsed arguments
    :type args: argparse.Namespace
    """
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
    fullpath=os.path.join(l.get_example_depot_location(),f'{fullname}.tgz')
    assert os.path.exists(fullpath),f'{fullpath}: not found.'
    shutil.copy(fullpath,'.')
    os.system(f'tar zxf {fullname}.tgz')
    if not kp:
        os.remove(f'{fullname}.tgz')

def pack_example(args):
    l=pfs.system()
    requires=['README.md','run.sh','lib/molecules']
    for r in requires:
        assert os.path.exists(r),f'{r}: not found.'
    ls=glob.glob('*.yaml')+glob.glob('*.yml')
    assert len(ls)>0,f'No yaml config file found.'
    bn=os.path.basename(os.getcwd())
    inspectname=bn.split('-')
    if len(inspectname)>0:
        firstfield=inspectname[0]
        if firstfield.isdigit():
            existing_n=int(firstfield)
            use_n=False
        else:
            use_n=True
    n=args.n
    existing_examples=l.get_example_names()
    numbers=[int(x.split('-')[0]) for x in existing_examples]
    assert not n in numbers,f'Choose an index other than {n}; an example with this index already exists.'
    if not use_n:
        assert not existing_n in numbers,f'Choose an index other than {n}; an example with this index already exists. To do this, rename the directory.'
    if n==-1:
        n=str(max(numbers)+1)
    if use_n:
        newname=f'{n}-{os.path.basename(os.getcwd())}'
    else:
        newname=f'{os.path.basename(os.getcwd())}'
    depot_location=l.get_example_depot_location()
    os.chdir('..')
    c=Command(f'tar --exclude="*/*/*/*/*" -zvcf {depot_location}/{newname}.tgz {bn}/README.md {bn}/run.sh {bn}/*.yaml {bn}/lib/molecules/')
    o,e=c.run()
    print(o)
    print(e)

def cli():
    """cli Command-line interface
    """
    l=pfs.lib_setup()
    example_names=l.get_example_names()
    example_ids=[x.split('-')[0] for x in example_names]
    # declare subcommands and their handlers
    commands={}
    commands['run']=run
    commands['parameterize']=parameterize
    commands['info']=info
    commands['plots']=plots
    commands['fetch-example']=fetch_example
    commands['input-check']=input_check
    commands['postsim']=postsim
    commands['analyze']=analyze
    commands['pack-example']=pack_example
    # declare subcommand helpstrings
    helps={}
    helps['run']='build a system using instructions in the config file and any required molecular structure inputs'
    helps['parameterize']='parameterize monomers and oligomer templates using instructions in the config file'
    helps['info']='print some information to the console'
    helps['plots']='generate some plots that summarize aspects of the current completed build'
    helps['fetch-example']='fetch and unpack example(s) from HTPolyNet.HTPolyNet.resources: '+', '.join([f'"{x}"' for x in l.get_example_names()])
    helps['pack-example']='pack current directory as a tarball to the resources/example_depot'
    helps['input-check']='reports number of atoms that would be in initial system based on config'
    helps['postsim']='perform specified postcure MD simulations on final results in one or more project directory'
    helps['analyze']='perform \'gmx <command>\' style analyses specified in the cfg file'
    parser=ap.ArgumentParser(description=textwrap.dedent(banner_message),formatter_class=ap.RawDescriptionHelpFormatter)
    # Subparsers, one per subcommand
    subparsers=parser.add_subparsers()
    subparsers.required=True
    command_parsers={}
    for k in commands:
        command_parsers[k]=subparsers.add_parser(k,help=helps[k])
        command_parsers[k].set_defaults(func=commands[k])
    ######## run ########
    command_parsers['run'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['run'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['run'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory.  Anything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['run'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['run'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['run'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['run'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['run'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['run'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
    ######## parameterize ########
    command_parsers['parameterize'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['parameterize'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['parameterize'].add_argument('-proj',type=str,default='next',help='project directory; "next" (default) generates next directory.  Anything other than "next": if it exists, "-restart" must be included as a parameter; if not, it is created as a new project')
    command_parsers['parameterize'].add_argument('-diag',type=str,default='htpolynet_runtime_diagnostics.log',help='diagnostic log file')
    command_parsers['parameterize'].add_argument('-restart',default=False,action='store_true',help='restart in latest proj dir')
    command_parsers['parameterize'].add_argument('--force-parameterization',default=False,action='store_true',help='force GAFF parameterization of any input mol2 structures')
    command_parsers['parameterize'].add_argument('--force-checkin',default=False,action='store_true',help='force check-in of any generated parameter files to the system library')
    command_parsers['parameterize'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['parameterize'].add_argument('--loglevel',type=str,default='debug',help='Log level for messages written to diagnostic log (debug|info)')
    ######## plots ########
    command_parsers['plots'].add_argument('source',type=str,choices=['diag','build','post'],default='build',help='source of data to plot')
    command_parsers['plots'].add_argument('--diags',type=str,default=[],nargs='+',help='names of diagnostic log files (1 or more)')
    command_parsers['plots'].add_argument('--proj',nargs='+',type=str,default=[],help='name of project director[y/ies]')
    command_parsers['plots'].add_argument('--cfg',type=str,nargs='+',default=[],help='name input config files')
    command_parsers['plots'].add_argument('--buildplot',type=str,nargs='+',default=['t'],choices=['t','g','n','c'],help='type of build plot to generate: t: traces (select using --traces); g: 2-D graph representations iteration by iteration; n: homo-N between crosslinks; c: cluster-size distributions')
    command_parsers['plots'].add_argument('--traces',type=str,nargs='+',default=['t','d','p'],choices=['t','d','p'],help='type of traces to plot from build: t: temperature; d: density; p: potential energy')
    command_parsers['plots'].add_argument('--n_points',type=int,nargs=2,default=[10,20],help='number of [cold-side,hot-side] data points in the Tg analysis to fit lines to')

    command_parsers['plots'].add_argument('--plotfile',type=str,default='',help='name of plot file to generate')
    command_parsers['plots'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['plots'].add_argument('--loglevel',type=str,default='info',help='Log level for messages written to diagnostic log (debug|info)')
    ######## fetch-example ########
    command_parsers['fetch-example'].add_argument('-n',type=str,choices=example_ids+['all'],help='number of example tarball to unpack from '+', '.join(example_names))
    command_parsers['fetch-example'].add_argument('-k',default=False,action='store_true',help='keep tarballs')
    ######## pack-example
    command_parsers['pack-example'].add_argument('-n',type=int,default=-1,help='desired index (integer) for this example')
    ######## input-check ########
    command_parsers['input-check'].add_argument('config',type=str,default=None,help='input configuration file in YAML format')
    command_parsers['input-check'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    ######## postsim ########
    command_parsers['postsim'].add_argument('-proj',type=str,default='',nargs='+',help='name of project directory')
    command_parsers['postsim'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['postsim'].add_argument('-ocfg',type=str,default='',help='original HTPolyNet config file used to generate project(s)')
    command_parsers['postsim'].add_argument('-cfg',type=str,default='',help='config file for specifying the MD simulations to perform')
    command_parsers['postsim'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['postsim'].add_argument('--loglevel',type=str,default='info',help='Log level for messages written to diagnostic log (debug|info)')

    ######## analyze ########
    command_parsers['analyze'].add_argument('-proj',type=str,default='',nargs='+',help='name of project directory')
    command_parsers['analyze'].add_argument('-lib',type=str,default='lib',help='local user library of molecular structures and parameterizations')
    command_parsers['analyze'].add_argument('-ocfg',type=str,default='',help='original HTPolyNet config file used to generate project(s)')
    command_parsers['analyze'].add_argument('-cfg',type=str,default='',help='config file for specifying the analyses to perform')
    command_parsers['analyze'].add_argument('--no-banner',default=False,action='store_true',help='turn off the banner')
    command_parsers['analyze'].add_argument('--loglevel',type=str,default='info',help='Log level for messages written to diagnostic log (debug|info)')
    args=parser.parse_args()
    args.func(args)
