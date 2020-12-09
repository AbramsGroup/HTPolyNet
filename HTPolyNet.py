#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""HTPolyNet is High-Throughput POLYmer NETwork utility to create atomistic models for MD Simulations 
created by the Abrams group at Drexel University"""

import sys
import os
import argparse

# Subsequent code changes to directory appropriate to the tools and invokes the main scripts

parser = argparse.ArgumentParser(description='Select Simulator: gromacs (default) or lammps')
parser.add_argument('--simulator', '-s', default="gromacs", help="Specify lammps, else default is gromacs")
parser.add_argument('--config', '-c', default="createtemplate", help="Only for lammps. Default: Create template. Otherwise specify config file. Ignored for gromacs.")
args = parser.parse_args()

script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
current_workdir = os.getcwd()

if(args.simulator=="gromacs"):

    print("Processing for GROMACS as simulator")

    try: 
        os.chdir(script_path + "/mh-cl") 
        os.system("python3 main.py")

    except:

        print("Error! Expected module files not found.")

elif(args.simulator=="lammps"):

    print("Processing for LAMMPS as simulator")
    try: 
        module_script_abs_path = os.path.dirname(os.path.realpath(sys.argv[0])) + "/lammps-module/htpnlammps_main.py"
s       os.system("python3 " + module_script_abs_path + " --config " + args.config) 
        
    except:
        print("Error! Expected module files not found.")
        
else:

    print("Invalid entry! Use HTPolyNet.py -h for assistance.")