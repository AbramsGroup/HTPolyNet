#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""HTPolyNet is High-Throughput POLYmer NETwork utility to create atomistic models for MD Simulations 
created by the Abrams group at Drexel University"""

import sys
import os
import argparse

# Subsequent code changes to directory appropriate to the tools and invokes the main scripts

parser = argparse.ArgumentParser(description='Select Simulator: gromacs (default) or lammps')
<<<<<<< HEAD
parser.add_argument('--simulator', '-s', default="gromacs", choices=['gromacs', 'lammps'], help="Specify lammps, else default is gromacs")
=======
parser.add_argument('--simulator', '-s', default="gromacs", help="Specify lammps, else default is gromacs")
>>>>>>> c6518ed5524a7269bc56af54fc99b227070888fc
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
<<<<<<< HEAD
        os.system("python3 " + module_script_abs_path + " --config " + args.config) 
=======
s       os.system("python3 " + module_script_abs_path + " --config " + args.config) 
>>>>>>> c6518ed5524a7269bc56af54fc99b227070888fc
        
    except:
        print("Error! Expected module files not found.")
        
else:

<<<<<<< HEAD
    print("Invalid entry! Use HTPolyNet.py -h for assistance.")
=======
    print("Invalid entry! Use HTPolyNet.py -h for assistance.")
>>>>>>> c6518ed5524a7269bc56af54fc99b227070888fc
