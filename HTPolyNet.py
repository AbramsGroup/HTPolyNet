#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""HTPolyNet is High-Throughput POLYmer NETwork utility to create atomistic models for MD Simulations 
created by the Abrams group at Drexel University
@author: ketankhare
@author: huang
"""

import sys
import os
import argparse


if __name__ == '__main__':

# Subsequent code changes to directory appropriate to the tools and invokes the main scripts

    parser = argparse.ArgumentParser(description='Select Simulator: gromacs (default) or lammps')
    parser.add_argument('--simulator', '-s', default="gromacs", choices=['gromacs', 'lammps'], help="Specify lammps, else default is gromacs")
    parser.add_argument('--config', '-c', default="createtemplate", help="Only for lammps. Default: Create template. Otherwise specify config file. Ignored for gromacs.")
    args = parser.parse_args()

    if(args.simulator=="gromacs"):

        print("Processing for GROMACS as simulator")

        try: 
            script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
            os.chdir(script_path + "/mh-cl") 
            #os.system("python3 main.py")
            sys.path.insert(0, os.path.abspath("../mh-cl"))
            import main as gromacs_module
            a = gromacs_module.main()
            a.preparePara()
            a.mainProcess(2)

        except:
            print("Error! Gromacs module calls have malfunctioned.")

    elif(args.simulator=="lammps"):

        print("Processing for LAMMPS as simulator")
        try: 
            module_script_abs_path = os.path.dirname(os.path.realpath(sys.argv[0])) + "/lammps-module/lammps_module_main.py"
            os.system("python3 " + module_script_abs_path + " --config " + args.config) 
        except:
            print("Error! LAMMPS module calls have malfunctioned!")
    else:

        print("Invalid entry! Use HTPolyNet.py -h for assistance.")
