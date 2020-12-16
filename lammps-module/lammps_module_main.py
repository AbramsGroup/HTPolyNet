#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Main file for the LAMMPS module of HTPolyNet utility

@author: KetanSKhare"""

import sys
import os
import argparse


'''Create template'''

def lammps_module_main():

    script_path = os.path.dirname(os.path.realpath(sys.argv[0]))
    os.chdir(script_path + "/../mh-cl") 
    sys.path.insert(0, os.path.abspath("../mh-cl"))
    
    import main as gromacs_module
    a = gromacs_module.main()
    a.preparePara("lammps")

if __name__ == '__main__':
    lammps_module_main()
