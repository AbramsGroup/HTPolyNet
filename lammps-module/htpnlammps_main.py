#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Main file for the LAMMPS module of HTPolyNet utility"""

import sys
import os
import argparse
import ConfigParser


'''Create template'''

if __name__ == '__main__':
    
    #Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--config')
    args=parser.parse_args()
    
    