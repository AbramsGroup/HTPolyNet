"""

.. module:: test_gromacs_get_energy_menu
   :synopsis: tests HTPolyNet.gromacs.gmx_energy_trace() function; only checks congruency of data types
   
.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import importlib.resources
import logging
logger=logging.getLogger(__name__)
from numpy import dtype
import os

from HTPolyNet.gromacs import gmx_energy_trace

class TestGmxEnergyTrace(unittest.TestCase):
    def setUp(self,loglevel='debug',logfile='testlog.log'):
        log_destination=str(importlib.resources.files('tests').joinpath(logfile))
        loglevel_numeric=getattr(logging,loglevel.upper())
        logging.basicConfig(filename=log_destination,filemode='a',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
        self.edr_extract_results={'items45':{'time(ps)': dtype('float64'), 'Bond': dtype('float64'), 'Angle': dtype('float64'), 'Proper-Dih.': dtype('float64'), 'Periodic-Improper-Dih.': dtype('float64'), 'LJ-14': dtype('float64'), 'Coulomb-14': dtype('float64'), 'LJ-(SR)': dtype('float64'), 'Disper.-corr.': dtype('float64'), 'Coulomb-(SR)': dtype('float64'), 'Coul.-recip.': dtype('float64'), 'Potential': dtype('float64'), 'Kinetic-En.': dtype('float64'), 'Total-Energy': dtype('float64'), 'Conserved-En.': dtype('float64'), 'Temperature': dtype('float64'), 'Pres.-DC': dtype('float64'), 'Pressure': dtype('float64'), 'Box-X': dtype('float64'), 'Box-Y': dtype('float64'), 'Box-Z': dtype('float64'), 'Volume': dtype('float64'), 'Density': dtype('float64'), 'pV': dtype('float64'), 'Enthalpy': dtype('float64'), 'Vir-XX': dtype('float64'), 'Vir-XY': dtype('float64'), 'Vir-XZ': dtype('float64'), 'Vir-YX': dtype('float64'), 'Vir-YY': dtype('float64'), 'Vir-YZ': dtype('float64'), 'Vir-ZX': dtype('float64'), 'Vir-ZY': dtype('float64'), 'Vir-ZZ': dtype('float64'), 'Pres-XX': dtype('float64'), 'Pres-XY': dtype('float64'), 'Pres-XZ': dtype('float64'), 'Pres-YX': dtype('float64'), 'Pres-YY': dtype('float64'), 'Pres-YZ': dtype('float64'), 'Pres-ZX': dtype('float64'), 'Pres-ZY': dtype('float64'), 'Pres-ZZ': dtype('float64'), '#Surf*SurfTen': dtype('float64'), 'T-System': dtype('float64'), 'Lamb-System': dtype('float64')},
        'items43':{'time(ps)': dtype('float64'), 'Angle': dtype('float64'), 'Proper-Dih.': dtype('float64'), 'Periodic-Improper-Dih.': dtype('float64'), 'LJ-14': dtype('float64'), 'Coulomb-14': dtype('float64'), 'LJ-(SR)': dtype('float64'), 'Coulomb-(SR)': dtype('float64'), 'Coul.-recip.': dtype('float64'), 'Potential': dtype('float64'), 'Kinetic-En.': dtype('float64'), 'Total-Energy': dtype('float64'), 'Conserved-En.': dtype('float64'), 'Temperature': dtype('float64'), 'Pressure': dtype('float64'), 'Constr.-rmsd': dtype('float64'), 'Box-X': dtype('float64'), 'Box-Y': dtype('float64'), 'Box-Z': dtype('float64'), 'Volume': dtype('float64'), 'Density': dtype('float64'), 'pV': dtype('float64'), 'Enthalpy': dtype('float64'), 'Vir-XX': dtype('float64'), 'Vir-XY': dtype('float64'), 'Vir-XZ': dtype('float64'), 'Vir-YX': dtype('float64'), 'Vir-YY': dtype('float64'), 'Vir-YZ': dtype('float64'), 'Vir-ZX': dtype('float64'), 'Vir-ZY': dtype('float64'), 'Vir-ZZ': dtype('float64'), 'Pres-XX': dtype('float64'), 'Pres-XY': dtype('float64'), 'Pres-XZ': dtype('float64'), 'Pres-YX': dtype('float64'), 'Pres-YY': dtype('float64'), 'Pres-YZ': dtype('float64'), 'Pres-ZX': dtype('float64'), 'Pres-ZY': dtype('float64'), 'Pres-ZZ': dtype('float64'), '#Surf*SurfTen': dtype('float64'), 'T-System': dtype('float64'), 'Lamb-System': dtype('float64')},
        'items31':{'time(ps)': dtype('float64'), 'Bond': dtype('float64'), 'Angle': dtype('float64'), 'Proper-Dih.': dtype('float64'), 'Periodic-Improper-Dih.': dtype('float64'), 'LJ-14': dtype('float64'), 'Coulomb-14': dtype('float64'), 'LJ-(SR)': dtype('float64'), 'Coulomb-(SR)': dtype('float64'), 'Coul.-recip.': dtype('float64'), 'Potential': dtype('float64'), 'Pressure': dtype('float64'), 'Vir-XX': dtype('float64'), 'Vir-XY': dtype('float64'), 'Vir-XZ': dtype('float64'), 'Vir-YX': dtype('float64'), 'Vir-YY': dtype('float64'), 'Vir-YZ': dtype('float64'), 'Vir-ZX': dtype('float64'), 'Vir-ZY': dtype('float64'), 'Vir-ZZ': dtype('float64'), 'Pres-XX': dtype('float64'), 'Pres-XY': dtype('float64'), 'Pres-XZ': dtype('float64'), 'Pres-YX': dtype('float64'), 'Pres-YY': dtype('float64'), 'Pres-YZ': dtype('float64'), 'Pres-ZX': dtype('float64'), 'Pres-ZY': dtype('float64'), 'Pres-ZZ': dtype('float64'), '#Surf*SurfTen': dtype('float64'), 'T-rest': dtype('float64')}}

    def test_edr(self):
        """
        Test of HTPolyNet.gromacs.gmx_energy_trace on three different edr files
        """
        for i,o in self.edr_extract_results.items():
            names=list(o.keys())
            logger.info(f'{i}')
            test_file_path_str = str(importlib.resources.files('tests.unit').joinpath(f'fixtures/{i}'))
            test_dataframe = gmx_energy_trace(test_file_path_str,names)
            res=test_dataframe.dtypes.to_dict()
            self.assertEqual(res,o)
    
    def test_avg(self):
        """
        Test of the report_averages keyword for HTPolyNet.gromacs.gmx_energy_trace
        """
        for i,o in self.edr_extract_results.items():
            names=list(o.keys())
            logger.info(f'{i}')
            test_file_path_str = str(importlib.resources.files('tests.unit').joinpath(f'fixtures/{i}'))
            test_dataframe = gmx_energy_trace(test_file_path_str,names,report_averages=True)
            all_there=[]
            for n in names[1:]:
                all_there.append(f'Running-average-{n}' in test_dataframe.columns)
                all_there.append(f'Rolling-10-average-{n}' in test_dataframe.columns)
            self.assertEqual(all(all_there),True)
    
    def test_keepfiles(self):
        """
        Test of the keep_files keyword for HTPolyNet.gromacs.gmx_energy_trace
        """
        for i,o in self.edr_extract_results.items():
            names=list(o.keys())
            logger.info(f'{i}')
            test_file_path_str = str(importlib.resources.files('tests.unit').joinpath(f'fixtures/{i}'))
            test_dataframe = gmx_energy_trace(test_file_path_str,names,keep_files=True)
            tmpfiles=[f'{i}-gmx.in',f'{i}-out.xvg']
            allthere=[]
            for t in tmpfiles:
                fn=str(importlib.resources.files('tests.unit').joinpath(f'fixtures/{t}'))
                if os.path.exists(fn):
                    allthere.append(True)
                    os.remove(fn)
            self.assertEqual(all(allthere),True)
