import unittest
import importlib.resources
import logging
logger=logging.getLogger(__name__)

from HTPolyNet.gromacs import get_energy_menu

class TestGetEnergyMenu(unittest.TestCase):
    def setUp(self,loglevel='debug',logfile='testlog.log'):
        log_destination=str(importlib.resources.files('tests').joinpath(logfile))
        loglevel_numeric=getattr(logging,loglevel.upper())
        logging.basicConfig(filename=log_destination,filemode='a',format='%(asctime)s %(name)s.%(funcName)s %(levelname)s> %(message)s',level=loglevel_numeric)
        self.edr_menu_results={'items45':{"Bond": "1", "Angle": "2", "Proper-Dih.": "3", "Periodic-Improper-Dih.": "4", "LJ-14": "5", "Coulomb-14": "6", "LJ-(SR)": "7", "Disper.-corr.": "8", "Coulomb-(SR)": "9", "Coul.-recip.": "10", "Potential": "11", "Kinetic-En.": "12", "Total-Energy": "13", "Conserved-En.": "14", "Temperature": "15", "Pres.-DC": "16", "Pressure": "17", "Box-X": "18", "Box-Y": "19", "Box-Z": "20", "Volume": "21", "Density": "22", "pV": "23", "Enthalpy": "24", "Vir-XX": "25", "Vir-XY": "26", "Vir-XZ": "27", "Vir-YX": "28", "Vir-YY": "29", "Vir-YZ": "30", "Vir-ZX": "31", "Vir-ZY": "32", "Vir-ZZ": "33", "Pres-XX": "34", "Pres-XY": "35", "Pres-XZ": "36", "Pres-YX": "37", "Pres-YY": "38", "Pres-YZ": "39", "Pres-ZX": "40", "Pres-ZY": "41", "Pres-ZZ": "42", "#Surf*SurfTen": "43", "T-System": "44", "Lamb-System": "45"},
        'items43':{"Angle": "1", "Proper-Dih.": "2", "Periodic-Improper-Dih.": "3", "LJ-14": "4", "Coulomb-14": "5", "LJ-(SR)": "6", "Coulomb-(SR)": "7", "Coul.-recip.": "8", "Potential": "9", "Kinetic-En.": "10", "Total-Energy": "11", "Conserved-En.": "12", "Temperature": "13", "Pressure": "14", "Constr.-rmsd": "15", "Box-X": "16", "Box-Y": "17", "Box-Z": "18", "Volume": "19", "Density": "20", "pV": "21", "Enthalpy": "22", "Vir-XX": "23", "Vir-XY": "24", "Vir-XZ": "25", "Vir-YX": "26", "Vir-YY": "27", "Vir-YZ": "28", "Vir-ZX": "29", "Vir-ZY": "30", "Vir-ZZ": "31", "Pres-XX": "32", "Pres-XY": "33", "Pres-XZ": "34", "Pres-YX": "35", "Pres-YY": "36", "Pres-YZ": "37", "Pres-ZX": "38", "Pres-ZY": "39", "Pres-ZZ": "40", "#Surf*SurfTen": "41", "T-System": "42", "Lamb-System": "43"},
        'items31':{"Bond": "1", "Angle": "2", "Proper-Dih.": "3", "Periodic-Improper-Dih.": "4", "LJ-14": "5", "Coulomb-14": "6", "LJ-(SR)": "7", "Coulomb-(SR)": "8", "Coul.-recip.": "9", "Potential": "10", "Pressure": "11", "Vir-XX": "12", "Vir-XY": "13", "Vir-XZ": "14", "Vir-YX": "15", "Vir-YY": "16", "Vir-YZ": "17", "Vir-ZX": "18", "Vir-ZY": "19", "Vir-ZZ": "20", "Pres-XX": "21", "Pres-XY": "22", "Pres-XZ": "23", "Pres-YX": "24", "Pres-YY": "25", "Pres-YZ": "26", "Pres-ZX": "27", "Pres-ZY": "28", "Pres-ZZ": "29", "#Surf*SurfTen": "30", "T-rest": "31"}}

    def test_edr(self):
        """
        Test of HTPolyNet.gromacs.get_energy_menu on three different edr files
        """
        for i,o in self.edr_menu_results.items():
            logger.info(f'{i}')
            test_file_path_str = str(importlib.resources.files('tests.unit').joinpath(f'fixtures/{i}'))
            test_data = get_energy_menu(test_file_path_str)
            self.assertEqual(test_data,o)
    
