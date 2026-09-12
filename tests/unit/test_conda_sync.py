"""

.. module:: test_conda_sync
   :synopsis: tests the release preflight that compares pyproject.toml
              against the conda-forge recipe

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import importlib.util
import sys
import unittest
from pathlib import Path

import pytest

# check-conda-sync.py is release tooling, run on a maintainer's machine, and it
# uses tomllib, which is stdlib only from Python 3.11.  htpolynet itself supports
# 3.10, and CI tests 3.10, so importing the script unconditionally made the
# whole 3.10 job fail at collection -- which also meant nothing else ran there.
if sys.version_info < (3, 11):
    pytest.skip('check-conda-sync.py needs tomllib (Python 3.11+)', allow_module_level=True)

_SCRIPT = Path(__file__).resolve().parents[2] / 'scripts' / 'check-conda-sync.py'
_spec = importlib.util.spec_from_file_location('check_conda_sync', _SCRIPT)
ccs = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(ccs)


RECIPE = """
package:
  name: htpolynet
requirements:
  host:
    - python {{ python_min }}
    - pip
  run:
    - python >={{ python_min }}
    - numpy >=1.24
    - matplotlib-base
    - rdkit >=2024.03
    - ambertools

test:
  imports:
    - htpolynet
"""


class TestFloorParsing(unittest.TestCase):
    def test_extracts_a_lower_bound(self):
        self.assertEqual(ccs._floor('numpy>=1.24'), '1.24')
        self.assertEqual(ccs._floor('python-graphviz >=0.20'), '0.20')

    def test_no_bound(self):
        self.assertIsNone(ccs._floor('matplotlib-base'))
        self.assertIsNone(ccs._floor('setuptools'))

    def test_ignores_an_upper_bound(self):
        # an upper bound in the recipe that pyproject lacks is a deliberate
        # conda-side pin, not drift
        self.assertEqual(ccs._floor('pandas >=2,<3'), '2')

    def test_jinja_pin_is_not_a_version(self):
        self.assertIsNone(ccs._floor('python >={{ python_min }}'))


class TestVersionComparison(unittest.TestCase):
    def test_zero_padded_segments_compare_equal(self):
        # conda and pip agree 2024.03 == 2024.3; a string compare does not,
        # and would have flagged the rdkit pin as drift forever
        self.assertFalse(ccs._floor_is_weaker('2024.03', '2024.3'))
        self.assertFalse(ccs._floor_is_weaker('2024.3', '2024.03'))

    def test_numeric_not_lexicographic(self):
        self.assertFalse(ccs._floor_is_weaker('1.24', '1.9'))
        self.assertTrue(ccs._floor_is_weaker('1.9', '1.24'))

    def test_absent_bound_is_weaker(self):
        self.assertTrue(ccs._floor_is_weaker(None, '3.5'))

    def test_equal_is_not_weaker(self):
        self.assertFalse(ccs._floor_is_weaker('3.5', '3.5'))

    def test_higher_is_not_weaker(self):
        # a stricter recipe is safe; only a looser one lets conda solve an
        # environment pip would have refused
        self.assertFalse(ccs._floor_is_weaker('4.0', '3.5'))

    def test_unequal_segment_counts(self):
        self.assertFalse(ccs._floor_is_weaker('2.0', '2'))
        self.assertFalse(ccs._floor_is_weaker('2', '2.0'))
        self.assertTrue(ccs._floor_is_weaker('2', '2.0.1'))


class TestRecipeParsing(unittest.TestCase):
    def test_names_and_floors(self):
        deps = ccs.parse_recipe_run_deps(RECIPE)
        self.assertEqual(deps['numpy'], '1.24')
        self.assertIsNone(deps['matplotlib-base'])
        self.assertEqual(deps['rdkit'], '2024.03')
        self.assertIn('ambertools', deps)

    def test_host_section_is_not_read_as_run(self):
        deps = ccs.parse_recipe_run_deps(RECIPE)
        self.assertNotIn('pip', deps)

    def test_stops_at_the_next_section(self):
        deps = ccs.parse_recipe_run_deps(RECIPE)
        self.assertNotIn('htpolynet', deps)


class TestCompare(unittest.TestCase):
    def test_missing_name(self):
        missing, extra, weaker = ccs.compare({'ycleptic': '2.4.1'}, {'numpy': '1.24'})
        self.assertIn('ycleptic', missing)

    def test_weaker_floor_is_reported(self):
        # the case grayskull caught on feedstock PR #21 and this script did not
        missing, extra, weaker = ccs.compare({'matplotlib': '3.5'},
                                             {'matplotlib-base': None})
        self.assertEqual(missing, set())
        self.assertEqual(weaker, [('matplotlib-base', None, '3.5')])

    def test_name_remap_is_applied_to_floors_too(self):
        _, _, weaker = ccs.compare({'graphviz': '0.20'}, {'python-graphviz': '0.1'})
        self.assertEqual(weaker, [('python-graphviz', '0.1', '0.20')])

    def test_equal_floors_are_silent(self):
        missing, extra, weaker = ccs.compare({'numpy': '1.24'}, {'numpy': '1.24'})
        self.assertEqual((missing, extra, weaker), (set(), set(), []))

    def test_conda_only_deps_are_not_extra(self):
        _, extra, _ = ccs.compare({'numpy': '1.24'},
                                  {'numpy': '1.24', 'ambertools': None,
                                   'python': None, 'graphviz': None})
        self.assertEqual(extra, set())

    def test_pyproject_dep_without_a_floor_is_not_judged(self):
        _, _, weaker = ccs.compare({'setuptools': None}, {'setuptools': None})
        self.assertEqual(weaker, [])
