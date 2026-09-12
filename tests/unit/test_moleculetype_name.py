"""

.. module:: test_moleculetype_name
   :synopsis: the whole-system moleculetype name htpolynet writes

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import tempfile
import os
from htpolynet.core.topology import Topology, _WHOLE_SYSTEM_MOLECULETYPE_


def _section_first_row(path, directive):
    """First data row under [ directive ] in a written top file."""
    with open(path) as f:
        lines = f.read().splitlines()
    for i, l in enumerate(lines):
        if l.strip() == f'[ {directive} ]':
            for row in lines[i + 1:]:
                s = row.strip()
                if s and not s.startswith(';'):
                    return s.split()
    return None


class TestWholeSystemMoleculetypeName(unittest.TestCase):
    def test_default_is_not_the_string_none(self):
        # 'None' read to users like an unset field or a bug, and it sits on the
        # exact line people inspect after GROMACS warns about inconsistent shifts
        self.assertEqual(_WHOLE_SYSTEM_MOLECULETYPE_, 'whole_system')
        T = Topology(system_name='htpolynet')
        self.assertEqual(T.D['moleculetype']['name'].iloc[0], 'whole_system')
        self.assertEqual(T.D['molecules']['Compound'].iloc[0], 'whole_system')

    def test_written_names_agree(self):
        # grompp refuses a [ molecules ] entry naming no [ moleculetype ]
        T = Topology(system_name='htpolynet')
        with tempfile.TemporaryDirectory() as td:
            p = os.path.join(td, 'sys.top')
            T.write_top(p)
            mt = _section_first_row(p, 'moleculetype')
            mol = _section_first_row(p, 'molecules')
        self.assertIsNotNone(mt)
        self.assertIsNotNone(mol)
        self.assertEqual(mt[0], mol[0])
        self.assertEqual(mt[0], 'whole_system')

    def test_an_existing_top_keeps_its_own_name(self):
        # topologies written before the rename carry 'None'; reading one back
        # must preserve that rather than silently rename half the file
        src = os.path.join(os.path.dirname(__file__), 'test_topology', 'test.top')
        T = Topology.read_top(src)
        self.assertEqual(str(T.D['moleculetype']['name'].iloc[0]), 'None')
        self.assertEqual(str(T.D['molecules']['Compound'].iloc[0]), 'None')
