"""

.. module:: test_relax_instrumentation
   :synopsis: tests the CURE.relax observability -- reactive-species mobility
              (Varshney's criterion) and per-stage density

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger=logging.getLogger(__name__)
import numpy as np
import pandas as pd
from htpolynet.cure.curecontroller import (CureController, minimum_image_displacements,
                                           varshney_criterion, _relax_stage_density)

class TestMinimumImageDisplacements(unittest.TestCase):
    def test_straight_displacement(self):
        box=np.eye(3)*5.0
        r=minimum_image_displacements([[0,0,0],[1,1,1]],[[0.3,0,0],[1,1,1.4]],box)
        self.assertTrue(np.allclose(r,[0.3,0.4]))

    def test_wrapped_atom_is_not_a_box_length_jump(self):
        # an atom at 4.9 that moves to 0.1 has moved 0.2 nm across the boundary,
        # not 4.8 nm back across the box.  Getting this wrong would report a
        # frozen network as maximally mobile.
        box=np.eye(3)*5.0
        r=minimum_image_displacements([[4.9,0,0]],[[0.1,0,0]],box)
        self.assertAlmostEqual(r[0],0.2,places=10)

    def test_accepts_a_box_diagonal(self):
        r=minimum_image_displacements([[4.9,0,0]],[[0.1,0,0]],np.array([5.0,5.0,5.0]))
        self.assertAlmostEqual(r[0],0.2,places=10)

    def test_unset_box_does_not_produce_nan(self):
        # a zero box would divide by zero and make every displacement nan; the
        # raw displacement is the right answer for a non-periodic snapshot
        r=minimum_image_displacements([[0,0,0]],[[0.3,0,0]],np.zeros((3,3)))
        self.assertAlmostEqual(r[0],0.3,places=10)

    def test_empty_input(self):
        r=minimum_image_displacements(np.empty((0,3)),np.empty((0,3)),np.eye(3)*5.0)
        self.assertEqual(r.shape,(0,))

    def test_shape_mismatch_asserts(self):
        with self.assertRaises(AssertionError):
            minimum_image_displacements([[0,0,0]],[[0,0,0],[1,1,1]],np.eye(3)*5.0)

class TestVarshneyCriterion(unittest.TestCase):
    def test_crossing_fraction_counts_the_capture_radius(self):
        s=varshney_criterion([0.1,0.4,0.6,0.9],0.5)
        self.assertEqual(s['n'],4)
        self.assertAlmostEqual(s['crossing_fraction'],0.5)
        self.assertAlmostEqual(s['max'],0.9)

    def test_boundary_is_inclusive(self):
        s=varshney_criterion([0.5],0.5)
        self.assertAlmostEqual(s['crossing_fraction'],1.0)

    def test_rmsd_is_not_the_mean(self):
        # rmsd weights the tail; a network with a few mobile species and many
        # frozen ones must not look uniformly mobile
        s=varshney_criterion([0.0,0.0,0.0,1.0],0.5)
        self.assertAlmostEqual(s['mean'],0.25)
        self.assertAlmostEqual(s['rmsd'],0.5)

    def test_fully_cured_box_reports_no_reactive_species(self):
        s=varshney_criterion([],0.5)
        self.assertEqual(s['n'],0)
        self.assertIsNone(s['rmsd'])
        self.assertIsNone(s['crossing_fraction'])

class TestReactivePositions(unittest.TestCase):
    def setUp(self):
        self.C=CureController.__new__(CureController)
        self.C.dicts={'controls':{'search_radius':0.5}}

    def _tc(self,z,pos,box=np.eye(3)*5.0):
        class _Coords:
            pass
        class _TC:
            pass
        c=_Coords()
        c.A=pd.DataFrame({
            'globalIdx':list(range(1,len(z)+1)),
            'z':z,
            'posX':[p[0] for p in pos],
            'posY':[p[1] for p in pos],
            'posZ':[p[2] for p in pos],
        })
        c.box=box
        tc=_TC()
        tc.Coordinates=c
        return tc

    def test_selects_only_still_reactive_atoms(self):
        tc=self._tc([1,0,2,0],[[0,0,0],[1,0,0],[2,0,0],[3,0,0]])
        idx,pos=self.C._reactive_positions(tc)
        self.assertEqual(list(idx),[1,3])
        self.assertEqual(pos.shape,(2,3))

    def test_none_when_nothing_is_reactive(self):
        tc=self._tc([0,0],[[0,0,0],[1,0,0]])
        self.assertIsNone(self.C._reactive_positions(tc))

    def test_none_when_z_is_absent(self):
        tc=self._tc([1,1],[[0,0,0],[1,0,0]])
        tc.Coordinates.A=tc.Coordinates.A.drop(columns=['z'])
        self.assertIsNone(self.C._reactive_positions(tc))

    def test_report_matches_atoms_by_globalIdx_not_row_order(self):
        # the snapshot is taken before the ladder and compared after it; the
        # pairing has to survive any reordering of the atom dataframe
        tc=self._tc([1,0,1],[[0,0,0],[9,9,9],[2,0,0]])
        snap=self.C._reactive_positions(tc)
        moved=self._tc([1,0,1],[[0.3,0,0],[9,9,9],[2.4,0,0]])
        moved.Coordinates.A=moved.Coordinates.A.iloc[::-1].reset_index(drop=True)
        with self.assertLogs('htpolynet.cure.curecontroller',level='INFO') as cm:
            self.C._report_reactive_mobility(moved,snap)
        msg='\n'.join(cm.output)
        self.assertIn('Reactive-species mobility over relax',msg)
        # displacements 0.3 and 0.4 -> rmsd 0.354, neither crosses 0.5
        self.assertIn('0.354',msg)
        self.assertIn('0% of 2',msg)

    def test_low_mobility_warns(self):
        tc=self._tc([1,1],[[0,0,0],[2,0,0]])
        snap=self.C._reactive_positions(tc)
        moved=self._tc([1,1],[[0.01,0,0],[2.01,0,0]])
        with self.assertLogs('htpolynet.cure.curecontroller',level='WARNING') as cm:
            self.C._report_reactive_mobility(moved,snap)
        self.assertIn('no longer satisfying the criterion',' '.join(cm.output))

    def test_mobile_network_does_not_warn(self):
        tc=self._tc([1,1],[[0,0,0],[2,0,0]])
        snap=self.C._reactive_positions(tc)
        moved=self._tc([1,1],[[0.9,0,0],[2.9,0,0]])
        with self.assertLogs('htpolynet.cure.curecontroller',level='INFO') as cm:
            self.C._report_reactive_mobility(moved,snap)
        self.assertFalse([r for r in cm.records if r.levelno>=logging.WARNING])

    def test_a_broken_snapshot_is_swallowed(self):
        # instrumentation must never be able to fail a build
        tc=self._tc([1,1],[[0,0,0],[2,0,0]])
        snap=(np.array([1,2]),np.empty((0,3)))
        self.C._report_reactive_mobility(tc,snap)

class TestRelaxStageDensity(unittest.TestCase):
    def test_missing_edr_returns_none(self):
        self.assertIsNone(_relax_stage_density('no-such-edr-anywhere'))

class TestRelaxIncrementDefault(unittest.TestCase):
    def test_relax_increment_is_usable_without_a_config(self):
        # _distance_attenuation derives its stage count as int(maxL/increment)
        # with no guard, so a 0.0 default raised ZeroDivisionError at the first
        # relax for any config that omitted the key.
        C=CureController()
        inc=C.dicts['relax']['increment']
        self.assertGreater(inc,0.0)
        self.assertAlmostEqual(inc,0.08)
        self.assertEqual(int(0.5/inc),6)

    def test_drag_increment_may_stay_zero(self):
        # drag's 0.0 is a sentinel, not a bug: the guard in __init__ leaves
        # dragging disabled unless a limit is also set, so nothing divides.
        C=CureController()
        self.assertEqual(C.dicts['drag']['increment'],0.0)
        self.assertFalse(C.dragging_enabled)
