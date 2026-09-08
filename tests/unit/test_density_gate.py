"""

.. module:: test_density_gate
   :synopsis: tests the densification density-convergence gate

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import logging
logger = logging.getLogger(__name__)
import numpy as np
from htpolynet.core.topocoord import TopoCoord


def drifting(lo, hi, n=600, seed=0):
    return np.linspace(lo, hi, n) + np.random.default_rng(seed).normal(0, 1.0, n)


def settled(mean=1100.0, sd=3.0, n=600, seed=1):
    return np.random.default_rng(seed).normal(mean, sd, n)


class _Harness(TopoCoord):
    """A TopoCoord that runs no MD: mdrun is recorded, densities are scripted.

    The gate's decisions are pure functions of the density traces, so the whole
    control flow is testable on a runner with no GROMACS.
    """
    def __init__(self, segments):
        self.calls = []
        self._segments = list(segments)
        self._n = 0

    def grompp_and_mdrun(self, **kwargs):
        self.calls.append(kwargs['out'])

    def _log_box(self):
        pass

    def _density_series(self, edr_pfx, gromacs_dict={}):
        if self._n >= len(self._segments):
            return None
        s = self._segments[self._n]
        self._n += 1
        return s


class TestConvergeDensity(unittest.TestCase):
    def test_extends_until_settled(self):
        h = _Harness([drifting(1080, 1098), settled()])
        with self.assertLogs('htpolynet.core.topocoord', level='INFO') as cm:
            extra = h._converge_density('dens', 'npt', {'tolerance': 1.0, 'max_repeats': 5},
                                        drifting(1000, 1080))
        self.assertEqual(extra, ['dens-converge-1-npt', 'dens-converge-2-npt'])
        self.assertEqual(h.calls, extra)
        self.assertIn('Density converged after 2 extension(s)', '\n'.join(cm.output))

    def test_already_settled_runs_nothing(self):
        h = _Harness([])
        with self.assertLogs('htpolynet.core.topocoord', level='INFO') as cm:
            extra = h._converge_density('dens', 'npt', {'tolerance': 1.0}, settled())
        self.assertEqual(extra, [])
        self.assertEqual(h.calls, [])
        self.assertIn('Density converged after 0 extension(s)', '\n'.join(cm.output))

    def test_ceiling_warns_and_says_not_to_trust_the_number(self):
        # the outcome that matters: a build that never settles must say so,
        # not silently report whichever density it happened to reach
        h = _Harness([drifting(1000, 1100, seed=2), drifting(1100, 1200, seed=3)])
        with self.assertLogs('htpolynet.core.topocoord', level='WARNING') as cm:
            extra = h._converge_density('dens', 'npt', {'tolerance': 0.01, 'max_repeats': 2},
                                        drifting(900, 1000))
        msg = '\n'.join(cm.output)
        self.assertEqual(len(extra), 2)
        self.assertIn('did NOT converge', msg)
        self.assertIn('unsettled box', msg)

    def test_max_repeats_is_a_hard_ceiling(self):
        h = _Harness([drifting(1000, 1100, seed=i) for i in range(20)])
        extra = h._converge_density('dens', 'npt', {'tolerance': 0.01, 'max_repeats': 3},
                                    drifting(900, 1000))
        self.assertEqual(len(extra), 3)
        self.assertEqual(len(h.calls), 3)

    def test_zero_ceiling_judges_but_never_extends(self):
        h = _Harness([])
        with self.assertLogs('htpolynet.core.topocoord', level='WARNING'):
            extra = h._converge_density('dens', 'npt', {'tolerance': 0.01, 'max_repeats': 0},
                                        drifting(900, 1000))
        self.assertEqual(extra, [])
        self.assertEqual(h.calls, [])

    def test_no_density_trace_does_not_gate(self):
        # gmx energy could not produce a Density column; refuse to gate rather
        # than loop forever or claim convergence
        h = _Harness([])
        with self.assertLogs('htpolynet.core.topocoord', level='WARNING') as cm:
            extra = h._converge_density('dens', 'npt', {'tolerance': 1.0}, None)
        self.assertEqual(extra, [])
        self.assertEqual(h.calls, [])
        self.assertIn('not gating', '\n'.join(cm.output))

    def test_defaults_are_applied_when_the_block_is_bare(self):
        h = _Harness([])
        with self.assertLogs('htpolynet.core.topocoord', level='INFO'):
            extra = h._converge_density('dens', 'npt', {}, settled())
        self.assertEqual(extra, [])
