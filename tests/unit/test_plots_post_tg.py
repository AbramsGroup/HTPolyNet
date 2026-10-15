"""

.. module:: test_plots_post_tg
   :synopsis: ``plots post`` says when a Tg depends on the fit window, and needs its config

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import os
import tempfile
import types
import unittest

import numpy as np
import pandas as pd

from htpolynet.analysis.plot import do_tg_plots, post_plots
from htpolynet.analysis.utils import tg_window_range

TEMPS = np.arange(300.0, 605.0, 5.0)          # 61 rungs, as in example 6's ladders
LADDER = {'warmup_ps': 0, 'ps_per_rise': 1, 'ps_per_run': 10, 'deltaT': 5}


def kinked(T, Tg=450.0, glassy=-0.20, rubbery=-0.35, noise=0.0, seed=0):
    rng = np.random.default_rng(seed)
    rho = np.where(T < Tg, 1200.0 + glassy * (T - 300.0),
                   1200.0 + glassy * (Tg - 300.0) + rubbery * (T - Tg))
    return rho + rng.normal(0.0, noise, T.shape)


def write_ladder(root, subdir, rho):
    """ladder.csv with one constant-T, constant-density block per rung."""
    rows = []
    step = LADDER['ps_per_rise'] + LADDER['ps_per_run']
    for k, (T, r) in enumerate(zip(TEMPS, rho)):
        for t in range(step):
            rows.append({'time(ps)': k * step + t, 'Temperature': T, 'Density': r,
                         'Volume': 1.0e5 / r})
    os.makedirs(os.path.join(root, subdir))
    pd.DataFrame(rows).to_csv(os.path.join(root, subdir, 'ladder.csv'), index=False)


class TestWindowRange(unittest.TestCase):
    def test_a_clean_kink_holds_its_tg_across_windows(self):
        lo, hi = tg_window_range(TEMPS, kinked(TEMPS))
        self.assertAlmostEqual(lo, 450.0, delta=1.0)
        self.assertAlmostEqual(hi, 450.0, delta=1.0)

    def test_nearly_parallel_lines_scatter_off_the_scale(self):
        # no visible transition: the two fitted lines are parallel up to noise
        lo, hi = tg_window_range(TEMPS, kinked(TEMPS, glassy=-0.30, rubbery=-0.30,
                                               noise=0.5, seed=0))
        self.assertTrue(lo < 300.0 or hi > 600.0, (lo, hi))

    def test_windows_that_do_not_fit_the_data_are_skipped(self):
        self.assertIsNone(tg_window_range(TEMPS[:15], kinked(TEMPS[:15])))


class TestDoTgPlots(unittest.TestCase):
    def run_tg(self, heat, cool):
        phases = [{'ladder': dict(LADDER, subdir='heat')},
                  {'ladder': dict(LADDER, subdir='cool')}]
        with tempfile.TemporaryDirectory() as d:
            write_ladder(d, 'heat', heat)
            write_ladder(d, 'cool', cool)
            cwd = os.getcwd()
            os.chdir(d)
            try:
                with self.assertLogs('htpolynet.analysis.plot', level='INFO') as logs:
                    result = do_tg_plots(phases, [d])
            finally:
                os.chdir(cwd)
        return result, '\n'.join(logs.output)

    def test_two_sound_ladders_report_their_window_range_and_no_warning(self):
        result, log = self.run_tg(kinked(TEMPS), kinked(TEMPS, Tg=440.0))
        self.assertAlmostEqual(result['heating'], 450.0, delta=1.0)
        self.assertIn('over fit windows', log)
        self.assertNotIn('WARNING', log)
        self.assertEqual(result['unquotable'], [])

    def test_an_unfittable_heating_ladder_is_named(self):
        bad = kinked(TEMPS, glassy=-0.30, rubbery=-0.30, noise=0.5, seed=0)
        result, log = self.run_tg(bad, kinked(TEMPS))
        self.assertEqual(result['unquotable'], ['heating'])
        self.assertIn('no heating Tg can be read off this ladder', log)

    def test_heating_and_cooling_far_apart_is_flagged(self):
        result, log = self.run_tg(kinked(TEMPS, Tg=422.0), kinked(TEMPS, Tg=478.0))
        self.assertIn('differ by 56 K', log)


class TestPostPlotsNeedsItsConfig(unittest.TestCase):
    def test_no_cfg_is_an_error_not_a_silent_success(self):
        args = types.SimpleNamespace(n_points=[10, 20], cfg=[], proj=['proj-0'])
        with self.assertRaisesRegex(ValueError, '--cfg'):
            post_plots(args)


if __name__ == '__main__':
    unittest.main()
