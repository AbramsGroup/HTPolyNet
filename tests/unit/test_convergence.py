"""

.. module:: test_convergence
   :synopsis: tests the autocorrelation-corrected convergence criterion

.. moduleauthor: Cameron F. Abrams, <cfa22@drexel.edu>

"""
import unittest
import numpy as np
from htpolynet.utils.convergence import (autocorrelation, corrected_sem,
                                         integrated_autocorrelation_time,
                                         series_converged)


def ar1(n, rho, sigma=1.0, mean=0.0, seed=0):
    """An AR(1) series, whose true tau_int is (1+rho)/(1-rho)."""
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.normal(0.0, sigma / np.sqrt(1 - rho ** 2))
    for i in range(1, n):
        x[i] = rho * x[i - 1] + rng.normal(0.0, sigma)
    return x + mean


class TestAutocorrelation(unittest.TestCase):
    def test_lag_zero_is_one(self):
        self.assertAlmostEqual(autocorrelation(ar1(500, 0.8))[0], 1.0)

    def test_flat_series_is_handled(self):
        # a perfectly flat trace has no variance to normalize by; it must not
        # divide by zero, and it is trivially converged
        r = autocorrelation(np.full(100, 3.0))
        self.assertEqual(r.tolist(), [1.0])

    def test_too_short_to_correlate(self):
        self.assertEqual(autocorrelation([1.0]).tolist(), [1.0])


class TestIntegratedAutocorrelationTime(unittest.TestCase):
    def test_white_noise_is_about_one(self):
        rng = np.random.default_rng(1)
        tau = integrated_autocorrelation_time(rng.normal(size=4000))
        self.assertLess(tau, 1.6)

    def test_recovers_the_ar1_value(self):
        # rho 0.8 -> tau_int 9.0; the estimator is noisy, so allow a wide band
        tau = integrated_autocorrelation_time(ar1(20000, 0.8, seed=2))
        self.assertGreater(tau, 6.0)
        self.assertLess(tau, 13.0)

    def test_never_below_one(self):
        self.assertGreaterEqual(integrated_autocorrelation_time([1.0, -1.0, 1.0, -1.0]), 1.0)


class TestCorrectedSem(unittest.TestCase):
    def test_correlated_series_gets_a_larger_sem(self):
        # this is the whole point: the naive SEM understates the error for a
        # correlated trace, so a gate built on it declares victory early
        x = ar1(4000, 0.8, seed=3)
        naive = x.std(ddof=1) / np.sqrt(x.size)
        sem, tau, n_eff = corrected_sem(x)
        self.assertGreater(sem, naive)
        self.assertGreater(tau, 5.0)
        self.assertLess(n_eff, x.size)

    def test_uncorrelated_series_is_close_to_naive(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=4000)
        naive = x.std(ddof=1) / np.sqrt(x.size)
        sem, _, _ = corrected_sem(x)
        self.assertLess(abs(sem - naive) / naive, 0.4)

    def test_a_drifting_window_inflates_tau_and_the_sem(self):
        # a ramp is correlated at every lag the window can see, so tau_int
        # grows with the window and the corrected sem stays large where the
        # naive one would shrink as 1/sqrt(n)
        long_ramp = np.linspace(0.0, 1.0, 400)
        sem, tau, n_eff = corrected_sem(long_ramp)
        naive = long_ramp.std(ddof=1) / np.sqrt(long_ramp.size)
        self.assertGreater(tau, 100.0)
        self.assertGreater(sem, 10 * naive)

    def test_n_effective_floors_near_three(self):
        # documents why corrected_sem has no "too short to judge" guard:
        # truncating the sum at lag n//2 caps tau_int near n/3, so n_effective
        # cannot fall below about 3 however correlated the series is.  A guard
        # on it would be a declaration that never fires; min_samples in
        # series_converged is the real protection.
        for n in (8, 20, 100):
            step = np.concatenate([np.zeros(n // 2), np.ones(n - n // 2)])
            self.assertGreater(corrected_sem(step)[2], 2.0)

    def test_fewer_than_two_samples(self):
        self.assertIsNone(corrected_sem([1.0])[0])


class TestSeriesConverged(unittest.TestCase):
    def test_settled_noise_converges(self):
        rng = np.random.default_rng(6)
        r = series_converged(rng.normal(1100.0, 3.0, size=2000), tolerance=1.0)
        self.assertTrue(r['converged'], r['reason'])
        self.assertAlmostEqual(r['mean'], 1100.0, delta=1.0)

    def test_a_rising_trace_does_not_converge(self):
        # the case the drift test exists for: a box still densifying can show a
        # small SEM while going nowhere near equilibrium
        # the ramp is gentle enough that the SEM test passes; only the drift
        # test can catch it, which is exactly why the drift test is there
        rng = np.random.default_rng(7)
        x = np.linspace(1100.0, 1102.0, 2000) + rng.normal(0, 1.0, size=2000)
        r = series_converged(x, tolerance=5.0)
        self.assertFalse(r['converged'])
        self.assertIn('drifting', r['reason'])

    def test_too_noisy_does_not_converge(self):
        rng = np.random.default_rng(8)
        r = series_converged(rng.normal(1100.0, 50.0, size=200), tolerance=0.1)
        self.assertFalse(r['converged'])
        self.assertIn('exceeds tolerance', r['reason'])

    def test_short_window_refuses_to_judge(self):
        r = series_converged([1100.0] * 10, tolerance=1.0, min_samples=50)
        self.assertFalse(r['converged'])
        self.assertIn('need 50', r['reason'])

    def test_flat_series_converges(self):
        r = series_converged(np.full(200, 1100.0), tolerance=1.0)
        self.assertTrue(r['converged'], r['reason'])
        self.assertEqual(r['sem'], 0.0)

    def test_reports_the_numbers_it_judged_on(self):
        rng = np.random.default_rng(9)
        r = series_converged(rng.normal(1100.0, 3.0, size=2000), tolerance=1.0)
        for k in ('mean', 'sem', 'tau_int', 'n', 'n_effective', 'drift'):
            self.assertIsNotNone(r[k], k)
        self.assertEqual(r['n'], 2000)

    def test_drift_tolerance_is_configurable(self):
        rng = np.random.default_rng(10)
        x = np.linspace(1100.0, 1100.6, 2000) + rng.normal(0, 1.0, size=2000)
        self.assertFalse(series_converged(x, tolerance=1.0, drift_sems=0.1)['converged'])
        self.assertTrue(series_converged(x, tolerance=1.0, drift_sems=50.0)['converged'])

    def test_settled_noise_is_not_rejected_by_the_drift_test(self):
        # regression: comparing the drift against the whole-window sem rather
        # than against the drift's own standard error rejected settled noise
        # about half the time
        n_reject = 0
        for seed in range(30):
            rng = np.random.default_rng(100 + seed)
            r = series_converged(rng.normal(1100.0, 3.0, size=2000), tolerance=1.0)
            if not r['converged']:
                n_reject += 1
        self.assertLessEqual(n_reject, 2, f'{n_reject}/30 settled traces rejected')
