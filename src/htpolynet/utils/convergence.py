"""Deciding when a fluctuating MD observable has settled.

Author: Cameron F. Abrams <cfa22@drexel.edu>

The question this answers is "has the box stopped changing?", and the naive
answer -- take the standard error of the mean over the trace -- is wrong for an
NPT density, which is autocorrelated over hundreds of steps.  Successive frames
are not independent samples, so a block-means SEM is optimistic by roughly a
factor of 1.5 and a gate built on it declares victory early.

The criterion here is the autocorrelation-corrected SEM,
``sigma/sqrt(N/tau_int)``, ported from pestifer's ``density_convergence``.  It
is size-aware for free: ``sigma/mean`` scales as ``1/sqrt(N_atoms)`` while
``tau_int`` is roughly size-independent, so the same tolerance means the same
thing for a small box and a large one.

Two things are deliberately *not* here.  Pestifer's chunking
(``next_chunk_steps``, ``is_patch_grid_crash``) exists because NAMD fixes its
patch and PME grids at the start of each ``run``; GROMACS rescales the box
within one ``mdrun`` and has no such failure mode, so that machinery is most of
pestifer's complexity and none of its value.  And a small SEM alone is not
taken as convergence: a trace can drift steadily while every short window looks
tight, so the residual drift across the window is tested as well.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)


def autocorrelation(x, max_lag=None):
    """Normalized autocorrelation of a series, lag 0 upward.

    Uses the biased estimator -- every lag divided by the full sample count
    rather than by the number of overlapping pairs -- because the unbiased one
    grows noisier with lag and can drive the integrated time negative.

    Args:
        x (array-like): the series
        max_lag (int): highest lag to compute, defaults to half the series

    Returns:
        numpy.ndarray: rho[k] for k = 0 .. max_lag, with rho[0] = 1
    """
    a = np.asarray(x, dtype=float).reshape(-1)
    n = a.size
    if n < 2:
        return np.ones(1)
    d = a - a.mean()
    denom = float((d ** 2).sum())
    if denom <= 0.0:                       # a perfectly flat series is converged
        return np.ones(1)
    if max_lag is None:
        max_lag = n // 2
    max_lag = int(min(max_lag, n - 1))
    return np.array([float((d[:n - k] * d[k:]).sum()) / denom for k in range(max_lag + 1)])


def integrated_autocorrelation_time(x):
    """Integrated autocorrelation time, by the initial-positive-sequence rule.

    ``tau_int = 1 + 2*sum(rho_k)``, summed until the first non-positive
    ``rho_k``; truncating there is what keeps the estimate from accumulating
    noise at long lag.  Never returns less than 1, which is the uncorrelated
    case.

    Args:
        x (array-like): the series

    Returns:
        float: tau_int, in samples
    """
    rho = autocorrelation(x)
    tau = 1.0
    for k in range(1, rho.size):
        if rho[k] <= 0.0:
            break
        tau += 2.0 * rho[k]
    return max(tau, 1.0)


def corrected_sem(x):
    """Standard error of the mean, corrected for autocorrelation.

    Args:
        x (array-like): the series

    Note that this cannot itself detect "the window is too short to judge".
    Truncating the autocorrelation sum at lag ``n//2`` caps ``tau_int`` at
    roughly ``n/3``, so ``n_effective`` floors near 3 even for a pure ramp or
    step -- a guard on it would never fire.  The protection against judging a
    short window is ``min_samples`` in :func:`series_converged`.

    Returns:
        tuple: (sem, tau_int, n_effective); sem is None only for fewer than two samples
    """
    a = np.asarray(x, dtype=float).reshape(-1)
    n = a.size
    if n < 2:
        return None, 1.0, float(n)
    tau = integrated_autocorrelation_time(a)
    n_eff = n / tau
    return float(a.std(ddof=1) / np.sqrt(n_eff)), tau, n_eff


def series_converged(x, tolerance, min_samples=50, drift_sems=2.0):
    """Decides whether a trailing window of an observable has settled.

    Both tests must pass.  The **precision** test asks whether the mean is
    known to better than ``tolerance``, using the autocorrelation-corrected
    SEM.  The **drift** test compares the first and second halves of the
    window and requires their difference to be within ``drift_sems`` standard
    errors, because a trace that is still climbing can show a small SEM in
    every short window while going nowhere near equilibrium.

    Args:
        x (array-like): the trailing window of the observable
        tolerance (float): largest acceptable corrected SEM, in the observable's units
        min_samples (int): refuse to judge a window shorter than this
        drift_sems (float): half-to-half difference allowed, in units of that difference's own standard error

    Returns:
        dict: 'converged' (bool), 'reason' (str), 'mean', 'sem', 'tau_int', 'n', 'n_effective', 'drift', 'sem_drift'
    """
    a = np.asarray(x, dtype=float).reshape(-1)
    out = {'converged': False, 'reason': '', 'mean': None, 'sem': None,
           'tau_int': None, 'n': int(a.size), 'n_effective': None, 'drift': None,
           'sem_drift': None}
    if a.size < max(int(min_samples), 4):
        out['reason'] = f'only {a.size} samples; need {max(int(min_samples), 4)}'
        return out
    sem, tau, n_eff = corrected_sem(a)
    out.update({'mean': float(a.mean()), 'sem': sem, 'tau_int': float(tau),
                'n_effective': float(n_eff)})
    half = a.size // 2
    first, second = a[:half], a[half:]
    drift = float(second.mean() - first.mean())
    out['drift'] = drift
    if sem > tolerance:
        out['reason'] = f'sem {sem:.3f} exceeds tolerance {tolerance:.3f}'
        return out
    # The drift is a difference of two half-window means, so its own standard
    # error is the quadrature sum of theirs -- about twice the whole-window
    # sem, since each half has half the samples.  Comparing the drift against
    # the whole-window sem instead would reject settled noise roughly half the
    # time, which is the difference between a gate and a coin toss.
    sem_first = corrected_sem(first)[0]
    sem_second = corrected_sem(second)[0]
    if sem_first is None or sem_second is None:
        sem_drift = 2.0 * sem
    else:
        sem_drift = float(np.hypot(sem_first, sem_second))
    out['sem_drift'] = sem_drift
    if abs(drift) > drift_sems * sem_drift:
        out['reason'] = (f'still drifting: halves differ by {drift:+.3f}, '
                         f'more than {drift_sems:g} x its standard error {sem_drift:.3f}')
        return out
    out['converged'] = True
    out['reason'] = (f'sem {sem:.3f} within {tolerance:.3f}, drift {drift:+.3f} '
                     f'within {drift_sems:g} x {sem_drift:.3f}')
    return out
