"""Background validation: the suggested range, and Fractal handling."""

import contextlib
import io
import os

import numpy as np

from venison.core import analysis as A
from venison.core.spectrum import spectrum

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "example_2dDEER.DTA")
X = np.linspace(0.2, 5, 100)


def test_suggests_first_well():
    # Early starts fit badly, a well around 1 µs, a hump, then a misleading
    # low tail from very late starts.
    rmsd = 1 + 3 * np.exp(-X / 0.3) + 0.3 * np.exp(-((X - 2.2) / 0.4) ** 2) - 0.5 * (X > 4)
    m, (lo, hi) = A.suggest_range(rmsd, dof=600)
    assert 0.5 < X[m] < 2.0
    assert lo <= m <= hi and X[hi] < 2.2


def test_minimum_on_first_trial_is_skipped():
    # Rises straight from the first trial, then a real well near 3 µs.
    rmsd = 2 + np.sin(X * 1.1 + 0.2)
    m, (lo, hi) = A.suggest_range(rmsd, dof=600)
    assert m > 0 and 2.5 < X[m] < 4.5
    assert lo <= m <= hi


def _prepared(model):
    s = spectrum(EXAMPLE)
    sig = A.complex_signal(s)
    c = s.raw_time[-1]
    A.apply_phase_and_cutoff(s, A.auto_phase(sig, s.raw_time, c), c)
    A.fit_zero_point(s, A.symmetry_point(s.raw_time))
    A.fit_background(s, s.background_start, model)
    A.select_alpha(s, A.compute_lcurve(s).loocv_index)
    return s


def test_fractal_dimension_held_fixed():
    s = _prepared("Fractal")
    d = s.bg_tab_fractal_d
    assert d is not None
    with contextlib.redirect_stdout(io.StringIO()):
        v = A.validate_background(s, n_trials=30, hold_fractal_d=True)
    assert s.validation_settings["held_d"] == d
    assert s.validation_settings["fallbacks"] == 0
    assert getattr(s, "fractal_fixed_d", None) is None  # not left on for later fits
    lo, hi = v.well
    assert lo <= v.min_index <= hi
    # Sweep honors the µs settings.
    assert v.fit_points[0] == np.float64(s.raw_time[A.nearest_index(s.raw_time, 0.2)])
    assert v.fit_points[-1] <= s.cutoff - 0.8 + 1e-9
