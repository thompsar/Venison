"""End-to-end run of the notebook workflow on the bundled example spectrum."""

import os

import numpy as np
import pytest

from venison.core import analysis as A
from venison.core.export import export_results
from venison.core.spectrum import spectrum

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "example_2dDEER.DTA")


@pytest.fixture(scope="module")
def spec():
    """Spectrum carried through every step with the notebook's default choices."""
    s = spectrum(EXAMPLE)
    signal = A.complex_signal(s)
    cutoff = s.raw_time[-1]
    A.apply_phase_and_cutoff(s, A.auto_phase(signal, s.raw_time, cutoff), cutoff)
    A.fit_zero_point(s, A.symmetry_point(s.raw_time))
    A.fit_background(s, s.background_start, "3D")
    return s


def test_load(spec):
    assert spec.real2D.shape[0] >= 1
    assert spec.slice_mask.all()


def test_phase_minimizes_imaginary(spec):
    signal = A.complex_signal(spec)
    _, imag = A.phased(signal, spec.phase)
    _, imag_off = A.phased(signal, spec.phase + 0.5)
    window = slice(len(signal) // 8, None)
    assert abs(imag[window].mean()) < abs(imag_off[window].mean())


def test_waveform_normalized(spec):
    assert spec.waveform.shape == spec.t.shape
    assert spec.waveform[0] == pytest.approx(1.0, abs=0.1)


def test_full_pipeline(spec, tmp_path):
    lc = A.compute_lcurve(spec)
    assert 0 <= lc.loocv_index < len(spec.alphas)
    A.select_alpha(spec, lc.loocv_index)

    v = A.validate_background(spec)
    lo, hi = v.well
    assert lo <= v.min_index <= max(hi, lo)
    # A local minimum (the first well), not the global one: RMSD keeps
    # falling in the tail where few points remain to fit.
    assert v.rmsd[v.min_index] <= v.rmsd[max(v.min_index - 1, 0)]
    assert v.rmsd[v.min_index] <= v.rmsd[min(v.min_index + 1, len(v.rmsd) - 1)]
    A.select_validated_background(spec, v.min_index, v.well)

    _, _, idx = A.stability_limits(spec)
    A.apply_stability_cutoff(spec, idx)
    assert np.trapezoid(spec.stable_solution) == pytest.approx(1.0)

    gm = A.GaussianModels(spec)
    gm.fit_to_distribution(0)
    gm.fit_to_distribution(1)
    gm.optimize_pending()
    assert gm.optimized[:2].all() and np.isfinite(gm.bics[:2]).all()
    assert gm.best in (0, 1)

    A.run_monte_carlo(spec, gm.best, A.default_rs_weights(spec, gm.best), iterations=400)
    assert spec.model_landscape.shape[1:] == (gm.best + 1, 3)
    lo_env, hi_env = A.best_monte_carlo_model(spec)
    assert np.all(lo_env <= hi_env)

    written = export_results(spec, str(tmp_path / "results"), str(tmp_path / "processed"))
    assert all(os.path.exists(p) for p in written)
    reloaded = spectrum(written[-1])
    np.testing.assert_allclose(reloaded.stable_solution, spec.stable_solution)
