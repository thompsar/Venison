"""The fast regularization used by the Stability tab must match the original."""

import os

import numpy as np

from venison.core import analysis as A
from venison.core.spectrum import spectrum

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "example_2dDEER.DTA")


def test_fast_regularization_matches_original():
    s = spectrum(EXAMPLE)
    sig = A.complex_signal(s)
    c = s.raw_time[-1]
    A.apply_phase_and_cutoff(s, A.auto_phase(sig, s.raw_time, c), c)
    A.fit_zero_point(s, A.symmetry_point(s.raw_time))
    A.fit_background(s, 1.0)
    s.alpha = s.alphas[26]
    s.solution = s.alpha_regularize(s.alpha)

    previous = None
    for idx in range(380, 300, -8):   # a drag: each step warm-starts from the last
        result = A.stability_cutoff(s, idx)
        assert result["instability_index"] == idx
        exact = s.waveform_regularize(result["stable_waveform"])
        exact = exact / np.trapezoid(exact)
        np.testing.assert_allclose(result["stable_solution"], exact, atol=1e-12)
        # Warm start from an unrelated solution must still give the exact answer.
        cold = A.regularize_fast(s, result["stable_waveform"], warm=previous)
        np.testing.assert_allclose(cold / np.trapezoid(cold), exact, atol=1e-12)
        s.stable_solution = result["stable_solution"]
        previous = np.roll(result["stable_solution"], 40)
