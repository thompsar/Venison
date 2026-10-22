"""The fast regularization used by the Stability tab must match the original."""

import pickle

import numpy as np

from conftest import prepare_example
from venison.core import analysis as A


def test_regularizer_cache_is_not_pickled():
    # Multiprocessing pickles the whole spectrum for every task.
    s = prepare_example(background_start=1.0)
    s.alpha = s.alphas[26]
    size = len(pickle.dumps(s))
    A.prepare_regularizer(s)
    assert len(pickle.dumps(s)) == size
    assert not any(k.startswith("_") for k in pickle.loads(pickle.dumps(s)).__dict__)


def test_fast_regularization_matches_original():
    s = prepare_example(background_start=1.0)
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
