"""Shared setup: the bundled example spectrum, prepared through the early steps."""

import os

import pytest

from venison.core import analysis as A
from venison.core.spectrum import spectrum

EXAMPLE = os.path.join(os.path.dirname(__file__), "..", "data", "raw", "example_2dDEER")


def prepare_example(model="3D", background_start=None, lcurve=False):
    """Example spectrum phased (auto), cutoff at the end, zero time fit, and
    background fit; optionally through the L-curve at the LOOCV alpha."""
    s = spectrum(EXAMPLE + ".DTA")
    signal = A.complex_signal(s)
    cutoff = s.raw_time[-1]
    A.apply_phase_and_cutoff(s, A.auto_phase(signal, s.raw_time, cutoff), cutoff)
    A.fit_zero_point(s, A.symmetry_point(s.raw_time))
    start = s.background_start if background_start is None else background_start
    A.fit_background(s, start, model)
    if lcurve:
        A.select_alpha(s, A.compute_lcurve(s).loocv_index)
    return s


@pytest.fixture(autouse=True)
def no_update_check(monkeypatch):
    """Tests never ask GitHub for the latest release."""
    monkeypatch.setenv("VENISON_NO_UPDATE_CHECK", "1")


@pytest.fixture
def example_bytes():
    """The example's raw .DSC and .DTA contents, for ingest tests."""
    with open(EXAMPLE + ".DSC", "rb") as dsc, open(EXAMPLE + ".DTA", "rb") as dta:
        return dsc.read(), dta.read()
