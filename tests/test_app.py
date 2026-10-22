"""App-level behavior, driven headlessly (no browser)."""

import asyncio
import types

import numpy as np
import pytest

import venison.app as app_module
from venison.app import VenisonApp, configure_panel
from venison.core import analysis as A
from venison.core.project import create_project
from venison.core.session import SpectrumSession


@pytest.fixture
def app(tmp_path, monkeypatch, example_bytes):
    monkeypatch.setattr(app_module, "CONFIG_PATH", tmp_path / "config.json")
    configure_panel()
    project = create_project(tmp_path, "Demo")
    project.add_spectrum("first", *example_bytes)
    project.add_spectrum("second", *example_bytes)
    a = VenisonApp()
    a.open_project(project.path)
    return a


def test_cannot_switch_spectrum_mid_computation(app):
    first = app.session
    app.busy = "Validating background…"
    app.open_spectrum("second")
    assert app.session is first and app._spectrum_select.value == "first"
    app.busy = ""
    app.open_spectrum("second")
    assert app.session.name == "second"


def test_failed_open_clears_picker(app, monkeypatch):
    def broken(*args, **kwargs):
        raise OSError("corrupt file")
    monkeypatch.setattr(app_module, "SpectrumSession", broken)
    app.open_spectrum("second")
    assert app.session is None and app._spectrum_select.value is None


def test_stability_refit_for_another_spectrum_is_dropped(app):
    tab = app.tabs[5]
    old = app.session
    spec = old.spec
    sig = A.complex_signal(spec)
    A.apply_phase_and_cutoff(spec, A.auto_phase(sig, spec.raw_time, spec.raw_time[-1]),
                             spec.raw_time[-1])
    A.fit_zero_point(spec, A.symmetry_point(spec.raw_time))
    A.fit_background(spec, 1.0)
    spec.alpha = spec.alphas[26]
    spec.solution = spec.alpha_regularize(spec.alpha)
    result = A.stability_cutoff(spec, 300)

    app.open_spectrum("second")           # switched while the refit was running
    tab._target, tab._save_pending = 300, True
    tab._show_refit(old, 300, result)
    assert not hasattr(app.session.spec, "stable_solution")
    assert not app.session.is_done("stability") and not tab._save_pending


def test_ingest_accepts_dsc_decoded_as_text(app, example_bytes):
    dsc, dta = example_bytes
    tab = app.tabs[0]
    tab.refresh()
    # FileDropper hands over text/* files as str.
    tab._on_drop(types.SimpleNamespace(new={"third.DSC": dsc.decode("utf-8"), "third.DTA": dta}))
    tab._on_add(None)
    assert (app.project.path / "raw_data" / "third.DSC").read_bytes() == dsc
    assert app.session.name == "third"
    assert np.isfinite(app.session.spec.real).all()


def test_result_discarded_if_spectrum_changed_mid_task(app, capsys):
    applied = []

    def work():
        # A switch slipping in before the result is applied (e.g. in the tick
        # between the task ending and its follow-up running).
        app.session = SpectrumSession(app.project, "second")
        return 42

    asyncio.run(app.run_task("Working", work, applied.append))
    assert applied == [] and "discarded" in capsys.readouterr().out
    assert not app._spectrum_select.disabled and not app.busy

    asyncio.run(app.run_task("Working", lambda: 7, applied.append))
    assert applied == [7]


def test_failed_refit_reruns_a_newer_slider_position(app):
    tab = app.tabs[5]
    reruns = []
    tab._update = lambda save: reruns.append(save)
    tab._target, tab._save_pending = 200, True     # released at 200 meanwhile
    tab._show_refit(app.session, 150, RuntimeError("singular"))
    assert reruns == [False] and tab._save_pending  # 200 gets refit and still saved
