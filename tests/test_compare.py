"""Compare tab: overlays spectra that are far enough along, per mode."""

import pytest
from bokeh.models import Legend

import venison.app as app_module
from venison.app import VenisonApp, configure_panel
from venison.core import analysis as A
from venison.core.project import create_project
from venison.core.session import STEPS, SpectrumSession


def _process(project, name, monte_carlo=False):
    """Take a spectrum through Stability (and optionally Monte Carlo), quickly."""
    session = SpectrumSession(project, name)
    s = session.spec
    signal = A.complex_signal(s)
    cutoff = s.raw_time[-1]
    A.apply_phase_and_cutoff(s, A.auto_phase(signal, s.raw_time, cutoff), cutoff)
    A.fit_zero_point(s, A.symmetry_point(s.raw_time))
    A.fit_background(s, 1.0)
    s.alpha = s.alphas[26]
    s.solution = s.alpha_regularize(s.alpha)
    A.apply_stability_cutoff(s, A.stability_limits(s)[2])
    session.done = set(STEPS[:STEPS.index("stability")])
    session.complete("stability")
    if monte_carlo:
        gm = A.GaussianModels(s)
        gm.fit_to_distribution(0)
        gm.optimize_pending()
        session.complete("gaussians")
        A.run_monte_carlo(s, 0, A.default_rs_weights(s, 0), iterations=200)
        session.complete("montecarlo")


@pytest.fixture
def app(tmp_path, monkeypatch, example_bytes):
    monkeypatch.setattr(app_module, "CONFIG_PATH", tmp_path / "config.json")
    configure_panel()
    project = create_project(tmp_path, "Demo")
    for name, sample in (("first", "WT"), ("second", "WT + ligand"), ("third", "")):
        project.add_spectrum(name, *example_bytes, sample=sample)
    _process(project, "first", monte_carlo=True)
    _process(project, "second")
    a = VenisonApp()
    a.open_project(project.path)
    return a


def _legend_labels(tab):
    figs = [obj.object for card in tab.main.objects[1:] if hasattr(card, "objects")
            for obj in card.objects if hasattr(getattr(obj, "object", None), "renderers")]
    labels = []
    for fig in figs:
        for legend in fig.select(type=Legend):
            labels.append([item.label.value for item in legend.items])
    return labels


def test_tikhonov_overlays_every_ready_spectrum(app):
    tab = next(t for t in app.tabs if t.title == "Compare")
    tab.refresh()
    waves, dists = _legend_labels(tab)
    assert waves == dists == ["first (WT)", "second (WT + ligand)"]  # third isn't ready
    sidebar = " ".join(str(getattr(o, "object", "")) for o in tab.sidebar.objects)
    assert "third" in sidebar and "Stability" in sidebar


def test_gaussian_mode_only_shows_monte_carlo_results(app):
    tab = next(t for t in app.tabs if t.title == "Compare")
    tab._mode = "Gaussian"
    tab.refresh()
    waves, dists = _legend_labels(tab)
    assert waves == dists == ["first (WT)"]


def test_hiding_a_spectrum(app):
    tab = next(t for t in app.tabs if t.title == "Compare")
    tab.refresh()
    tab._on_pick([{"name": "first"}, {"name": "second"}], ["second"])
    assert _legend_labels(tab) == [["second (WT + ligand)"]] * 2
