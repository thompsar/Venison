"""Compare tab: overlay several spectra of the project, waveforms and distances."""

from __future__ import annotations

import html
from dataclasses import dataclass

import numpy as np
import panel as pn

from venison.core import analysis as A
from venison.core.session import STEP_LABELS, SpectrumSession
from venison.ui.common import (DISTANCE_LABEL, TIME_LABEL, Tab, card, empty_state, hint, line,
                               page_header, plot, section, sidebar_widget, zero_line)

# One color per spectrum; Okabe–Ito, ordered so neighbours stay distinct with
# common color blindness.
SERIES_COLORS = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#56B4E9", "#E69F00", "#26211E"]

MODES = {"Tikhonov": "stability", "Gaussian": "montecarlo"}  # mode -> step it needs


@dataclass
class Curves:
    t: np.ndarray
    data: np.ndarray       # stable waveform
    fit: np.ndarray
    r: np.ndarray
    distribution: np.ndarray   # area-normalized over r


def curves(spec, mode: str) -> Curves:
    if mode == "Tikhonov":
        sol = spec.stable_solution
        return Curves(spec.t, spec.stable_waveform, spec.stable_fit, spec.r,
                      sol / np.trapezoid(sol, spec.r))
    if getattr(spec, "best_model_fit", None) is None:  # older analyses: compute it now
        A.best_monte_carlo_model(spec)
    return Curves(spec.t, spec.stable_waveform, spec.best_model_fit, spec.r,
                  spec.best_model_solution)


def label_for(info: dict) -> str:
    extra = ", ".join(v for v in (info.get("sample"), info.get("condition")) if v)
    return f"{info['name']} ({extra})" if extra else info["name"]


class CompareTab(Tab):
    title = "Compare"
    icon = "compare"
    needs_spectrum = False

    def __init__(self, app):
        super().__init__(app)
        self._mode = "Tikhonov"
        self._hidden: set[str] = set()
        # (name, mode) -> (workspace mtime, Curves). Only the plotted arrays are
        # kept (a few kB each), not the whole saved analysis (10-20 MB each).
        self._cache: dict[tuple[str, str], tuple[float, Curves]] = {}

    # ── Data ─────────────────────────────────────────────────────────────────

    def _curves_for(self, name: str, mode: str) -> Curves:
        """Plotted curves for a spectrum; the open one comes straight from memory."""
        session = self.app.session
        if session is not None and session.name == name:
            return curves(session.spec, mode)
        npz, _ = self.app.project.workspace_paths(name)
        mtime = npz.stat().st_mtime
        cached = self._cache.get((name, mode))
        if cached is None or cached[0] != mtime:
            spec = SpectrumSession(self.app.project, name).spec  # loaded, used, dropped
            cached = self._cache[(name, mode)] = (mtime, curves(spec, mode))
        return cached[1]

    # ── View ─────────────────────────────────────────────────────────────────

    def build(self):
        project = self.app.project
        if project is None:
            return [empty_state("No project open", "Open a project in the Project tab.")], []
        spectra = project.spectra()
        needed = MODES[self._mode]
        ready = [s for s in spectra if needed in s["steps"]]
        waiting = [s for s in spectra if needed not in s["steps"]]
        colors = {s["name"]: SERIES_COLORS[i % len(SERIES_COLORS)]
                  for i, s in enumerate(spectra)}  # stable per spectrum across modes

        mode = sidebar_widget(pn.widgets.RadioButtonGroup(
            options=list(MODES), value=self._mode, button_type="default"))
        mode.param.watch(self._on_mode, "value")
        side = [section("Show"), mode,
                hint("<b>Tikhonov:</b> stable waveform, Tikhonov fit and stable P(r). "
                     "<b>Gaussian:</b> best Monte Carlo model fit and its P(r).")]
        if ready:
            pick = sidebar_widget(pn.widgets.CheckBoxGroup(
                options={label_for(s): s["name"] for s in ready},
                value=[s["name"] for s in ready if s["name"] not in self._hidden],
                inline=False))
            pick.param.watch(lambda e: self._on_pick(ready, e.new), "value")
            side += [section("Spectra"), pick]
        if waiting:
            items = "".join(f"<li>{html.escape(s['name'])}</li>" for s in waiting)
            side += [section("Not ready"),
                     hint(f"These need <b>{STEP_LABELS[needed]}</b> first:"
                          f"<ul style='margin:4px 0 0;padding-left:18px'>{items}</ul>")]

        header = page_header("Compare spectra",
                             f"{len(ready)} of {len(spectra)} spectra in this project "
                             f"have a {self._mode} result.")
        shown = [s for s in ready if s["name"] not in self._hidden]
        if not shown:
            body = "Analyze at least one spectrum through the " \
                   f"{STEP_LABELS[needed]} step." if not ready else "All spectra are hidden."
            return [header, card("", hint(body))], side

        p_wave = plot(TIME_LABEL, "Form factor", height=360)
        zero_line(p_wave)
        p_dist = plot(DISTANCE_LABEL, "P(r), area-normalized", height=360)
        failed = []
        for s in shown:
            try:
                c = self._curves_for(s["name"], self._mode)
            except Exception as exc:  # unreadable or incomplete saved analysis
                failed.append(f"{s['name']}: {exc}")
                continue
            color, label = colors[s["name"]], label_for(s)
            line(p_wave, c.t, c.data, color, width=1, line_alpha=0.3)
            line(p_wave, c.t, c.fit, color, width=2.2, legend_label=label)
            line(p_dist, c.r, c.distribution, color, width=2.2, legend_label=label)
        for p in (p_wave, p_dist):
            if p.legend:
                p.legend.label_text_font_size = "11px"
                p.legend.click_policy = "hide"
        main = [header,
                card(f"Waveforms (faint) and {self._mode} fits", p_wave),
                card("Distance distributions", p_dist)]
        if failed:
            main.insert(1, hint("Couldn't load: " + "; ".join(html.escape(f) for f in failed)))
        side.append(hint("Click a legend entry to hide that spectrum's curve."))
        return main, side

    # ── Callbacks ────────────────────────────────────────────────────────────

    def _on_mode(self, event):
        self._mode = event.new
        self.refresh()

    def _on_pick(self, ready, chosen):
        self._hidden = {s["name"] for s in ready} - set(chosen)
        self.refresh()
