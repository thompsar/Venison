"""Background tab: zero time, V(t=0) normalization, and background fit."""

from __future__ import annotations

import numpy as np
import panel as pn

from venison.core import analysis as A
from venison.ui.common import (COLORS, TIME_LABEL, Tab, card, grid, hint, line, plot, points,
                               section, sidebar_widget, vline, zero_line)


_EQUATIONS = (
    "<b>3D:</b> B(t) = b·exp(−t / c)<br>"
    "<b>Fractal:</b> B(t) = b·exp(−t<sup>d/3</sup> / c)<br>"
    "d is the dimensionality of the spin distribution: 3 for spins spread evenly "
    "through solution (identical to the 3D model), lower for e.g. membranes. "
    "c is the decay constant (µs for 3D)."
)


class BackgroundTab(Tab):
    title = "Background"
    icon = "background"
    step = "background"
    requires = "kernel"

    def build(self):
        spec = self.spec
        lo, hi = A.zero_window_bounds(spec)
        options = A.background_start_options(spec)
        if not self.session.is_done(self.step):
            # First visit (or an earlier step changed): fit with the current
            # settings, falling back to the notebook defaults.
            spec.zero_window = int(np.clip(getattr(spec, "zero_window",
                                                   A.symmetry_point(spec.raw_time)), lo, hi))
            spec.normalizer = str(getattr(spec, "normalizer", A.NORMALIZERS[0]))
            start = options[A.nearest_index(options, spec.background_start)]
            self._fit(start, str(getattr(spec, "bgmodel_name", "3D")))

        w = int(spec.zero_window)
        p_zero = plot(TIME_LABEL, height=280)
        self._zero_src = line(p_zero, spec.raw_time[1:w], spec.real[1:w], width=1.5)
        self._zero_fit = line(p_zero, spec.raw_time[1:w],
                              A.gaussian(spec.raw_time[1:w], *spec.zero_fit), COLORS["fit"], width=2)
        self._zero_line = vline(p_zero, spec.zeropoint, color=COLORS["alt"])
        self._zero_pt = points(p_zero, spec.zeropoint, spec.zeroamp, color="#C0392B", size=8)

        p_bg = plot(TIME_LABEL, height=280)
        line(p_bg, spec.raw_time, spec.real)
        # 3D fit over the same range, shown for comparison when Fractal is selected.
        self._bg3d_src = line(p_bg, spec.raw_time, spec.raw_time * np.nan, COLORS["muted"],
                              width=1.5, line_dash="dashed")
        self._bg_src = line(p_bg, spec.raw_time, spec.background, COLORS["fit"], width=2)
        vline(p_bg, spec.cutoff, color="#C0392B", dash="solid")
        self._bg_line = vline(p_bg, spec.background_start, color=COLORS["data"])

        p_wave = plot(TIME_LABEL, "Form factor", height=300)
        zero_line(p_wave)
        self._wave_src = line(p_wave, spec.t, spec.waveform)

        main = [grid(card("Zero time (Gaussian fit to the echo maximum)", p_zero),
                     card("Background fit (start in blue, cutoff in red; dashed: 3D fit "
                          "for comparison)", p_bg)),
                card("Background-corrected waveform", p_wave)]

        self._window = sidebar_widget(pn.widgets.IntSlider(
            name="Fit window (points)", start=lo, end=hi, value=w))
        self._norm = sidebar_widget(pn.widgets.Select(
            name="V(t=0) from", options=list(A.NORMALIZERS), value=spec.normalizer))
        self._model = sidebar_widget(pn.widgets.Select(
            name="Background model", options=list(A.BG_MODELS),
            value=str(getattr(spec, "bgmodel_name", "3D"))))
        dt = float(spec.raw_time[1] - spec.raw_time[0])
        self._start = sidebar_widget(pn.widgets.FloatSlider(
            name="Background start (µs)", start=float(options[0]), end=float(options[-1]),
            step=dt, value=float(options[A.nearest_index(options, spec.background_start)]),
            format="0.000"))
        # Sliders refit live while dragging, like the notebook; the result is
        # saved (and later steps cleared) once, on release.
        for widget in (self._window, self._start):
            widget.param.watch(lambda e: self._on_change(save=False), "value")
            widget.param.watch(lambda e: self._on_change(save=True), "value_throttled")
        for widget in (self._norm, self._model):
            widget.param.watch(lambda e: self._on_change(save=True), "value")

        self._fit_info = pn.pane.HTML("", margin=(2, 8), sizing_mode="stretch_width")
        side = [section("Zero time"), self._window, self._norm,
                hint("The Gaussian fit window starts at the second point."),
                section("Background"), self._model, self._start,
                hint("Fit from the start to the cutoff. Refine it on the Validation tab."),
                hint(_EQUATIONS), self._fit_info]
        if self.session.is_done("validation"):
            side.append(hint("<b>Note:</b> the plots show the background chosen on the "
                             "Validation tab; the parameters above are for the start set "
                             "here. Changing anything here clears the validation choice."))
        self._show_fit_info()
        return main, side

    def _show_fit_info(self):
        spec = self.spec
        options = A.background_start_options(spec)
        start = options[A.nearest_index(options, self._start.value)]
        try:
            s = A.summarize_background(spec, start, self._model.value)
        except Exception:  # fit failure mid-drag: keep the last numbers
            return
        cols = [s] + ([s.alt] if s.alt else [])
        head = "".join(f"<th style='text-align:right;padding:0 0 0 10px'>{c.model}</th>"
                       for c in cols)

        def row(label, values):
            cells = "".join(f"<td style='text-align:right;padding:0 0 0 10px;"
                            f"font-variant-numeric:tabular-nums'>{v}</td>" for v in values)
            return f"<tr><td style='color:#7A706A'>{label}</td>{cells}</tr>"

        rows = [row("b", [f"{c.b:.3f}" for c in cols]),
                row("c", [f"{c.c:.3g}" for c in cols]),
                row("d", [f"{c.d:.2f}" for c in cols]),
                row("fit RMSD", [f"{c.rmsd:.2e}" for c in cols])]
        warn = ""
        if s.model == "Fractal" and min(abs(s.d - 1), abs(s.d - 5)) < 0.1:
            warn = (f"<div class='vn-hint' style='padding:4px 0 0'><b>Warning:</b> d is at or near "
                    f"its limit (1–5); the Fractal model may not suit this range.</div>")
        elif self._model.value == "Fractal" and s.model == "3D":
            warn = ("<div class='vn-hint' style='padding:4px 0 0'><b>Note:</b> the Fractal "
                    "fit failed for this range, so 3D is used.</div>")
        self._fit_info.object = (
            "<table style='font-size:12.5px;border-collapse:collapse;width:100%'>"
            f"<tr><th></th>{head}</tr>{''.join(rows)}</table>{warn}")
        self._bg3d_src.data = {"x": spec.raw_time,
                               "y": s.alt_curve if s.alt_curve is not None
                               else spec.raw_time * np.nan}

    def _fit(self, start, model, save=True):
        spec = self.spec
        spec.zero_fit = A.fit_zero_point(spec, int(spec.zero_window), spec.normalizer)
        A.fit_background(spec, start, model)
        spec.bg_tab_start = float(start)
        if save:
            self.complete(zero_window=int(spec.zero_window), normalizer=spec.normalizer,
                          background_model=model, background_start=float(start),
                          zeropoint=float(spec.zeropoint))
        else:
            self.invalidate()  # clears later steps' dots now; only writes when they change

    def _on_change(self, save):
        spec = self.spec
        spec.zero_window = int(self._window.value)
        spec.normalizer = self._norm.value
        options = A.background_start_options(spec)
        start = options[A.nearest_index(options, self._start.value)]
        try:
            self._fit(start, self._model.value, save=save)
        except Exception as exc:  # e.g. curve_fit not converging for a tiny window
            if save:  # mid-drag failures just skip a frame
                self.app.notify(f"Fit failed: {exc}", "error")
            return
        w = int(spec.zero_window)
        x = spec.raw_time[1:w]
        self._zero_src.data = {"x": x, "y": spec.real[1:w]}
        self._zero_fit.data = {"x": x, "y": A.gaussian(x, *spec.zero_fit)}
        self._zero_line.location = float(spec.zeropoint)
        self._zero_pt.data = {"x": [spec.zeropoint], "y": [spec.zeroamp]}
        self._bg_src.data = {"x": spec.raw_time, "y": spec.background}
        self._bg_line.location = float(spec.background_start)
        self._wave_src.data = {"x": spec.t, "y": spec.waveform}
        self._show_fit_info()
