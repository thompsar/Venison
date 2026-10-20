"""Validation tab: background start validation by reconstruction RMSD."""

from __future__ import annotations

import numpy as np
import panel as pn
from bokeh.models import ColumnDataSource

from venison.core import analysis as A
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid,
                               hint, line, plot, points, section, sidebar_widget, vline,
                               zero_line)

_ABOUT = (
    "Tries many background start positions. For each one the background is refit, the "
    "waveform regularized at your alpha, and the raw echo rebuilt; the RMSD says how well "
    "that rebuild matches the data.<br><br>"
    "Starts that are too early fold dipolar signal into the background. Very late starts "
    "can score well only because the error is soaked up as long-distance P(r). The "
    "suggested range is the bottom of the first RMSD well plus its statistically "
    "equivalent neighbours; the blue band on the distance plot is the spread across it."
)


class ValidationTab(Tab):
    title = "Validation"
    icon = "validation"
    step = "validation"
    requires = "lcurve"

    def build(self):
        spec = self.spec
        side = [hint(_ABOUT)]
        done = self.session.is_done(self.step) and hasattr(spec, "validation_rmsd")
        main = self._build_plots() if done else [
            card("", hint("Set the sweep in the sidebar and click <b>Run validation</b>."))]
        if done:
            side += self._selection_controls()
        return main, side + self._run_controls()

    # ── Sidebar ──────────────────────────────────────────────────────────────

    def _run_controls(self):
        spec = self.spec
        last = dict(A.VALIDATION_DEFAULTS, **getattr(spec, "validation_settings", {}) or {})
        self._earliest = sidebar_widget(pn.widgets.FloatInput(
            name="Earliest start (µs)", value=float(last["earliest_start"]), step=0.05, start=0))
        self._shortest = sidebar_widget(pn.widgets.FloatInput(
            name="Shortest background fit (µs)", value=float(last["shortest_fit"]), step=0.1,
            start=0.05))
        self._trials = sidebar_widget(pn.widgets.IntInput(
            name="Trials", value=int(last["n_trials"]), start=10, step=10))
        out = [section("Sweep"), self._earliest, self._shortest, self._trials,
               hint("The latest trial starts this far before the cutoff, so every background "
                    "is fit to at least this much of the trace.")]

        self._hold_d = None
        if getattr(spec, "bgmodel_name", "3D") == "Fractal":
            d = A.background_tab_fractal_d(spec)
            if d is None:
                out.append(hint("<b>Note:</b> the Fractal fit on the Background tab failed and "
                                "fell back to 3D, so there is no dimensionality to hold."))
            else:
                self._hold_d = sidebar_widget(pn.widgets.Checkbox(
                    name=f"Hold Fractal dimension at d = {d:.2f}",
                    value=bool(last["hold_fractal_d"])))
                out += [self._hold_d,
                        hint("Left free, d swings from trial to trial and makes the RMSD "
                             "curve jump between bands.")]
                if min(abs(d - 1), abs(d - 5)) < 1e-3:
                    out.append(hint(f"<b>Warning:</b> d = {d:.2f} is at its limit (1–5); the "
                                    "Fractal model may not suit this trace."))
        run = button("Run validation", "primary")
        run.on_click(self._on_run)
        return out + [run]

    def _selection_controls(self):
        spec = self.spec
        v = A.validation_of(spec)
        x = v.fit_points
        dx = float(x[1] - x[0]) if len(x) > 1 else 0.01
        lo, hi = (int(i) for i in getattr(spec, "validation_range", v.well))
        sel = int(getattr(spec, "validation_index", v.min_index))
        self._range = sidebar_widget(pn.widgets.RangeSlider(
            name="Range (µs)", start=float(x[0]), end=float(x[-1]), step=dx,
            value=(float(x[lo]), float(x[hi])), format="0.00"))
        self._range.param.watch(lambda e: self._on_range(save=False), "value")
        self._range.param.watch(lambda e: self._on_range(save=True), "value_throttled")
        self._pick = sidebar_widget(pn.widgets.FloatSlider(
            name="Selected start (µs)", start=float(x[0]), end=float(x[-1]), step=dx,
            value=float(x[sel]), format="0.00"))
        self._pick.param.watch(lambda e: self._show(self._index(e.new)), "value")
        self._pick.param.watch(lambda e: self._select(self._index(e.new)), "value_throttled")
        to_min = button("Select range minimum")
        to_min.on_click(lambda e: self._select(self._range_min()))
        reset = button("Reset to suggestion")
        reset.on_click(self._on_reset)
        out = [section("Selection"), self._range, self._pick, pn.Row(to_min, reset, margin=0)]
        notes = []
        settings = getattr(spec, "validation_settings", {}) or {}
        if settings.get("fallbacks"):
            notes.append(f"<b>{settings['fallbacks']} of {len(x)}</b> trials couldn't fit a "
                         "Fractal background and used 3D instead.")
        if settings.get("held_d") is not None:
            notes.append(f"Fractal dimension held at d = {settings['held_d']:.2f}.")
        return out + [hint(n) for n in notes]

    # ── Plots ────────────────────────────────────────────────────────────────

    def _build_plots(self):
        spec = self.spec
        v = A.validation_of(spec)
        lo, hi = (int(i) for i in getattr(spec, "validation_range", v.well))
        sel = int(getattr(spec, "validation_index", v.min_index))
        x = v.fit_points

        # Log axis, everything visible: the curve spans orders of magnitude and
        # the context outside the range is what makes the well readable.
        p_rmsd = plot("Background fit start (µs)", "Reconstruction RMSD", height=300,
                      y_axis_type="log")
        points(p_rmsd, x, v.rmsd, color=COLORS["muted"], size=4)
        self._well = ColumnDataSource({"x": x[lo:hi + 1], "y": v.rmsd[lo:hi + 1]})
        p_rmsd.line("x", "y", source=self._well, line_color=COLORS["alt"], line_width=4,
                    line_alpha=0.8)
        m = A.validation_range_min(v, lo, hi + 1)
        self._min = points(p_rmsd, x[m], v.rmsd[m], color="#C0392B", size=12,
                           fill_alpha=0, line_width=2)
        self._sel = points(p_rmsd, x[sel], v.rmsd[sel], color=COLORS["ink"], size=7)
        if hasattr(spec, "bg_tab_start"):
            vline(p_rmsd, spec.bg_tab_start, color=COLORS["fit"])

        p_bg = plot(TIME_LABEL, height=300)
        line(p_bg, spec.raw_time, spec.real, COLORS["muted"])
        line(p_bg, spec.raw_time, spec.background_before_validation, COLORS["fit"],
             line_dash="dashed", legend_label="Background tab")
        self._bg = line(p_bg, spec.raw_time, spec.backgrounds[sel], COLORS["alt"], width=2,
                        legend_label="Selected")
        p_bg.legend.label_text_font_size = "11px"

        p_wave = plot(TIME_LABEL, "Form factor", height=300)
        zero_line(p_wave)
        self._wave = line(p_wave, spec.t, spec.waveforms[sel], COLORS["muted"])
        self._fit = line(p_wave, spec.t, spec.validation_tikhonovfits[sel], COLORS["fit"], width=2)

        p_dist = plot(DISTANCE_LABEL, "P(r)", height=300)
        env_lo, env_hi = A.validation_envelope(spec, lo, hi + 1)
        self._env = ColumnDataSource({"x": spec.r, "y1": env_lo, "y2": env_hi})
        p_dist.varea("x", "y1", "y2", source=self._env, fill_color=COLORS["sky"],
                     fill_alpha=0.35, level="underlay")
        self._dist = line(p_dist, spec.r, spec.validation_solutions[sel], COLORS["data"], width=2)

        return [grid(card("RMSD by background start (green: range, red: its minimum, "
                          "dashed: Background tab)", p_rmsd),
                     card("Backgrounds", p_bg),
                     card("Selected waveform and fit", p_wave),
                     card("Distance distribution (band: spread over the range)", p_dist))]

    # ── Selection ────────────────────────────────────────────────────────────

    def _index(self, start: float) -> int:
        return A.nearest_index(self.spec.background_fit_points, start)

    def _range_indices(self) -> tuple[int, int]:
        lo, hi = (self._index(v) for v in self._range.value)
        return lo, hi

    def _range_min(self) -> int:
        lo, hi = self._range_indices()
        return A.validation_range_min(A.validation_of(self.spec), lo, hi + 1)

    def _show(self, sel):
        spec = self.spec
        v = A.validation_of(spec)
        self._sel.data = {"x": [v.fit_points[sel]], "y": [v.rmsd[sel]]}
        self._bg.data = {"x": spec.raw_time, "y": spec.backgrounds[sel]}
        self._wave.data = {"x": spec.t, "y": spec.waveforms[sel]}
        self._fit.data = {"x": spec.t, "y": spec.validation_tikhonovfits[sel]}
        self._dist.data = {"x": spec.r, "y": spec.validation_solutions[sel]}

    def _select(self, sel):
        sel = int(sel)
        x = self.spec.background_fit_points
        if self._index(self._pick.value) != sel:
            self._pick.value = float(x[sel])
        self._show(sel)
        lo, hi = self._range_indices()
        A.select_validated_background(self.spec, sel, (lo, hi))
        self.complete(background_start=float(self.spec.background_start),
                      validation_range=[float(x[lo]), float(x[hi])])

    def _on_range(self, save):
        spec = self.spec
        v = A.validation_of(spec)
        lo, hi = self._range_indices()
        self._well.data = {"x": v.fit_points[lo:hi + 1], "y": v.rmsd[lo:hi + 1]}
        m = self._range_min()
        self._min.data = {"x": [v.fit_points[m]], "y": [v.rmsd[m]]}
        env_lo, env_hi = A.validation_envelope(spec, lo, hi + 1)
        self._env.data = {"x": spec.r, "y1": env_lo, "y2": env_hi}
        if save:
            self._select(m)
        else:
            self._show(m)

    def _on_reset(self, _event):
        v = A.validation_of(self.spec)
        x = v.fit_points
        self._range.value = (float(x[v.well[0]]), float(x[v.well[1]]))
        self._on_range(save=True)
        self._select(v.min_index)

    async def _on_run(self, _event):
        spec = self.spec
        kwargs = {"earliest_start": float(self._earliest.value),
                  "shortest_fit": float(self._shortest.value),
                  "n_trials": int(self._trials.value),
                  "hold_fractal_d": bool(self._hold_d.value) if self._hold_d is not None else False}
        if kwargs["earliest_start"] >= float(spec.cutoff) - kwargs["shortest_fit"]:
            self.app.notify("The earliest start must leave room for the shortest background "
                            "fit before the cutoff.", "warning")
            return

        def work():
            # Validate against the Background tab's choice, not a previous
            # validation's selection, and keep it for comparison in the plots.
            spec.background_start = getattr(spec, "bg_tab_start", spec.background_start)
            v = A.validate_background(spec, **kwargs)
            spec.background_before_validation = np.array(spec.background)
            return v

        def then(v):
            spec.validation_range = np.array(v.well)
            spec.validation_index = v.min_index
            A.select_validated_background(spec, v.min_index, v.well)
            self.complete(**kwargs, background_start=float(spec.background_start),
                          suggested_range=[float(v.fit_points[i]) for i in v.well])
            self.refresh()

        await self.app.run_task("Validating background", work, then)
