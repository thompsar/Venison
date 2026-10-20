"""Stability tab: drop the unreliable long-distance tail (background supplementation)."""

from __future__ import annotations

import numpy as np
import panel as pn
from bokeh.models import BoxAnnotation, ColumnDataSource

from venison.core import analysis as A
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid, hint,
                               line, plot, section, sidebar_widget, vline, zero_line)


class StabilityTab(Tab):
    title = "Stability"
    icon = "stability"
    step = "stability"
    requires = "validation"

    def build(self):
        spec = self.spec
        width_cut, dist_cut, default_idx = A.stability_limits(spec)
        if not self.session.is_done(self.step):
            self._apply(int(getattr(spec, "instability_index", default_idx)))
        idx = int(spec.instability_index)
        renormalized = self._renormalized(idx)

        ymax = float(spec.solution.max() * 1.1)
        p_sol = plot(DISTANCE_LABEL, "P(r)", height=300, y_range=(0, ymax))
        lo, hi = (int(v) for v in spec.validation_range)
        env_lo, env_hi = A.validation_envelope(spec, lo, hi + 1)
        p_sol.varea("x", "y1", "y2", source=ColumnDataSource(
            {"x": spec.r, "y1": env_lo, "y2": env_hi}),
            fill_color=COLORS["sky"], fill_alpha=0.3, level="underlay")
        for left, right, color in ((None, width_cut, "#3D8B5A"), (width_cut, dist_cut, "#E6B400"),
                                   (dist_cut, None, "#C0392B")):
            p_sol.add_layout(BoxAnnotation(left=left, right=right, fill_color=color,
                                           fill_alpha=0.08, level="underlay"))
        line(p_sol, spec.r, spec.solution, COLORS["data"], width=2)
        self._cut = vline(p_sol, spec.r[idx])

        p_filt = plot(DISTANCE_LABEL, "P(r), normalized", height=300)
        self._renorm = line(p_filt, spec.r, renormalized, COLORS["muted"],
                            legend_label="Before")
        self._stable = line(p_filt, spec.r, spec.stable_solution, COLORS["data"], width=2,
                            legend_label="Stable")
        p_filt.legend.label_text_font_size = "11px"

        p_wave = plot(TIME_LABEL, "Form factor", height=280)
        zero_line(p_wave)
        self._swave = line(p_wave, spec.t, spec.stable_waveform, COLORS["muted"])
        self._sfit = line(p_wave, spec.t, spec.stable_fit, COLORS["fit"], width=2)

        main = [grid(card("Reliability: shape (green), mean only (yellow), none (red)", p_sol),
                     card("Before and after removing the tail", p_filt)),
                card("Stable waveform and fit", p_wave)]

        dr = float(spec.r[1] - spec.r[0])
        self._slider = sidebar_widget(pn.widgets.FloatSlider(
            name="Instability point (nm)", start=float(spec.r[0]), end=float(spec.r[-1]),
            step=dr, value=float(spec.r[idx]), format="0.00"))
        # Live while dragging (~0.1 s per re-regularization); saved on release.
        self._slider.param.watch(lambda e: self._update(save=False), "value")
        self._slider.param.watch(lambda e: self._update(save=True), "value_throttled")
        reset = button("Reset to default")
        reset.on_click(lambda e: self._set(default_idx))
        side = [section("Cutoff"), self._slider, reset,
                hint(f"For this trace length, distribution shape is reliable to "
                     f"{width_cut:.2f} nm and mean distance to {dist_cut:.2f} nm "
                     "(Jeschke 2012). The default cutoff sits midway.")]
        return main, side

    def _renormalized(self, idx):
        mask = np.zeros_like(self.spec.r)
        mask[idx:] = 1
        sol = self.spec.solution
        return sol / np.trapezoid(sol - sol * mask)

    def _apply(self, idx, save=True):
        A.apply_stability_cutoff(self.spec, idx)
        # Gaussian models were fit to the previous stable distribution.
        self.spec.gauss_fitted = None
        if save:
            self.complete(instability_point=float(self.spec.r[idx]))
        else:
            self.invalidate()

    def _set(self, idx):
        self._slider.value = float(self.spec.r[idx])
        self._update(save=True)

    def _update(self, save):
        spec = self.spec
        idx = A.nearest_index(spec.r, self._slider.value)
        self._cut.location = float(spec.r[idx])
        self._apply(idx, save=save)
        self._renorm.data = {"x": spec.r, "y": self._renormalized(idx)}
        self._stable.data = {"x": spec.r, "y": spec.stable_solution}
        self._swave.data = {"x": spec.t, "y": spec.stable_waveform}
        self._sfit.data = {"x": spec.t, "y": spec.stable_fit}
