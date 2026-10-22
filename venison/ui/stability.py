"""Stability tab: drop the unreliable long-distance tail (background supplementation)."""

from __future__ import annotations

import threading
from functools import partial

import panel as pn
from bokeh.models import BoxAnnotation, ColumnDataSource, Label

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
        # Factor the regularization system now, so the first drag tick doesn't.
        threading.Thread(target=A.prepare_regularizer, args=(spec,), daemon=True).start()

        ymax = float(spec.solution.max() * 1.1)
        p_sol = plot(DISTANCE_LABEL, "P(r)", height=300, y_range=(0, ymax))
        lo, hi = (int(v) for v in spec.validation_range)
        env_lo, env_hi = A.validation_envelope(spec, lo, hi + 1)
        p_sol.varea("x", "y1", "y2", source=ColumnDataSource(
            {"x": spec.r, "y1": env_lo, "y2": env_hi}),
            fill_color=COLORS["sky"], fill_alpha=0.3, level="underlay")
        # Zones are labeled in text; tints only reinforce the warning (faint
        # tints can't be told apart reliably with color blindness on their own).
        r0, r1 = float(spec.r[0]), float(spec.r[-1])
        for left, right, label, color, alpha in (
                (r0, width_cut, "shape reliable", None, 0),
                (width_cut, dist_cut, "mean only", "#F0E442", 0.25),
                (dist_cut, r1, "unreliable", "#D55E00", 0.18)):
            left, right = max(left, r0), min(right, r1)
            if right <= left:
                continue
            if color:
                p_sol.add_layout(BoxAnnotation(left=left, right=right, fill_color=color,
                                               fill_alpha=alpha, level="underlay"))
            p_sol.add_layout(Label(x=(left + right) / 2, y=ymax * 0.98, text=label,
                                   text_align="center", text_baseline="top",
                                   text_font_size="10.5px", text_color=COLORS["ink"],
                                   text_alpha=0.75))
        for boundary in (width_cut, dist_cut):
            if r0 < boundary < r1:
                vline(p_sol, boundary, color=COLORS["muted"], dash="dotted", width=1)
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

        main = [grid(card("Distance reliability for this trace length", p_sol),
                     card("Before and after removing the tail", p_filt)),
                card("Stable waveform and fit", p_wave)]

        dr = float(spec.r[1] - spec.r[0])
        self._slider = sidebar_widget(pn.widgets.FloatSlider(
            name="Instability point (nm)", start=float(spec.r[0]), end=float(spec.r[-1]),
            step=dr, value=float(spec.r[idx]), format="0.00"))
        # Live while dragging; saved on release.
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
        return A.renormalized_solution(self.spec, idx)

    def _apply(self, idx, save=True):
        A.apply_stability_cutoff(self.spec, idx)
        self._commit(idx, save)

    def _commit(self, idx, save):
        # Gaussian models were fit to the previous stable distribution.
        self.spec.gauss_fitted = None
        if save:
            self.complete(instability_point=float(self.spec.r[idx]))
        else:
            self.invalidate()

    def _set(self, idx):
        self._slider.value = float(self.spec.r[idx])
        self._update(save=True)

    # A drag sends ticks faster than a re-regularization (~40 ms) finishes, and
    # refitting on every tick inline queued them up so the plots trailed the
    # slider by seconds. Instead each tick only records where the slider is;
    # one worker thread refits to the latest position (skipping any it has
    # been overtaken on) and hands the result back to the event loop.
    def _update(self, save):
        spec = self.spec
        idx = A.nearest_index(spec.r, self._slider.value)
        self._cut.location = float(spec.r[idx])
        self._renorm.data = {"x": spec.r, "y": self._renormalized(idx)}
        self._target = idx
        self._save_pending = getattr(self, "_save_pending", False) or save
        if not getattr(self, "_refitting", False):
            self._refitting = True
            threading.Thread(target=self._refit_latest, args=(self.session, pn.state.curdoc),
                             daemon=True).start()

    def _refit_latest(self, session, doc):
        spec = session.spec
        while True:
            idx = self._target
            try:
                result = A.stability_cutoff(spec, idx)
            except Exception as exc:
                result = exc
            if idx == self._target:
                break
        finish = partial(self._show_refit, session, idx, result)
        if doc is not None:
            doc.add_next_tick_callback(finish)
        else:  # no browser session (tests)
            finish()

    def _show_refit(self, session, idx, result):
        self._refitting = False
        if session is not self.session:
            # A different spectrum was opened while this refit ran; it was
            # computed for the old one, so drop it rather than save it there.
            self._save_pending = False
            return
        spec = session.spec
        if isinstance(result, Exception):
            if self._save_pending:  # mid-drag failures just skip a frame
                self.app.notify(f"Stability refit failed: {result}", "error")
            self._save_pending = False
            return
        for key, value in result.items():
            setattr(spec, key, value)
        self._stable.data = {"x": spec.r, "y": spec.stable_solution}
        self._swave.data = {"x": spec.t, "y": spec.stable_waveform}
        self._sfit.data = {"x": spec.t, "y": spec.stable_fit}
        if idx != self._target:  # the slider moved on while this was shown
            self._update(save=False)
            return
        save, self._save_pending = self._save_pending, False
        self._commit(idx, save)
