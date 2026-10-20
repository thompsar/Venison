"""Validation tab: background start validation by reconstruction RMSD."""

from __future__ import annotations

import numpy as np
import panel as pn
from bokeh.models import ColumnDataSource

from venison.core import analysis as A
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid,
                               hint, line, plot, points, section, sidebar_widget, zero_line)


class ValidationTab(Tab):
    title = "Validation"
    icon = "validation"
    step = "validation"
    requires = "lcurve"

    def build(self):
        self._offset = sidebar_widget(pn.widgets.IntInput(
            name="End offset (points before cutoff)", value=100, start=1))
        run = button("Run validation", "primary")
        run.on_click(self._on_run)
        run_side = [section("Run"), self._offset, run,
                    hint("Refits the background from ~100 start positions and regularizes each "
                         "waveform at the selected alpha.")]
        spec = self.spec
        if not self.session.is_done(self.step) or not hasattr(spec, "validation_rmsd"):
            return ([card("", hint("Click <b>Run validation</b> in the sidebar to start."))],
                    run_side)

        v = A.validation_of(spec)
        lo, hi = (int(x) for x in getattr(spec, "validation_range", v.well))
        sel = int(getattr(spec, "validation_index", v.min_index))
        x = v.fit_points

        p_rmsd = plot("Background fit start (µs)", "Reconstruction RMSD", height=300)
        points(p_rmsd, x, v.rmsd, color=COLORS["muted"], size=5)
        self._well = ColumnDataSource({"x": x[lo:hi + 1], "y": v.rmsd[lo:hi + 1]})
        p_rmsd.line("x", "y", source=self._well, line_color=COLORS["ink"], line_width=3)
        m = A.validation_range_min(v, lo, hi + 1)
        self._min = points(p_rmsd, x[m], v.rmsd[m], color="#C0392B", size=11,
                           fill_alpha=0, line_width=2)
        self._sel = points(p_rmsd, x[sel], v.rmsd[sel], color=COLORS["alt"], size=8)
        self._p_rmsd = p_rmsd
        self._zoom(lo, hi)

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
                     fill_alpha=0.3, level="underlay")
        self._dist = line(p_dist, spec.r, spec.validation_solutions[sel], COLORS["data"], width=2)

        main = [grid(card("RMSD by background start (well in bold, minimum circled)", p_rmsd),
                     card("Backgrounds", p_bg),
                     card("Selected waveform and fit", p_wave),
                     card("Distance distribution (envelope over the range)", p_dist))]

        n = len(v.rmsd)
        self._range = sidebar_widget(pn.widgets.IntRangeSlider(
            name="Backgrounds in range", start=0, end=n - 1, value=(lo, hi)))
        self._range.param.watch(self._on_range, "value_throttled")
        self._pick = sidebar_widget(pn.widgets.IntSlider(
            name="Selected background", start=lo, end=max(hi, lo + 1), value=sel))
        self._pick.param.watch(self._on_pick_drag, "value")
        self._pick.param.watch(self._on_pick, "value_throttled")
        to_min = button("Select range minimum")
        to_min.on_click(self._on_to_min)
        self._info = pn.pane.HTML(self._info_html(sel), margin=0)
        side = [section("Selection"), self._range, self._pick, to_min, self._info,
                hint("The range starts as the F-test well around the RMSD minimum.")] + run_side
        return main, side

    def _info_html(self, sel):
        start = self.spec.background_fit_points[sel]
        return f'<div class="vn-kv"><b>Start</b>{start:.3f} µs</div>'

    def _zoom(self, lo, hi):
        rmsd = A.validation_of(self.spec).rmsd[lo:hi + 1]
        pad = (rmsd.max() - rmsd.min()) * 0.1 or rmsd.max() * 0.02
        self._p_rmsd.y_range.start = float(rmsd.min() - pad)
        self._p_rmsd.y_range.end = float(rmsd.max() + pad)

    def _show(self, sel):
        spec = self.spec
        v = A.validation_of(spec)
        self._sel.data = {"x": [v.fit_points[sel]], "y": [v.rmsd[sel]]}
        self._bg.data = {"x": spec.raw_time, "y": spec.backgrounds[sel]}
        self._wave.data = {"x": spec.t, "y": spec.waveforms[sel]}
        self._fit.data = {"x": spec.t, "y": spec.validation_tikhonovfits[sel]}
        self._dist.data = {"x": spec.r, "y": spec.validation_solutions[sel]}
        self._info.object = self._info_html(sel)

    def _select(self, sel):
        sel = int(sel)
        if self._pick.value != sel:
            self._pick.value = sel
        self._show(sel)
        A.select_validated_background(self.spec, sel, self._range.value)
        self.complete(background_start=float(self.spec.background_start),
                      validation_range=list(self._range.value))

    def _on_to_min(self, _event):
        lo, hi = self._range.value
        self._select(A.validation_range_min(A.validation_of(self.spec), lo, hi + 1))

    def _on_pick_drag(self, event):
        self._show(int(event.new))

    def _on_pick(self, event):
        self._select(event.new)

    def _on_range(self, event):
        lo, hi = event.new
        spec = self.spec
        v = A.validation_of(spec)
        self._well.data = {"x": v.fit_points[lo:hi + 1], "y": v.rmsd[lo:hi + 1]}
        m = A.validation_range_min(v, lo, hi + 1)
        self._min.data = {"x": [v.fit_points[m]], "y": [v.rmsd[m]]}
        env_lo, env_hi = A.validation_envelope(spec, lo, hi + 1)
        self._env.data = {"x": spec.r, "y1": env_lo, "y2": env_hi}
        self._zoom(lo, hi)
        self._pick.start, self._pick.end = lo, max(hi, lo + 1)
        self._select(m)

    async def _on_run(self, _event):
        spec = self.spec
        offset = int(self._offset.value)

        def work():
            # Validate against the Background tab's choice, not a previous
            # validation's selection, and keep it for comparison in the plots.
            spec.background_start = getattr(spec, "bg_tab_start", spec.background_start)
            v = A.validate_background(spec, end_offset=offset)
            spec.background_before_validation = np.array(spec.background)
            return v

        def then(v):
            spec.validation_range = np.array(v.well)
            A.select_validated_background(spec, v.min_index, v.well)
            self.complete(end_offset=offset, background_start=float(spec.background_start),
                          validation_range=list(v.well))
            self.refresh()

        await self.app.run_task("Validating background", work, then)
