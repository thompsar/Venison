"""Tikhonov tab: L-curve and regularization parameter choice."""

from __future__ import annotations

import panel as pn

from venison.core import analysis as A
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid,
                               hint, line, plot, points, section, sidebar_widget, zero_line)


class TikhonovTab(Tab):
    title = "Tikhonov"
    icon = "lcurve"
    step = "lcurve"
    requires = "background"

    def build(self):
        run = button("Compute L-curve", "primary")
        run.on_click(self._on_compute)
        if not self.session.is_done(self.step):
            return ([card("", hint(
                "Regularizes the waveform at each of the 51 alpha values "
                f"(10<sup>-4</sup> to 10<sup>6</sup>) and plots the L-curve. "
                "Click <b>Compute L-curve</b> in the sidebar to start."))],
                [section("L-curve"), run])

        spec = self.spec
        lc = A.lcurve_of(spec)
        idx = int(getattr(spec, "alpha_index", lc.loocv_index))

        p_l = plot("Residual norm", "Solution norm", height=300,
                   x_axis_type="log", y_axis_type="log")
        points(p_l, lc.residual_norm, lc.solution_norm, color=COLORS["data"], size=5, alpha=0.7)
        points(p_l, lc.residual_norm[lc.loocv_index], lc.solution_norm[lc.loocv_index],
               color="#C0392B", size=10, marker="circle", fill_alpha=0, line_width=2,
               legend_label="LOOCV")
        self._sel = points(p_l, lc.residual_norm[idx], lc.solution_norm[idx],
                           color=COLORS["ink"], size=8, legend_label="Selected")
        p_l.legend.location = "top_right"
        p_l.legend.label_text_font_size = "11px"

        p_res = plot(TIME_LABEL, "Residual", height=300)
        zero_line(p_res)
        self._res = line(p_res, spec.t, spec.lcurve_waveform - spec.tikhonovfits[idx], COLORS["muted"])

        p_fit = plot(TIME_LABEL, "Form factor", height=300)
        zero_line(p_fit)
        line(p_fit, spec.t, spec.lcurve_waveform, COLORS["muted"])
        self._fit = line(p_fit, spec.t, spec.tikhonovfits[idx], COLORS["fit"], width=2)

        p_dist = plot(DISTANCE_LABEL, "P(r)", height=300)
        line(p_dist, spec.r, spec.solutions[lc.loocv_index], "#C0392B", width=1,
             line_dash="dashed", legend_label="LOOCV")
        self._dist = line(p_dist, spec.r, spec.solutions[idx], COLORS["data"], width=2,
                          legend_label="Selected")
        p_dist.legend.label_text_font_size = "11px"

        main = [grid(card("L-curve", p_l), card("Fit residual", p_res),
                     card("Waveform and Tikhonov fit", p_fit),
                     card("Distance distribution", p_dist))]

        self._alpha = sidebar_widget(pn.widgets.IntSlider(
            name="Alpha index", start=0, end=len(spec.alphas) - 1, value=idx))
        self._alpha.param.watch(self._on_alpha_drag, "value")
        self._alpha.param.watch(self._on_alpha_commit, "value_throttled")
        self._alpha_text = pn.pane.HTML(self._alpha_html(idx), margin=0)
        loocv = button("Reset to LOOCV")
        loocv.on_click(lambda e: self._set_alpha(lc.loocv_index))
        side = [section("Regularization"), self._alpha, self._alpha_text, loocv,
                hint("The red circle is the leave-one-out cross-validation choice."),
                section("L-curve"), run]
        return main, side

    def _alpha_html(self, idx):
        lc = A.lcurve_of(self.spec)
        tag = " (LOOCV)" if idx == lc.loocv_index else ""
        return f'<div class="vn-kv"><b>Alpha</b>{self.spec.alphas[idx]:.4g}{tag}</div>'

    def _show(self, idx):
        spec = self.spec
        lc = A.lcurve_of(spec)
        self._sel.data = {"x": [lc.residual_norm[idx]], "y": [lc.solution_norm[idx]]}
        self._res.data = {"x": spec.t, "y": spec.lcurve_waveform - spec.tikhonovfits[idx]}
        self._fit.data = {"x": spec.t, "y": spec.tikhonovfits[idx]}
        self._dist.data = {"x": spec.r, "y": spec.solutions[idx]}
        self._alpha_text.object = self._alpha_html(idx)

    def _on_alpha_drag(self, event):
        self._show(int(event.new))

    def _on_alpha_commit(self, event):
        self._set_alpha(int(event.new))

    def _set_alpha(self, idx):
        if self._alpha.value != idx:
            self._alpha.value = idx
        A.select_alpha(self.spec, idx)
        self._show(idx)
        self.complete(alpha=float(self.spec.alpha))

    async def _on_compute(self, _event):
        spec = self.spec
        def then(lc):
            A.select_alpha(spec, lc.loocv_index)
            self.complete(alpha=float(spec.alpha), loocv_index=lc.loocv_index)
            self.refresh()

        await self.app.run_task("Computing L-curve", lambda: A.compute_lcurve(spec), then)
