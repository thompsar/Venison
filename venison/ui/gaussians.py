"""Gaussians tab: 1–5 Gaussian models of the stable distance distribution."""

from __future__ import annotations

import html

import numpy as np
import panel as pn
from bokeh.models import Span

from venison.core import analysis as A
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid, hint,
                               line, plot, points, section, sidebar_widget, zero_line)


class GaussiansTab(Tab):
    title = "Gaussians"
    icon = "gaussians"
    step = "gaussians"
    requires = "stability"

    def build(self):
        spec = self.spec
        self.gm = A.GaussianModels(spec, reset=getattr(spec, "gauss_fitted", None) is None)
        gm = self.gm
        self.model = int(getattr(self, "model", gm.best if np.isfinite(gm.bics).any() else 0))
        self.model = min(self.model, gm.ngauss - 1)

        p_dist = plot(DISTANCE_LABEL, "P(r)", height=300)
        line(p_dist, spec.r, spec.stable_solution, COLORS["muted"])
        self._seeds = points(p_dist, [], [], color="#C0392B", size=9)
        self._model_dist = line(p_dist, spec.r, spec.r * 0, COLORS["fit"], width=2)

        stable_rmsd = float(np.sqrt(np.mean((spec.stable_waveform - spec.stable_fit) ** 2)) * 1e2)
        p_rmsd = plot("Number of Gaussians", "RMSD (×100)", height=300,
                      x_range=(0.5, gm.ngauss + 0.5), y_range=(0, 1.5 * stable_rmsd))
        self._p_rmsd, self._stable_rmsd = p_rmsd, stable_rmsd
        p_rmsd.add_layout(Span(location=stable_rmsd, dimension="width", line_dash="dashed",
                               line_color=COLORS["ink"], line_width=1))
        self._rmsd_line = line(p_rmsd, [], [], COLORS["ink"], width=1)
        self._rmsd_pts = points(p_rmsd, [], [], color=COLORS["data"], size=8)
        self._best = points(p_rmsd, [], [], color="#C0392B", size=16, fill_alpha=0, line_width=2)

        p_wave = plot(TIME_LABEL, "Form factor", height=300)
        zero_line(p_wave)
        line(p_wave, spec.t, spec.stable_waveform, COLORS["muted"], line_alpha=0.8)
        line(p_wave, spec.t, spec.stable_fit, COLORS["data"], legend_label="Tikhonov")
        self._model_wave = line(p_wave, spec.t, spec.t * np.nan, COLORS["fit"], width=2,
                                legend_label="Gaussian model")
        p_wave.legend.label_text_font_size = "11px"

        p_res = plot(TIME_LABEL, "Residual", height=300)
        zero_line(p_res)
        self._res = line(p_res, spec.t, spec.stable_waveform - spec.stable_fit, COLORS["muted"])

        main = [grid(card("Stable distribution with seeds (red) and model", p_dist),
                     card("Fit RMSD by model (dashed: Tikhonov; circled: lowest BIC)", p_rmsd),
                     card("Waveform fits", p_wave),
                     card("Residual", p_res))]

        self._picker = sidebar_widget(pn.widgets.RadioButtonGroup(
            options={f"{m + 1}G": m for m in range(gm.ngauss)}, value=self.model,
            button_type="default"))
        self._picker.param.watch(self._on_model, "value")
        self._centers = pn.Column(margin=0, sizing_mode="stretch_width")
        init = button("Initialize Gaussian(s)", "primary")
        init.on_click(self._on_init)
        fit = button("Fit waveform")
        fit.on_click(self._on_fit)
        self._table = pn.pane.HTML("", margin=(4, 8), sizing_mode="stretch_width")
        side = [section("Model"), self._picker, section("Seed centers"), self._centers,
                init, fit,
                hint("<b>Initialize</b> fits the selected model to the distance distribution "
                     "from the seed centers. <b>Fit waveform</b> then refines every "
                     "initialized model against the waveform."),
                section("Models"), self._table]
        self._show_model()
        self._show_scores()
        return main, side

    # ── Rendering ────────────────────────────────────────────────────────────

    def _center_sliders(self):
        spec, gm, m = self.spec, self.gm, self.model
        dr = float(spec.r[1] - spec.r[0])
        sliders = []
        for i, r in enumerate(gm.model(m)[:, 0]):
            s = sidebar_widget(pn.widgets.FloatSlider(
                name=f"Gaussian {i + 1} (nm)", start=float(spec.r[0]), end=float(spec.r[-1]),
                step=dr, value=float(np.clip(r, spec.r[0], spec.r[-1])), format="0.00"))
            s.param.watch(lambda e, i=i: self._on_center(i, e.new), "value_throttled")
            sliders.append(s)
        self._centers[:] = sliders

    def _show_model(self):
        spec, gm, m = self.spec, self.gm, self.model
        centers = gm.model(m)[:, 0]
        nearest = [A.nearest_index(spec.r, c) for c in centers]
        self._seeds.data = {"x": centers, "y": spec.stable_solution[nearest]}
        if gm.fitted[m]:
            dist = gm.distribution(m)
            wave = gm.waveform(m)
            self._model_dist.data = {"x": spec.r, "y": dist}
            self._model_wave.data = {"x": spec.t, "y": wave}
            self._res.data = {"x": spec.t, "y": spec.stable_waveform - wave}
        else:
            self._model_dist.data = {"x": spec.r, "y": spec.r * np.nan}
            self._model_wave.data = {"x": spec.t, "y": spec.t * np.nan}
            self._res.data = {"x": spec.t, "y": spec.stable_waveform - spec.stable_fit}
        self._center_sliders()

    def _show_scores(self):
        gm = self.gm
        n = np.where(gm.fitted)[0]
        self._rmsd_line.data = {"x": n + 1, "y": gm.rmsds[n]}
        self._rmsd_pts.data = {"x": n + 1, "y": gm.rmsds[n]}
        top = max([self._stable_rmsd, *gm.rmsds[n]])
        self._p_rmsd.y_range.end = float(top * (1.1 if n.size else 1.5))
        if n.size:
            self._best.data = {"x": [gm.best + 1], "y": [gm.rmsds[gm.best]]}
        rows = []
        for m in range(gm.ngauss):
            state = ("optimized" if gm.optimized[m] else "initialized" if gm.fitted[m]
                     else "—")
            bic = f"{gm.bics[m]:.1f}" if gm.fitted[m] else ""
            star = " ★" if gm.fitted[m] and m == gm.best else ""
            rows.append(f"<tr><td>{m + 1}G{star}</td><td>{html.escape(state)}</td>"
                        f"<td style='text-align:right'>{bic}</td></tr>")
        self._table.object = (
            "<table style='width:100%;font-size:12.5px;border-collapse:collapse'>"
            "<tr style='color:#7A706A'><td>Model</td><td>State</td>"
            "<td style='text-align:right'>BIC</td></tr>" + "".join(rows) + "</table>")

    def _sync_step(self):
        if self.gm.optimized.any():
            self.complete(models=[int(m) + 1 for m in np.where(self.gm.optimized)[0]],
                          best=self.gm.best + 1)
        else:
            self.invalidate()
            self.session.save()

    # ── Callbacks ────────────────────────────────────────────────────────────

    def _on_model(self, event):
        self.model = int(event.new)
        self._show_model()

    def _on_center(self, i, value):
        self.gm.set_center(self.model, i, value)
        self._show_model()

    def _on_init(self, _event):
        try:
            self.gm.fit_to_distribution(self.model)
        except Exception as exc:
            self.app.notify(f"Initialization failed: {exc}", "error")
            return
        self._show_model()
        self._show_scores()
        self._sync_step()

    async def _on_fit(self, _event):
        gm = self.gm
        if not (gm.fitted & ~gm.optimized).any():
            self.app.notify("Initialize at least one model that hasn't been fit yet.", "warning")
            return
        def then(_):
            self._show_model()
            self._show_scores()
            self._sync_step()

        await self.app.run_task("Fitting Gaussian models", gm.optimize_pending, then)
