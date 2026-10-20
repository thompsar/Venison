"""Results tab: best Monte Carlo model, confidence tables, and export."""

from __future__ import annotations

import html
import os

import numpy as np
import panel as pn
from bokeh.models import ColumnDataSource

from venison.core import analysis as A
from venison.core.export import CONFIDENCE_THRESHOLDS, export_results, gaussian_table
from venison.ui.common import (COLORS, DISTANCE_LABEL, TIME_LABEL, Tab, button, card, grid, hint,
                               line, plot, section, zero_line)
from venison.ui.project_tab import _reveal

_POP_COLORS = ["#0072B2", "#009E73", "#CC79A7", "#56B4E9", "#D55E00"]


class ResultsTab(Tab):
    title = "Results"
    icon = "results"
    requires = "montecarlo"

    def build(self):
        spec = self.spec
        env_lo, env_hi = A.best_monte_carlo_model(spec)

        p_fit = plot(TIME_LABEL, "Form factor", height=320)
        zero_line(p_fit)
        p_fit.varea("x", "y1", "y2", source=ColumnDataSource(
            {"x": spec.t, "y1": env_lo, "y2": env_hi}),
            fill_color=COLORS["fit"], fill_alpha=0.35, level="underlay",
            legend_label="95% band")
        line(p_fit, spec.t, spec.stable_waveform, COLORS["muted"], legend_label="Stable waveform")
        line(p_fit, spec.t, spec.best_model_fit, COLORS["data"], width=2, legend_label="Best model")
        p_fit.legend.label_text_font_size = "11px"

        p_pop = plot(DISTANCE_LABEL, "P(r)", height=320)
        line(p_pop, spec.r, spec.stable_solution / np.trapezoid(spec.stable_solution, spec.r),
             COLORS["muted"], legend_label="Tikhonov (stable)")
        for i, curve in enumerate(A.population_curves(spec)):
            line(p_pop, spec.r, curve, _POP_COLORS[i % len(_POP_COLORS)],
                 legend_label=f"Population {i + 1}")
        line(p_pop, spec.r, spec.best_model_solution, COLORS["fit"], width=2.5,
             legend_label="Sum")
        p_pop.legend.label_text_font_size = "11px"

        tables = [self._table_html(t) for t in sorted(CONFIDENCE_THRESHOLDS, reverse=True)]
        main = [grid(card("Best Monte Carlo fit", p_fit), card("Distance populations", p_pop)),
                card("Parameters (best value, lower and upper bound per F-test cutoff)",
                     pn.pane.HTML("".join(tables), sizing_mode="stretch_width",
                                  margin=(4, 14, 10, 14)))]

        export = button("Export results", "primary")
        export.on_click(self._on_export)
        reveal = button("Show results folder")
        reveal.on_click(lambda e: _reveal(self.app.project.path / "results"))
        self._written = pn.pane.HTML("", margin=0, sizing_mode="stretch_width")
        status = self.app.project.spectra()
        exported = any(s["name"] == self.session.name and s["exported"] for s in status)
        side = [section("Export"), export, reveal,
                hint("Writes the results and model tables (CSV) to <code>results/</code> and a "
                     "reloadable spectrum (.npz) to <code>processed/</code>."
                     + (" <b>Already exported</b>; exporting again overwrites." if exported
                        else "")),
                self._written]
        return main, side

    def _table_html(self, threshold):
        df = gaussian_table(self.spec, threshold)
        head = "".join(f"<th style='text-align:right;padding:2px 10px'>{html.escape(c)}</th>"
                       for c in ("r (nm)", "σ (nm)", "Mole fraction"))
        rows = []
        for pop, row in df.iterrows():
            cells = "".join(
                f"<td style='text-align:right;padding:2px 10px;font-variant-numeric:tabular-nums'>"
                f"{v[0]:.3f} <span style='color:#7A706A'>[{v[1]:.3f}, {v[2]:.3f}]</span></td>"
                for v in row)
            rows.append(f"<tr><td style='padding:2px 10px'>{pop}</td>{cells}</tr>")
        return (f"<div style='font-size:12.5px;font-weight:650;margin:10px 0 4px'>"
                f"F-test cutoff {threshold:.2f}</div>"
                f"<table style='font-size:12.5px;border-collapse:collapse'>"
                f"<tr style='color:#7A706A'><th style='text-align:left;padding:2px 10px'>Pop.</th>"
                f"{head}</tr>{''.join(rows)}</table>")

    def _on_export(self, _event):
        project = self.app.project
        try:
            written = export_results(self.spec, str(project.path / "results"),
                                     str(project.path / "processed"), name=self.session.name)
        except Exception as exc:
            self.app.notify(f"Export failed: {exc}", "error")
            return
        project.log("export", spectrum=self.session.name,
                    files=[os.path.relpath(p, project.path) for p in written])
        self._written.object = "<ul style='font-size:12px;padding-left:22px;margin:4px 0'>" + "".join(
            f"<li>{html.escape(os.path.relpath(p, project.path))}</li>" for p in written) + "</ul>"
        self.app.mark_dirty(after=self)
        self.app.notify(f"Exported {self.session.name}", "success")
