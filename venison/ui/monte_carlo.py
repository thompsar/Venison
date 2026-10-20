"""Monte Carlo tab: error surface of the chosen Gaussian model."""

from __future__ import annotations

import numpy as np
import panel as pn

from venison.core import analysis as A
from venison.ui.common import (CONFIDENCE_COLORS, COLORS, Tab, button, card, hint, plot, section,
                               sidebar_widget)
from bokeh.models import ColumnDataSource, Span

_PARAMS = ("r (nm)", "σ (nm)", "Mole fraction")


class MonteCarloTab(Tab):
    title = "Monte Carlo"
    icon = "montecarlo"
    step = "montecarlo"
    requires = "gaussians"

    def build(self):
        spec = self.spec
        gm = A.GaussianModels(spec, reset=False)
        optimized = [int(m) for m in np.where(gm.optimized)[0]]
        default = int(getattr(spec, "mc_model_id", gm.best))
        if default not in optimized:
            default = gm.best if gm.best in optimized else optimized[0]

        self._model = sidebar_widget(pn.widgets.Select(
            name="Model", options={f"{m + 1} Gaussian{'s' if m else ''}": m for m in optimized},
            value=default))
        self._model.param.watch(lambda e: self._render_weights(), "value")
        self._weights = pn.Column(margin=0, sizing_mode="stretch_width")
        reset = button("Reset weights")
        reset.on_click(lambda e: self._render_weights(reset=True))
        self._iterations = sidebar_widget(pn.widgets.IntInput(
            name="Iterations", value=int(getattr(spec, "mc_iterations", 5000)), start=100,
            step=1000))
        run = button("Run Monte Carlo", "primary")
        run.on_click(self._on_run)
        self._render_weights()
        side = [section("Search"), self._model, section("Step weights (r, σ)"), self._weights,
                reset, self._iterations, run,
                hint("Random walk around the fitted model; each population's r and σ steps are "
                     "scaled by its weight. Larger weights explore further. Defaults favor "
                     "small populations.")]

        if not self.session.is_done(self.step) or not hasattr(spec, "model_landscape"):
            return [card("", hint("Choose a model and click <b>Run Monte Carlo</b> in the "
                                  "sidebar."))], side
        return self._plots(), side

    def _render_weights(self, reset: bool = False):
        spec = self.spec
        m = int(self._model.value)
        saved = getattr(spec, "mc_rs_weights", None)
        if reset or saved is None or len(saved) != m + 1 or getattr(spec, "mc_model_id", -1) != m:
            weights = A.default_rs_weights(spec, m)
        else:
            weights = saved
        self._weights[:] = [sidebar_widget(pn.widgets.IntInput(
            name=f"Population {i + 1}", value=int(w), start=1)) for i, w in enumerate(weights)]

    def _plots(self):
        spec = self.spec
        scaled = A.mole_fraction_landscape(spec)
        stats = spec.landscape_statistics
        npop = scaled.shape[1]

        # Each sample is drawn once, colored by the tightest F-statistic level
        # it meets (the notebook overplotted every level, so samples were
        # drawn up to three times), sorted so tighter samples land on top,
        # and rendered with WebGL: tens of thousands of points per plot
        # otherwise stall the browser for seconds.
        levels = sorted(CONFIDENCE_COLORS)                    # 0.6, 0.75, 0.95
        shown = np.where(stats <= levels[-1])[0]
        shown = shown[np.argsort(-stats[shown])]
        tightest = np.asarray(levels)[np.searchsorted(levels, stats[shown])]
        colors = [CONFIDENCE_COLORS[l] for l in tightest]
        labels = [f"F statistic ≤ {l}" for l in tightest]

        def scatter(p, x, y, legend=False):
            src = ColumnDataSource({"x": x, "y": y, "color": colors, "label": labels})
            p.scatter("x", "y", source=src, color="color", size=2,
                      **({"legend_group": "label"} if legend else {}))

        p_surface = plot("Gaussian center (nm)", "Gaussian width (nm)", height=460,
                         output_backend="webgl")
        for pop in range(npop):
            scatter(p_surface, scaled[shown, pop, 0], scaled[shown, pop, 1], legend=pop == 0)
        p_surface.legend.label_text_font_size = "11px"

        nparams = 3 if npop > 1 else 2  # a single population's mole fraction is always 1
        rows = []
        for pop in range(npop):
            cells = []
            for j in range(nparams):
                p = plot(_PARAMS[j], "F statistic" if j == 0 else "", height=200,
                         y_range=(0, 1), title=f"Population {pop + 1}", output_backend="webgl")
                p.title.text_font_size = "11px"
                scatter(p, scaled[shown, pop, j], stats[shown])
                p.add_layout(Span(location=0.5, dimension="width", line_dash="dashed",
                                  line_color=COLORS["ink"], line_width=1))
                cells.append(pn.pane.Bokeh(p, sizing_mode="stretch_width"))
            rows.append(pn.Row(*cells, sizing_mode="stretch_width"))
        n = len(stats)
        return [card(f"Error surface: {npop} Gaussian{'s' if npop > 1 else ''}, {n:,} samples "
                     f"({len(shown):,} with F statistic ≤ {levels[-1]} shown)", p_surface),
                card("Parameter projections", *rows)]

    async def _on_run(self, _event):
        spec = self.spec
        m = int(self._model.value)
        weights = np.array([w.value for w in self._weights], dtype=float)
        iterations = int(self._iterations.value)
        def then(_):
            spec.mc_rs_weights = weights
            spec.mc_iterations = iterations
            self.complete(model=m + 1, iterations=iterations, weights=weights.tolist())
            self.refresh()

        await self.app.run_task(f"Running Monte Carlo ({iterations:,} iterations)",
                                lambda: A.run_monte_carlo(spec, m, weights, iterations), then)
