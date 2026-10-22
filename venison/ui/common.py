"""Shared UI building blocks: plot factory, cards, sidebar labels, tab base class."""

from __future__ import annotations

import html
from typing import TYPE_CHECKING

import numpy as np
import panel as pn
from bokeh.models import ColumnDataSource, Span
from bokeh.plotting import figure

if TYPE_CHECKING:
    from venison.app import VenisonApp

# Okabe–Ito: distinguishable for common color-vision deficiencies.
COLORS = {
    "data": "#0072B2",
    "imag": "#D55E00",
    "fit": "#E69F00",
    "alt": "#009E73",
    "accent": "#CC79A7",
    "sky": "#56B4E9",
    "muted": "#9A928C",
    "ink": "#26211E",
    "highlight": "#D55E00",   # the algorithm's pick (LOOCV, RMSD minimum, best BIC)
}
# Tight to loose reads dark to faded; distinct under common color blindness.
CONFIDENCE_COLORS = {0.6: "#0072B2", 0.75: "#D55E00", 0.95: "#C9C3BD"}

TIME_LABEL = "Evolution time (µs)"
DISTANCE_LABEL = "Distance (nm)"


def plot(x_label: str = "", y_label: str = "", height: int = 300, **kwargs):
    p = figure(height=height, sizing_mode="stretch_width",
               tools="pan,box_zoom,wheel_zoom,reset,save", toolbar_location="right", **kwargs)
    p.toolbar.logo = None
    p.toolbar.autohide = True
    p.outline_line_alpha = 0
    p.border_fill_alpha = 0
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_alpha = 0.5
    p.ygrid.grid_line_dash = [6, 4]
    for axis, label in ((p.xaxis, x_label), (p.yaxis, y_label)):
        axis.axis_label = label
        axis.axis_label_text_font_style = "normal"
        axis.axis_label_text_font_size = "12px"
        axis.major_label_text_font_size = "11px"
        axis.axis_line_color = "#BDB5AE"
        axis.major_tick_line_color = "#BDB5AE"
        axis.minor_tick_line_color = None
    return p


def line(p, x, y, color=COLORS["data"], width=1.5, **kwargs) -> ColumnDataSource:
    src = ColumnDataSource({"x": np.asarray(x), "y": np.asarray(y)})
    p.line("x", "y", source=src, line_color=color, line_width=width, **kwargs)
    return src


def points(p, x, y, color=COLORS["data"], size=6, **kwargs) -> ColumnDataSource:
    src = ColumnDataSource({"x": np.atleast_1d(x), "y": np.atleast_1d(y)})
    p.scatter("x", "y", source=src, color=color, size=size, **kwargs)
    return src


def vline(p, x, color=COLORS["ink"], dash="dashed", width=1.2) -> Span:
    span = Span(location=float(x), dimension="height", line_color=color,
                line_dash=dash, line_width=width)
    p.add_layout(span)
    return span


def zero_line(p) -> Span:
    span = Span(location=0, dimension="width", line_color="#8D847E",
                line_dash="dashed", line_width=1)
    p.add_layout(span)
    return span


def card(title: str, *objects, **kwargs) -> pn.Column:
    items = [pn.pane.HTML(f'<div class="vn-card-title">{html.escape(title)}</div>',
                          margin=0, sizing_mode="stretch_width")] if title else []
    for obj in objects:
        if hasattr(obj, "renderers"):  # a Bokeh figure
            obj = pn.pane.Bokeh(obj, sizing_mode="stretch_width", margin=(0, 6, 6, 6))
        items.append(obj)
    kwargs.setdefault("sizing_mode", "stretch_width")
    return pn.Column(*items, css_classes=["vn-card"], margin=(0, 0, 14, 0), **kwargs)


def grid(*cards, ncols: int = 2) -> pn.GridBox:
    return pn.GridBox(*cards, ncols=ncols, sizing_mode="stretch_width")


def section(text: str) -> pn.pane.HTML:
    return pn.pane.HTML(f'<span class="vn-section">{html.escape(text)}</span>',
                        sizing_mode="stretch_width", margin=0)


def hint(text: str) -> pn.pane.HTML:
    return pn.pane.HTML(f'<div class="vn-hint">{text}</div>',
                        sizing_mode="stretch_width", margin=0)


def page_header(title: str, subtitle: str = "") -> pn.pane.HTML:
    sub = f'<div class="vn-page-sub">{subtitle}</div>' if subtitle else ""
    return pn.pane.HTML(f'<h1 class="vn-page-title">{html.escape(title)}</h1>{sub}',
                        sizing_mode="stretch_width", margin=(4, 0, 14, 0))


def empty_state(title: str, body: str) -> pn.pane.HTML:
    return pn.pane.HTML(
        f'<div class="vn-empty"><h2>{html.escape(title)}</h2><p>{body}</p></div>',
        sizing_mode="stretch_width")


def sidebar_widget(widget):
    widget.sizing_mode = "stretch_width"
    widget.margin = (4, 8, 2, 8)
    return widget


def button(name: str, kind: str = "default", **kwargs) -> pn.widgets.Button:
    return sidebar_widget(pn.widgets.Button(name=name, button_type=kind, **kwargs))


class Tab:
    """One rail tab: a main area and a sidebar, rebuilt when the open spectrum
    changes. Subclasses implement `build()`; `requires` names the workflow
    step that must be complete first (None: just needs an open spectrum)."""

    title = ""
    icon = ""
    step: str | None = None       # the workflow step this tab completes
    requires: str | None = None
    needs_spectrum = True
    has_sidebar = True

    def __init__(self, app: "VenisonApp"):
        self.app = app
        self.main = pn.Column(sizing_mode="stretch_width", margin=(18, 22, 30, 22))
        self.sidebar = pn.Column(sizing_mode="stretch_width", margin=0)
        self.dirty = True

    @property
    def session(self):
        return self.app.session

    @property
    def spec(self):
        return self.app.session.spec

    def refresh(self) -> None:
        self.dirty = False
        title = pn.pane.HTML(f'<span class="vn-sidebar-title">{html.escape(self.title)}</span>',
                             sizing_mode="stretch_width", margin=0)
        if self.needs_spectrum and self.session is None:
            self.main[:] = [empty_state(
                "No spectrum open",
                "Open a project and pick a spectrum in the Project tab.")]
            self.sidebar[:] = [title]
            return
        if self.requires and not self.session.is_done(self.requires):
            from venison.core.session import STEP_LABELS
            tab = self.app.tab_for_step(self.requires)
            self.main[:] = [empty_state(
                f"{STEP_LABELS[self.requires]} comes first",
                f"Finish the <b>{html.escape(tab.title)}</b> tab for "
                f"<b>{html.escape(self.session.name)}</b>, then come back here.")]
            self.sidebar[:] = [title]
            return
        main, side = self.build()
        self.main[:] = list(main)
        self.sidebar[:] = [title] + list(side)

    def build(self) -> tuple[list, list]:
        raise NotImplementedError

    def complete(self, **log) -> None:
        self.session.complete(self.step, **log)
        self.app.on_step_completed(self)

    def invalidate(self) -> None:
        """This tab's inputs changed: its step and everything after need a rerun."""
        if self.step:
            self.session.invalidate_from(self.step)
            self.app.on_step_completed(self)


def stale_banner() -> pn.pane.HTML:
    return pn.pane.HTML(
        '<div class="vn-stale">An earlier step changed since this was computed. '
        "Rerun this step to bring it up to date.</div>",
        sizing_mode="stretch_width", margin=(0, 0, 12, 0))
