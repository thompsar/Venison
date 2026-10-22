"""Console tab: everything printed by the app and its analysis code."""

from __future__ import annotations

import html

import panel as pn

from venison.core import log_capture
from venison.ui.common import Tab, card, page_header


class ConsoleTab(Tab):
    title = "Console"
    icon = "console"
    needs_spectrum = False
    has_sidebar = False

    def __init__(self, app):
        super().__init__(app)
        self._pane = pn.pane.HTML("", sizing_mode="stretch_width", margin=0)
        self._seen = -1
        self._callback = None

    def build(self):
        self._update()
        if self._callback is None:
            self._callback = pn.state.add_periodic_callback(self._poll, period=1000)
        return [page_header("Console", "Output from the app and its computations."),
                card("", self._pane)], []

    def _poll(self):
        if self.app.tabs[self.app.active_tab] is self:
            self._update()

    def _update(self):
        if log_capture.version == self._seen:
            return
        self._seen = log_capture.version
        lines = log_capture.lines()[-500:] or ["(nothing yet)"]
        self._pane.object = (f'<pre class="vn-console">'
                             f'{html.escape(chr(10).join(lines))}</pre>')
