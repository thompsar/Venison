"""Main window: top bar, icon tab rail, per-tab sidebar, and content area."""

from __future__ import annotations

import asyncio
import html
import json
import os
import threading
import traceback
from functools import partial
from importlib.resources import files as res_files
from pathlib import Path

import panel as pn
import param

from venison import __version__
from venison.core.project import Project, default_projects_root
from venison.core.session import SpectrumSession
from venison.ui.icons import icon_svg

CONFIG_PATH = Path(os.environ.get("VENISON_CONFIG") or Path.home() / ".venison_config.json")


def _load_config() -> dict:
    try:
        return json.loads(CONFIG_PATH.read_text())
    except (FileNotFoundError, json.JSONDecodeError):
        return {}


def _save_config(config: dict) -> None:
    try:
        CONFIG_PATH.write_text(json.dumps(config, indent=2))
    except OSError as exc:
        print(f"Could not save settings to {CONFIG_PATH}: {exc}")


def _rail_template(tabs) -> str:
    parts = ['<div id="rail" class="vn-rail">']
    for i, tab in enumerate(tabs):
        if not tab.has_sidebar and i == len(tabs) - 1:
            parts.append('<div style="flex:1"></div>')
        elif i == 1:
            parts.append('<div class="vn-rail-sep"></div>')
        parts.append(f'<button class="vn-rail-btn" type="button" data-idx="{i}" '
                     f'title="{html.escape(tab.title)}">{icon_svg(tab.icon)}</button>')
    parts.append("</div>")
    return "".join(parts)


def make_rail_widget(tabs):
    """Build the rail as a ReactiveHTML class so its buttons are static HTML
    from the first paint. `done` is a comma list of tab indices showing a
    completion dot; `click_seq` increments on every click so re-clicking the
    active tab (to toggle the sidebar) still notifies Python."""

    sync = """
        rail.querySelectorAll("[data-idx]").forEach(function(btn) {
            var idx = btn.getAttribute("data-idx");
            btn.classList.toggle("vn-active", parseInt(idx) === data.active);
            btn.classList.toggle("vn-done", (data.done || "").split(",").indexOf(idx) !== -1);
        });
    """

    class RailWidget(pn.reactive.ReactiveHTML):
        active = param.Integer(default=0)
        done = param.String(default="")
        last_click = param.Integer(default=-1)
        click_seq = param.Integer(default=0)

        _template = _rail_template(tabs)
        _scripts = {
            "render": """
                rail.addEventListener("click", function(e) {
                    var btn = e.target.closest("[data-idx]");
                    if (!btn) return;
                    data.last_click = parseInt(btn.getAttribute("data-idx"));
                    data.click_seq = data.click_seq + 1;
                });
            """ + sync,
            "active": sync,
            "done": sync,
        }

    return RailWidget


class VenisonApp(param.Parameterized):
    project = param.Parameter(default=None, doc="Open Project, or None")
    session = param.Parameter(default=None, doc="SpectrumSession being analyzed, or None")
    active_tab = param.Integer(default=0)
    busy = param.String(default="", doc="Label of the running background task")

    def __init__(self, **params):
        super().__init__(**params)
        from venison.ui.tabs import TAB_CLASSES

        self.config = _load_config()
        self.projects_root = self.config.get("projects_root") or default_projects_root()
        self.tabs = [cls(self) for cls in TAB_CLASSES]
        self._build_top_bar()
        self._build_body()
        last = self.config.get("last_project")
        if last and Path(last, "raw_data").is_dir():
            try:
                self.open_project(last)
            except Exception as exc:
                print(f"Could not reopen {last}: {exc}")
        self.switch_tab(0)
        self._check_for_update()

    # ── Layout ───────────────────────────────────────────────────────────────

    def _build_top_bar(self):
        wordmark = pn.pane.HTML('<div class="vn-wordmark">Venison</div>',
                                margin=(0, 6, 0, 0), align="center")
        version = pn.pane.HTML(f'<div class="vn-version">v{__version__}</div>', align="center",
                               margin=0)
        self._project_btn = pn.widgets.Button(
            name="No project ▾", button_type="light", align="center", margin=(0, 4),
            description="Switch project", stylesheets=[_CHROMELESS_BUTTON])
        self._project_btn.on_click(lambda e: self.show_project_browser())
        self._spectrum_select = pn.widgets.Select(
            options=[], width=230, align="center", margin=(0, 4), visible=False)
        self._spectrum_select.param.watch(self._on_spectrum_picked, "value")
        self._status = pn.pane.HTML("", align="center", margin=(0, 4))
        self.param.watch(self._render_status, "busy")
        center = pn.Row(self._project_btn, self._spectrum_select, margin=0,
                        styles={"position": "absolute", "left": "50%",
                                "transform": "translateX(-50%)", "top": "8px"})
        self.top_bar = pn.Row(wordmark, version, center, pn.HSpacer(), self._status,
                              css_classes=["vn-topbar"], sizing_mode="stretch_width",
                              height=52, margin=0)

    def _build_body(self):
        self.rail = make_rail_widget(self.tabs)(width=64, sizing_mode="stretch_height", margin=0)
        self.rail.param.watch(lambda e: self.switch_tab(self.rail.last_click, toggle=True),
                              "click_seq")
        self.sidebar = pn.Column(width=290, sizing_mode="stretch_height", margin=0,
                                 css_classes=["vn-sidebar"])
        self.content = pn.Column(sizing_mode="stretch_both", margin=0,
                                 styles={"overflow-y": "auto", "min-height": "0"})
        self._update_banner = pn.Row(visible=False, sizing_mode="stretch_width", margin=0,
                                     css_classes=["vn-update"])
        self.layout = pn.Column(
            self.top_bar,
            self._update_banner,
            pn.Row(self.rail, self.sidebar, self.content, sizing_mode="stretch_both", margin=0),
            sizing_mode="stretch_both", margin=0, styles={"height": "100vh"},
        )

    # ── Update notice ────────────────────────────────────────────────────────

    def _check_for_update(self):
        """Ask GitHub for the latest release off the event loop (a few
        seconds at most, once per process) and show a banner if it's newer."""
        from venison.core.updates import newer_release

        doc = pn.state.curdoc

        def work():
            latest = newer_release(__version__)
            if latest:
                show = partial(self._show_update, latest)
                doc.add_next_tick_callback(show) if doc is not None else show()

        threading.Thread(target=work, daemon=True).start()

    def _show_update(self, latest: str):
        from venison.core.updates import INSTALL_LINE

        copy = pn.widgets.Button(name="Copy", button_type="light", margin=(0, 4),
                                 align="center", width=70)
        copy.js_on_click(args={"line": INSTALL_LINE}, code="""
            navigator.clipboard.writeText(line);
            cb_obj.label = "Copied";
        """)
        dismiss = pn.widgets.Button(name="Dismiss", button_type="light", margin=(0, 8, 0, 4),
                                    align="center", width=80)
        dismiss.on_click(lambda e: setattr(self._update_banner, "visible", False))
        message = pn.pane.HTML(
            f'<div class="vn-update-text"><b>Venison {html.escape(latest)} is available</b> '
            f'(you have {__version__}). To update, quit Venison and paste into Terminal: '
            f'<code>{html.escape(INSTALL_LINE)}</code></div>',
            sizing_mode="stretch_width", margin=(0, 0, 0, 14), align="center")
        self._update_banner[:] = [message, copy, dismiss]
        self._update_banner.visible = True

    def _render_status(self, _event=None, message: str = ""):
        text = self.busy or message
        cls = "vn-status vn-busy" if self.busy else "vn-status"
        self._status.object = f'<span class="{cls}">{html.escape(text)}</span>' if text else ""

    # ── Tabs ─────────────────────────────────────────────────────────────────

    def switch_tab(self, idx: int, toggle: bool = False):
        if not 0 <= idx < len(self.tabs):
            return
        tab = self.tabs[idx]
        if toggle and idx == self.active_tab and tab.has_sidebar:
            self.sidebar.visible = not self.sidebar.visible
            return
        self.active_tab = idx
        self.rail.active = idx
        if tab.dirty:
            self._refresh(tab)
        self.content[:] = [tab.main]
        self.sidebar[:] = [tab.sidebar]
        self.sidebar.visible = tab.has_sidebar

    def tab_for_step(self, step: str):
        return next(t for t in self.tabs if t.step == step)

    def _refresh(self, tab):
        try:
            tab.refresh()
        except Exception as exc:
            traceback.print_exc()
            tab.dirty = True
            self.notify(f"{tab.title}: {exc}", "error")

    def mark_dirty(self, after=None):
        """Mark tabs for rebuild on next visit (all, or those after `after`)."""
        start = self.tabs.index(after) + 1 if after is not None else 0
        for tab in self.tabs[start:]:
            tab.dirty = True
        self.tabs[0].dirty = True  # the project dashboard shows progress
        self._sync_rail()

    def on_step_completed(self, tab):
        self.mark_dirty(after=tab)

    def _sync_rail(self):
        done = []
        if self.session is not None:
            done = [str(i) for i, t in enumerate(self.tabs)
                    if t.step and self.session.is_done(t.step)]
        self.rail.done = ",".join(done)

    # ── Projects & spectra ───────────────────────────────────────────────────

    def set_projects_root(self, root: str):
        self.projects_root = root
        self.config["projects_root"] = root
        _save_config(self.config)

    def show_project_browser(self):
        self.tabs[0].show_browser()
        self.switch_tab(0)

    def _refuse_while_busy(self) -> bool:
        """Switching project or spectrum mid-computation would hand the result
        to the wrong spectrum; refuse and put the picker back."""
        if not self.busy:
            return False
        self.notify("Wait for the current computation to finish before switching.", "warning")
        self._syncing_select = True
        try:
            self._spectrum_select.value = self.session.name if self.session else None
        finally:
            self._syncing_select = False
        return True

    def open_project(self, path):
        if self._refuse_while_busy():
            return
        project = Project(path)
        self.project = project
        self.config["last_project"] = str(project.path)
        _save_config(self.config)
        self._project_btn.name = f"{project.name} ▾"
        self.tabs[0].browsing = False
        self.reload_spectra(select=project.last_spectrum)
        self.mark_dirty()
        if self.active_tab == 0:
            self.switch_tab(0)

    def reload_spectra(self, select: str | None = None):
        names = self.project.spectrum_names() if self.project else []
        current = self.session.name if self.session else None
        if select not in names:
            select = None
        target = select or (current if current in names else (names[0] if names else None))
        self._syncing_select = True
        try:
            self._spectrum_select.options = names
            self._spectrum_select.value = target
        finally:
            self._syncing_select = False
        self._spectrum_select.visible = bool(names)
        if target != current or self.project is not getattr(self.session, "project", None):
            self.open_spectrum(target)

    def _on_spectrum_picked(self, event):
        if getattr(self, "_syncing_select", False):
            return
        if event.new and (self.session is None or event.new != self.session.name):
            self.open_spectrum(event.new)

    def open_spectrum(self, name: str | None):
        if self._refuse_while_busy():
            return
        if name is None:
            self.session = None
        else:
            try:
                self.session = SpectrumSession(self.project, name)
            except Exception as exc:
                traceback.print_exc()
                self.notify(f"Could not open {name}: {exc}", "error")
                self.session = None
                name = None  # clear the picker, so picking it again retries
            else:
                self.project.last_spectrum = name
        self._syncing_select = True
        try:
            self._spectrum_select.value = name
        finally:
            self._syncing_select = False
        self.mark_dirty()
        self.switch_tab(self.active_tab)

    # ── Background work & messages ───────────────────────────────────────────

    async def run_task(self, label: str, fn, then=None):
        """Run `fn` off the event loop behind a busy overlay, then `then(result)`.

        `then` runs as a next-tick callback because code after an `await`
        doesn't hold the Bokeh document lock, and editing plot data without
        it raises. Errors are reported and `then` is skipped.
        """
        if self.busy:
            self.notify("Another computation is still running.", "warning")
            return
        doc = pn.state.curdoc
        session = self.session
        self.busy = f"{label}…"
        self.content.loading = True
        self.sidebar.styles = {"pointer-events": "none", "opacity": "0.6"}
        self._spectrum_select.disabled = self._project_btn.disabled = True
        print(f"{label}…")
        try:
            result = await asyncio.get_running_loop().run_in_executor(None, fn)
        except Exception as exc:
            traceback.print_exc()
            self.notify(f"{label} failed: {exc}", "error")
            return
        finally:
            self.busy = ""
            self.content.loading = False
            # Explicit values: an emptied dict doesn't clear styles in the browser.
            self.sidebar.styles = {"pointer-events": "auto", "opacity": "1"}
            self._spectrum_select.disabled = self._project_btn.disabled = False
        if then is not None:
            def finish():
                # `busy` is already clear by now, so a switch could slip in
                # before this tick; never apply a result to another spectrum.
                if self.session is not session:
                    self.notify(f"{label}: discarded, the open spectrum changed.", "warning")
                    return
                try:
                    then(result)
                except Exception as exc:
                    traceback.print_exc()
                    self.notify(f"{label}: {exc}", "error")
            if doc is not None:
                doc.add_next_tick_callback(finish)
            else:
                finish()

    def notify(self, message: str, kind: str = "info"):
        print(f"[{kind}] {message}")
        notifications = getattr(pn.state, "notifications", None)
        if notifications is not None:
            getattr(notifications, kind if kind in ("info", "success", "warning", "error")
                    else "info")(message, duration=6000 if kind == "error" else 3000)


_CHROMELESS_BUTTON = """
:host .bk-btn {
  background: transparent !important; border: none !important; box-shadow: none !important;
  font-size: 13.5px !important; font-weight: 650 !important; color: #26211E !important;
  padding: 6px 10px !important; border-radius: 6px !important;
}
:host .bk-btn:hover { background: rgba(38, 33, 30, 0.05) !important; }
"""


def configure_panel() -> None:
    """Register extensions and CSS. Must run before the server starts, so the
    page template ships their JS/CSS."""
    css = (res_files("venison.assets") / "layout.css").read_text(encoding="utf-8")
    pn.extension("tabulator", "filedropper", notifications=True, raw_css=[css])


def create_app() -> pn.viewable.Viewable:
    """Build one window's app (called per browser session)."""
    if pn.state.notifications is not None:
        pn.state.notifications.position = "bottom-right"
    return VenisonApp().layout
