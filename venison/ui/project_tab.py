"""Project tab: the project browser, the open project's spectra, and ingest."""

from __future__ import annotations

import html
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import panel as pn

from venison.core.project import create_project, list_projects, validate_name
from venison.core.session import STEP_LABELS, STEPS
from venison.ui.common import Tab, button, card, hint, page_header, section, sidebar_widget

_PILL = ("display:inline-block;padding:1px 11px;border-radius:10px;font-weight:600;"
         "font-size:12px;cursor:pointer;")
_ACTION = f"<span style='{_PILL}background:rgba(154,75,44,.1);color:#9A4B2C'>{{}}</span>"
_CURRENT = f"<span style='{_PILL}background:#9A4B2C;color:#fff'>{{}}</span>"
_OK = "<span style='color:#3D8B5A;font-weight:600'>✓</span>"

_TABLE_CSS = [""".tabulator { border: none !important; font-size: 13px; }
.tabulator .tabulator-header { border-bottom: 1px solid #E6E1DC !important; background: #fff !important; }
.tabulator .tabulator-header .tabulator-col { background: #fff !important; }
.tabulator-row { border-bottom: 1px solid #F1EDE9 !important; }"""]


def _reveal(path: Path) -> None:
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif sys.platform == "win32":
        os.startfile(str(path))  # noqa: S606 - local desktop folder
    else:
        subprocess.Popen(["xdg-open", str(path)])


def progress_text(steps: list[str]) -> str:
    if not steps:
        return "Not started"
    last = max(steps, key=STEPS.index)
    return f"{len(steps)}/{len(STEPS)} · {STEP_LABELS[last]}"


class ProjectTab(Tab):
    title = "Project"
    icon = "project"
    needs_spectrum = False

    def __init__(self, app):
        super().__init__(app)
        self.browsing = False
        self._staged: dict[str, dict[str, bytes]] = {}
        self._confirm_remove = False

    def show_browser(self):
        self.browsing = True
        self.dirty = True

    def build(self):
        if self.app.project is None or self.browsing:
            return self._build_browser()
        return self._build_dashboard()

    # ── Browser ──────────────────────────────────────────────────────────────

    def _build_browser(self):
        app = self.app
        summaries = list_projects(app.projects_root)
        current = str(app.project.path) if app.project else ""
        self._paths = {s.name: s.path for s in summaries}

        new_name = pn.widgets.TextInput(placeholder="New project name", width=280)
        create = pn.widgets.Button(name="Create project", button_type="primary")

        def on_create(_e):
            try:
                proj = create_project(app.projects_root, new_name.value)
            except (ValueError, FileExistsError, OSError) as exc:
                app.notify(str(exc), "error")
                return
            app.open_project(proj.path)
            app.notify(f"Created {proj.name}. Add spectra below.", "success")

        create.on_click(on_create)

        def frame(filter_text=""):
            rows = []
            for s in summaries:
                if filter_text and filter_text.lower() not in s.name.lower():
                    continue
                is_open = s.path == current
                rows.append({
                    "": (_CURRENT if is_open else _ACTION).format("Return" if is_open else "Open"),
                    "Project": html.escape(s.name),
                    "Spectra": s.n_spectra,
                    "Exported": s.n_exported,
                    "Last worked on": s.last_modified.strftime("%Y-%m-%d %H:%M")
                    if s.last_modified else "",
                    "_name": s.name,
                })
            return pd.DataFrame(rows, columns=["", "Project", "Spectra", "Exported",
                                               "Last worked on", "_name"])

        table = pn.widgets.Tabulator(
            frame(), show_index=False, disabled=True, selectable=False, hidden_columns=["_name"],
            formatters={"": {"type": "html"}, "Project": {"type": "html"}},
            widths={"": 90, "Spectra": 90, "Exported": 90, "Last worked on": 160},
            sizing_mode="stretch_width", layout="fit_columns", stylesheets=_TABLE_CSS,
            margin=(4, 8, 8, 8))

        def on_click(event):
            if event.column != "":
                return
            name = table.value.iloc[event.row]["_name"]
            if self._paths.get(name) == current:
                self.browsing = False
                self.dirty = True
                app.switch_tab(0)
            else:
                app.open_project(self._paths[name])

        table.on_click(on_click)
        search = pn.widgets.TextInput(placeholder="Filter projects", width=240, margin=(10, 8, 0, 14))
        search.param.watch(lambda e: setattr(table, "value", frame(e.new)), "value_input")

        listing = card("", pn.Row(search, margin=0), table) if summaries else card(
            "", hint("No projects here yet. Create one above, or point the sidebar at the "
                     "folder that holds your projects."))
        main = [
            page_header("Projects", html.escape(app.projects_root)),
            pn.Row(new_name, create, margin=(0, 0, 14, 0)),
            listing,
        ]

        root_input = sidebar_widget(pn.widgets.TextInput(value=app.projects_root))
        use_root = button("Use this folder")

        def on_root(_e):
            path = os.path.expanduser(root_input.value.strip())
            if not os.path.isdir(path):
                try:
                    os.makedirs(path)
                except OSError as exc:
                    app.notify(f"Can't use {path}: {exc}", "error")
                    return
            app.set_projects_root(path)
            self.dirty = True
            app.switch_tab(0)

        use_root.on_click(on_root)
        side = [section("Projects folder"), root_input, use_root,
                hint("Each project is a folder in here holding its spectra and results.")]
        if app.project is not None:
            back = button(f"Back to {app.project.name}", "primary")

            def on_back(_e):
                self.browsing = False
                self.dirty = True
                app.switch_tab(0)

            back.on_click(on_back)
            side = [section("Open project"), back] + side
        return main, side

    # ── Dashboard ────────────────────────────────────────────────────────────

    def _build_dashboard(self):
        app = self.app
        project = app.project
        spectra = project.spectra()
        current = app.session.name if app.session else None

        rows = [{
            "": (_CURRENT if s["name"] == current else _ACTION).format(
                "Current" if s["name"] == current else "Analyze"),
            "Spectrum": s["name"],
            "Sample": s.get("sample", ""),
            "Condition": s.get("condition", ""),
            "Notes": s.get("notes", ""),
            "Progress": progress_text(s["steps"]),
            "Exported": _OK if s["exported"] else "",
        } for s in spectra]
        df = pd.DataFrame(rows, columns=["", "Spectrum", "Sample", "Condition", "Notes",
                                         "Progress", "Exported"])
        table = pn.widgets.Tabulator(
            df, show_index=False, selectable=False, sizing_mode="stretch_width",
            layout="fit_columns", stylesheets=_TABLE_CSS, margin=(4, 8, 8, 8),
            formatters={"": {"type": "html"}, "Exported": {"type": "html"}},
            editors={c: None for c in ("", "Spectrum", "Progress", "Exported")},
            widths={"": 96, "Progress": 170, "Exported": 84})

        def on_click(event):
            if event.column == "":
                name = table.value.iloc[event.row]["Spectrum"]
                if name != current:
                    app.open_spectrum(name)
                app.switch_tab(1)

        def on_edit(event):
            name = table.value.iloc[event.row]["Spectrum"]
            project.update_labels(name, **{event.column.lower(): event.value})

        table.on_click(on_click)
        table.on_edit(on_edit)

        spectra_card = card(
            "Spectra",
            hint("Click <b>Analyze</b> to work on a spectrum; the analysis tabs on the left "
                 "follow the spectrum picked here or in the top bar. Sample, Condition and "
                 "Notes are editable.") if spectra else hint(
                "No spectra yet. Add .DSC/.DTA pairs below."),
            *([table] if spectra else []))

        main = [page_header(project.name, html.escape(str(project.path))),
                spectra_card, self._build_ingest()]

        side = [section("Project"),
                pn.pane.HTML(f'<div class="vn-kv"><b>Spectra</b>{len(spectra)}</div>'
                             f'<div class="vn-kv"><b>Exported</b>'
                             f'{sum(s["exported"] for s in spectra)}</div>', margin=0)]
        switch = button("Switch project…")
        switch.on_click(lambda e: app.show_project_browser())
        reveal = button("Show folder")
        reveal.on_click(lambda e: _reveal(project.path))
        side += [switch, reveal]
        if app.session is not None:
            side += self._current_spectrum_section(spectra)
        return main, side

    def _current_spectrum_section(self, spectra):
        app = self.app
        name = app.session.name
        info = next((s for s in spectra if s["name"] == name), {})
        kv = "".join(f'<div class="vn-kv"><b>{k}</b>{html.escape(str(v) or "—")}</div>' for k, v in (
            ("Name", name), ("Sample", info.get("sample", "")),
            ("Condition", info.get("condition", "")),
            ("Progress", progress_text(info.get("steps", [])))))
        remove = button("Remove from project…", "danger" if self._confirm_remove else "default")
        if self._confirm_remove:
            remove.name = f"Confirm: delete {name}'s raw files"

        def on_remove(_e):
            if not self._confirm_remove:
                self._confirm_remove = True
                self.refresh()
                return
            self._confirm_remove = False
            app.project.remove_spectrum(name)
            app.reload_spectra()
            app.notify(f"Removed {name}", "info")
            self.refresh()

        remove.on_click(on_remove)
        out = [section("Current spectrum"), pn.pane.HTML(kv, margin=0), remove]
        if self._confirm_remove:
            out.append(hint("Deletes its raw files and saved analysis from this project. "
                            "Exported results are kept."))
        return out

    # ── Ingest ───────────────────────────────────────────────────────────────

    def _build_ingest(self):
        app = self.app
        self._dropper_holder = pn.Column(sizing_mode="stretch_width", margin=(6, 10, 4, 10))
        self._staging = pn.pane.HTML("", sizing_mode="stretch_width", margin=(0, 14))
        self._labels = {k: pn.widgets.TextInput(name=k.title(), placeholder="optional",
                                                sizing_mode="stretch_width")
                        for k in ("sample", "condition", "notes")}
        self._replace = pn.widgets.Checkbox(name="Replace spectra that already exist "
                                                 "(discards their saved analysis)")
        self._add = pn.widgets.Button(name="Add spectra", button_type="primary", disabled=True)
        self._clear = pn.widgets.Button(name="Clear", button_type="light")
        self._add.on_click(self._on_add)
        self._clear.on_click(lambda e: self._reset_staging())
        self._replace.param.watch(lambda e: self._render_staging(), "value")
        self._new_dropper()
        self._render_staging()
        return card(
            "Add spectra",
            hint("Drop Bruker .DSC and .DTA files (both files of each pair). Labels entered "
                 "here apply to every spectrum in this batch and can be edited later."),
            self._dropper_holder, self._staging,
            pn.Row(*self._labels.values(), sizing_mode="stretch_width", margin=(0, 8)),
            pn.Row(self._add, self._clear, self._replace, margin=(4, 8, 12, 8)),
        )

    def _new_dropper(self):
        dropper = pn.widgets.FileDropper(
            multiple=True, layout="compact", sizing_mode="stretch_width",
            stylesheets=[".filepond--file-action-button { display: none !important; }"],
            styles={"border": "2px dashed #D9D2CB", "border-radius": "8px", "background": "#FCFBFA"})
        dropper.param.watch(self._on_drop, "value")
        self._dropper_holder[:] = [dropper]

    def _on_drop(self, event):
        rejected = []
        for fname, data in (event.new or {}).items():
            stem, ext = os.path.splitext(os.path.basename(fname))
            ext = ext.upper()
            if ext not in (".DSC", ".DTA"):
                rejected.append(f"{fname}: not a .DSC or .DTA file")
                continue
            try:
                stem = validate_name(stem)
            except ValueError:
                rejected.append(f"{fname}: unsupported characters in the file name")
                continue
            self._staged.setdefault(stem, {})[ext] = data
        for msg in rejected:
            self.app.notify(msg, "warning")
        self._render_staging()

    def _ready(self) -> list[str]:
        existing = set(self.app.project.spectrum_names())
        return [n for n, f in self._staged.items()
                if len(f) == 2 and (n not in existing or self._replace.value)]

    def _render_staging(self):
        if not self._staged:
            self._staging.object = ""
            self._add.disabled = True
            return
        existing = set(self.app.project.spectrum_names())
        ready = self._ready()
        items = []
        for name, files in sorted(self._staged.items()):
            missing = {".DSC", ".DTA"} - set(files)
            if missing:
                state = f"waiting for {', '.join(sorted(missing))}"
            elif name in existing:
                state = "will replace existing" if self._replace.value else "already in project, skipped"
            else:
                state = "ready"
            color = "#3D8B5A" if name in ready else "#7A706A"
            items.append(f"<li><b>{html.escape(name)}</b> "
                         f"<span style='color:{color}'>{state}</span></li>")
        self._staging.object = (f"<ul style='font-size:13px;margin:6px 0;padding-left:18px'>"
                                f"{''.join(items)}</ul>")
        self._add.disabled = not ready
        self._add.name = f"Add {len(ready)} spectr{'um' if len(ready) == 1 else 'a'}"

    def _reset_staging(self):
        self._staged = {}
        self._new_dropper()
        self._render_staging()

    def _on_add(self, _event):
        app = self.app
        labels = {k: w.value.strip() for k, w in self._labels.items()}
        added = []
        for name in self._ready():
            files = self._staged[name]
            try:
                app.project.add_spectrum(name, files[".DSC"], files[".DTA"],
                                         overwrite=self._replace.value, **labels)
            except Exception as exc:
                app.notify(f"Could not add {name}: {exc}", "error")
                continue
            added.append(name)
        if added:
            app.notify(f"Added {', '.join(added)}", "success")
            # Re-open even if it was current: a replaced spectrum starts over.
            app.session = None
            app.reload_spectra(select=added[-1])
        self._staged = {}
        self.dirty = True
        app.switch_tab(0)
