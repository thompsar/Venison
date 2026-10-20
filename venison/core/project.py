"""Projects on disk.

A project is a folder holding any number of spectra, analyzed one at a time:

    <root>/<project>/
        raw_data/<name>.DSC, <name>.DTA     ingested Bruker files
        raw_data/spectra.json               per-spectrum labels (sample, condition, notes)
        workspace/<name>.npz, <name>.json   in-progress analysis state + step status
        processed/<name>.npz                exported, reloadable spectrum (global analysis input)
        results/<name>_*.csv                exported tables
        analysis_log.json                   append-only record of analysis actions

Listing never writes; only explicit actions (create, ingest, edit) do.
"""

from __future__ import annotations

import json
import os
import re
import tempfile
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

SUBDIRS = ("raw_data", "workspace", "processed", "results")
MANIFEST = "spectra.json"
LOG = "analysis_log.json"
LABEL_FIELDS = ("sample", "condition", "notes")

_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _.+-]{0,99}$")


def default_projects_root() -> str:
    return os.environ.get("VENISON_PROJECTS_ROOT") or str(Path.home() / "venison_projects")


def validate_name(name: str) -> str:
    name = (name or "").strip()
    if not _NAME_RE.match(name) or name.endswith("."):
        raise ValueError(
            "Names may use letters, digits, spaces and _ . + - and must start "
            "with a letter or digit.")
    return name


def _write_json(path: Path, payload) -> None:
    """Atomic write so a crash never leaves a half-written file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp = tempfile.mkstemp(dir=path.parent, suffix=".tmp")
    with os.fdopen(fd, "w") as f:
        json.dump(payload, f, indent=2, default=str)
    os.replace(tmp, path)


def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def is_project_dir(path) -> bool:
    return (Path(path) / "raw_data").is_dir()


@dataclass(frozen=True)
class ProjectSummary:
    name: str
    path: str
    n_spectra: int
    n_exported: int
    last_modified: datetime | None


def list_projects(root) -> list[ProjectSummary]:
    root = Path(root)
    if not root.is_dir():
        return []
    out = []
    for p in sorted(root.iterdir(), key=lambda p: p.name.lower()):
        if not is_project_dir(p):
            continue
        proj = Project(p)
        names = proj.spectrum_names()
        mtimes = [f.stat().st_mtime for sub in SUBDIRS if (p / sub).is_dir()
                  for f in (p / sub).iterdir()]
        out.append(ProjectSummary(
            name=p.name, path=str(p.resolve()), n_spectra=len(names),
            n_exported=sum((p / "processed" / f"{n}.npz").exists() for n in names),
            last_modified=datetime.fromtimestamp(max(mtimes)) if mtimes else None,
        ))
    return out


def create_project(root, name) -> "Project":
    name = validate_name(name)
    path = Path(root) / name
    if path.exists():
        raise FileExistsError(f"A project named {name!r} already exists")
    for sub in SUBDIRS:
        (path / sub).mkdir(parents=True)
    proj = Project(path)
    _write_json(proj.manifest_path, {"created": _now(), "spectra": []})
    proj.log("create_project")
    return proj


class Project:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.name = self.path.name

    # ── Manifest ─────────────────────────────────────────────────────────────

    @property
    def manifest_path(self) -> Path:
        return self.path / "raw_data" / MANIFEST

    def _manifest(self) -> dict:
        try:
            with open(self.manifest_path) as f:
                manifest = json.load(f)
        except FileNotFoundError:
            manifest = {"spectra": []}
        # Raw files dropped in by hand still count, with empty labels.
        known = {s["name"] for s in manifest["spectra"]}
        for dta in sorted((self.path / "raw_data").glob("*.DTA")):
            if dta.stem not in known and dta.with_suffix(".DSC").exists():
                manifest["spectra"].append({"name": dta.stem, "added": None,
                                            **{k: "" for k in LABEL_FIELDS}})
        manifest["spectra"] = [s for s in manifest["spectra"]
                               if self.raw_path(s["name"]).exists()]
        return manifest

    def spectra(self) -> list[dict]:
        """One dict per spectrum: name, labels, and analysis status."""
        rows = []
        for s in self._manifest()["spectra"]:
            status = self.read_status(s["name"])
            rows.append({**s, "steps": status.get("done", []),
                         "exported": (self.path / "processed" / f"{s['name']}.npz").exists(),
                         "updated": status.get("updated")})
        return rows

    def spectrum_names(self) -> list[str]:
        return [s["name"] for s in self._manifest()["spectra"]]

    def raw_path(self, name: str) -> Path:
        return self.path / "raw_data" / f"{name}.DTA"

    def add_spectrum(self, name: str, dsc: bytes, dta: bytes, *, overwrite=False, **labels) -> None:
        name = validate_name(name)
        if self.raw_path(name).exists() and not overwrite:
            raise FileExistsError(f"{name} is already in this project")
        raw = self.path / "raw_data"
        (raw / f"{name}.DSC").write_bytes(dsc)
        (raw / f"{name}.DTA").write_bytes(dta)
        manifest = self._manifest()
        manifest["spectra"] = [s for s in manifest["spectra"] if s["name"] != name]
        manifest["spectra"].append({"name": name, "added": _now(),
                                    **{k: str(labels.get(k, "") or "") for k in LABEL_FIELDS}})
        _write_json(self.manifest_path, manifest)
        if overwrite:
            self.clear_workspace(name)
        self.log("ingest", spectrum=name, **labels)

    def update_labels(self, name: str, **labels) -> None:
        manifest = self._manifest()
        for s in manifest["spectra"]:
            if s["name"] == name:
                s.update({k: str(v) for k, v in labels.items() if k in LABEL_FIELDS})
        _write_json(self.manifest_path, manifest)

    def remove_spectrum(self, name: str) -> None:
        """Remove a spectrum's raw files and working state. Exports are kept."""
        for p in (self.path / "raw_data" / f"{name}.DSC", self.raw_path(name)):
            p.unlink(missing_ok=True)
        self.clear_workspace(name)
        manifest = self._manifest()
        _write_json(self.manifest_path, manifest)
        self.log("remove_spectrum", spectrum=name)

    # ── Workspace (in-progress analysis) ─────────────────────────────────────

    def workspace_paths(self, name: str) -> tuple[Path, Path]:
        ws = self.path / "workspace"
        return ws / f"{name}.npz", ws / f"{name}.json"

    def read_status(self, name: str) -> dict:
        try:
            with open(self.workspace_paths(name)[1]) as f:
                return json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            return {}

    def write_status(self, name: str, status: dict) -> None:
        _write_json(self.workspace_paths(name)[1], {**status, "updated": _now()})

    def clear_workspace(self, name: str) -> None:
        for p in self.workspace_paths(name):
            p.unlink(missing_ok=True)

    # ── Log ──────────────────────────────────────────────────────────────────

    def log(self, step: str, **data) -> None:
        path = self.path / LOG
        try:
            with open(path) as f:
                entries = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError):
            entries = []
        entries.append({"timestamp": _now(), "step": step, **data})
        _write_json(path, entries)

    @property
    def last_spectrum(self) -> str | None:
        name = self._manifest().get("last_spectrum")
        return name if name in self.spectrum_names() else None

    @last_spectrum.setter
    def last_spectrum(self, name: str) -> None:
        manifest = self._manifest()
        manifest["last_spectrum"] = name
        _write_json(self.manifest_path, manifest)
