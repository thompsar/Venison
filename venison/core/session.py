"""One spectrum's analysis in progress: the `spectrum` object plus which
workflow steps are complete, saved to the project's workspace after each step.
"""

from __future__ import annotations

import numpy as np

from venison.core.analysis import BG_MODELS
from venison.core.project import Project
from venison.core.spectrum import spectrum

# Workflow order. Changing a step invalidates every step after it.
STEPS = ("kernel", "background", "lcurve", "validation", "stability", "gaussians", "montecarlo")

STEP_LABELS = {
    "kernel": "Phase & cutoff",
    "background": "Background",
    "lcurve": "L-curve",
    "validation": "Validation",
    "stability": "Stability",
    "gaussians": "Gaussian models",
    "montecarlo": "Monte Carlo",
}

# Attributes that can't round-trip through savez; rebuilt on load instead.
_UNSAVED = ("bgmodel",)


class SpectrumSession:
    def __init__(self, project: Project, name: str):
        self.project = project
        self.name = name
        self.raw_path = str(project.raw_path(name))
        npz, _ = project.workspace_paths(name)
        status = project.read_status(name)
        self.done: set[str] = set()
        if npz.exists():
            try:
                self.spec = spectrum(str(npz))
                self.done = set(status.get("done", [])) & set(STEPS)
            except Exception as exc:  # corrupt or incompatible workspace: start over
                print(f"Could not restore saved analysis for {name} ({exc}); starting fresh.")
                self.spec = spectrum(self.raw_path)
        else:
            self.spec = spectrum(self.raw_path)
        self.spec.path = self.raw_path
        if hasattr(self.spec, "bgmodel_name"):
            self.spec.bgmodel = BG_MODELS[str(self.spec.bgmodel_name)]

    def is_done(self, step: str) -> bool:
        return step in self.done

    def complete(self, step: str, **log) -> None:
        """Mark `step` done, invalidate everything after it, save, and log."""
        self.done = {s for s in self.done if STEPS.index(s) < STEPS.index(step)} | {step}
        self.save()
        self.project.log(step, spectrum=self.name, **log)

    def invalidate_from(self, step: str) -> None:
        """Mark `step` and everything after it as needing a rerun."""
        before = len(self.done)
        self.done = {s for s in self.done if STEPS.index(s) < STEPS.index(step)}
        if len(self.done) != before:
            self.project.write_status(self.name, {"done": sorted(self.done, key=STEPS.index)})

    def save(self) -> None:
        npz, _ = self.project.workspace_paths(self.name)
        npz.parent.mkdir(parents=True, exist_ok=True)
        state = {k: v for k, v in vars(self.spec).items()
                 if k not in _UNSAVED and not k.startswith("_")}  # "_": caches
        tmp = npz.with_suffix(".tmp.npz")
        np.savez(tmp, **state)
        tmp.replace(npz)
        self.project.write_status(self.name, {"done": sorted(self.done, key=STEPS.index)})
