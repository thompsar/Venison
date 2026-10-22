"""Write analysis results: the per-spectrum CSVs and the reloadable .npz."""

from __future__ import annotations

import os

import numpy as np
import pandas as pd

from venison.core.analysis import mole_fraction_landscape, population_curves

RESULT_ARRAYS = [
    "full_t", "raw_time", "real", "imag", "background", "t", "waveform", "tikhonovfit",
    "stable_waveform", "unstable_waveform", "stable_fit", "best_model_fit",
    "r", "solution", "stable_solution", "best_model_solution",
]

CONFIDENCE_THRESHOLDS = (0.6, 0.75, 0.95)


def results_frame(spec) -> pd.DataFrame:
    """Every trace side by side; columns differ in length, so they're padded."""
    out = pd.DataFrame({name: pd.Series(getattr(spec, name)) for name in RESULT_ARRAYS})
    out["Rescaled Tik"] = pd.Series(
        spec.stable_solution / np.trapezoid(spec.stable_solution, spec.r))
    for p, curve in enumerate(population_curves(spec), start=1):
        out[f"Population: {p}"] = pd.Series(curve)
    return out


def gaussian_table(spec, threshold: float) -> pd.DataFrame:
    """Best-model parameters with min/max bounds over models under `threshold`."""
    scaled = mole_fraction_landscape(spec)
    best = scaled[np.argmin(spec.landscape_statistics)]
    within = scaled[spec.landscape_statistics < threshold]
    columns = ["[r (nm), lb, ub]", "[σ (nm), lb, ub]", "[Mol Frac, lb, ub]"]
    rows = []
    for i in range(best.shape[0]):
        rows.append([np.round([best[i, j], within[:, i, j].min(), within[:, i, j].max()],
                              decimals=3) for j in range(3)])
    return pd.DataFrame(rows, index=range(1, best.shape[0] + 1), columns=columns)


def export_results(spec, results_dir: str, processed_dir: str, name: str | None = None) -> list[str]:
    """Write the notebook's outputs; returns the paths written."""
    name = name or os.path.splitext(os.path.basename(spec.path))[0]
    os.makedirs(results_dir, exist_ok=True)
    os.makedirs(processed_dir, exist_ok=True)
    written = []
    for thresh in CONFIDENCE_THRESHOLDS:
        path = os.path.join(results_dir, f"{name}_{int(thresh * 100)}cutoff_model.csv")
        gaussian_table(spec, thresh).to_csv(path)
        written.append(path)
    path = os.path.join(results_dir, f"{name}_results.csv")
    results_frame(spec).to_csv(path)
    written.append(path)
    # Same layout the global-analysis notebook loads via `spectrum(<npz>)`.
    path = os.path.join(processed_dir, f"{name}.npz")
    np.savez(path, **spec.__getstate__())  # without in-memory caches
    written.append(path)
    return written
