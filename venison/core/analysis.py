"""Analysis steps that used to live inside notebook cells.

Each function operates on a `spectrum` instance (mutating it the same way the
notebook did) and returns whatever the UI needs to draw. Nothing here imports
a plotting or widget library, so every step can be tested headless.

Numerics follow the notebook. Where the notebook had a bug, the fix is marked
with a `NOTE(notebook):` comment explaining the difference.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize
from scipy.signal import argrelextrema
from scipy.stats import f as ftest

from venison.core.background_models import bg3D, bgFractal
from venison.core.gradient_descent import gradient_descent

BG_MODELS = {"3D": bg3D, "Fractal": bgFractal}

NORMALIZERS = (
    "Gaussian Fit",
    "Gaussian Center, Data Maximum Amp",
    "Data Max Center, Data Maximum Amp",
)


def nearest_index(axis, value) -> int:
    return int(np.abs(np.asarray(axis) - value).argmin())


# ── Slices ───────────────────────────────────────────────────────────────────

def set_slice_mask(spec, mask) -> None:
    """Replace the slice mask and re-average. At least one slice must remain."""
    mask = np.asarray(mask, dtype=bool)
    if not mask.any():
        raise ValueError("At least one slice must be kept")
    spec.slice_mask = mask
    spec.avg_slices()


# ── Phase & cutoff ───────────────────────────────────────────────────────────

def complex_signal(spec) -> np.ndarray:
    """Unphased, slice-averaged complex signal.

    Re-averages from the slice mask rather than reading `spec.real`, which
    the phasing step overwrites with the phased trace.
    """
    spec.avg_slices()
    return spec.real + 1j * spec.imaginary


def auto_phase(signal, raw_time, cutoff) -> float:
    """Phase that minimizes the mean imaginary part between len/8 and cutoff."""
    stop = nearest_index(raw_time, cutoff)
    window = signal[len(signal) // 8:stop]
    return float(optimize.fmin(
        lambda x: abs((window * np.exp(1j * x)).imag.mean()), 1, disp=False)[0])


def phased(signal, phase) -> tuple[np.ndarray, np.ndarray]:
    rotated = signal * np.exp(1j * phase)
    peak = np.max(rotated.real)
    return rotated.real / peak, rotated.imag / peak


def apply_phase_and_cutoff(spec, phase, cutoff, *, rebuild_kernel=True) -> None:
    """Commit phase + cutoff to the spectrum and (re)build the kernel.

    The kernel only depends on the cutoff (tmax), so callers may skip the
    expensive rebuild when only the phase changed.
    """
    signal = complex_signal(spec)
    spec.phase = float(phase)
    spec.tmax = cutoff
    spec.cutoff = cutoff
    if rebuild_kernel or not hasattr(spec, "full_kernel"):
        spec.build_kernel()
    spec.real, spec.imag = phased(signal, spec.phase)


# ── Zero time & background ───────────────────────────────────────────────────

def gaussian(x, r, s, a):
    return a * np.exp(-((x - r) ** 2 / (2 * s ** 2)))


def symmetry_point(raw_time) -> int:
    """Index of the point mirroring raw_time[0] about t=0, +1 for slicing."""
    return int(np.argmax(raw_time >= -raw_time[0]) + 1)


def zero_window_bounds(spec) -> tuple[int, int]:
    start = 4 + int(np.argmax(spec.real == spec.real.max()))
    return start, max(start, round(len(spec.raw_time) / 30))


def background_start_options(spec) -> np.ndarray:
    lo = 4 + int(np.argmax(spec.real == spec.real.max()))
    hi = nearest_index(spec.raw_time, spec.cutoff) - 20
    return spec.raw_time[lo:hi]


def fit_zero_point(spec, window: int, normalizer: str = NORMALIZERS[0]) -> np.ndarray:
    """Fit a Gaussian to the first `window` points to locate t=0 and V(0)."""
    gaussfit = optimize.curve_fit(
        gaussian, spec.raw_time[1:window], spec.real[1:window])[0]
    spec.zeropoint, _, spec.zeroamp = gaussfit
    if normalizer == NORMALIZERS[1]:
        spec.zeroamp = spec.real[nearest_index(spec.raw_time, spec.zeropoint)]
    elif normalizer == NORMALIZERS[2]:
        spec.zeropoint = spec.raw_time[np.argmax(spec.real == spec.real.max())]
        spec.zeroamp = spec.real.max()
    return gaussfit


def fit_background(spec, background_start, model_name: str = "3D") -> np.ndarray:
    """Fit the background from `background_start` to cutoff and correct."""
    spec.bgmodel = BG_MODELS[model_name]
    spec.bgmodel_name = model_name
    spec.background_start = background_start
    spec.update_ranges()
    spec.fit_background(spec.bgmodel, spec.background_range)
    return spec.background_correct()


# ── Tikhonov / L-curve ───────────────────────────────────────────────────────

@dataclass
class LCurve:
    residual_norm: np.ndarray
    solution_norm: np.ndarray
    loocv_index: int


def compute_lcurve(spec) -> LCurve:
    spec.generate_lcurve()
    residual = np.linalg.norm(spec.tikhonovfits - spec.waveform, axis=1) ** 2
    solution = np.linalg.norm(np.dot(spec.Lmatrix, spec.solutions.T).T, axis=1) ** 2
    # Carried over from the notebook, which flags this LOOCV estimate as
    # unverified.
    loocv = int((residual / (np.sum(1 - spec.hMats, axis=1)) ** 2).argmin())
    return LCurve(residual, solution, loocv)


def select_alpha(spec, index: int) -> None:
    spec.alpha_index = int(index)
    spec.alpha = spec.alphas[index]
    spec.solution = spec.solutions[index]
    spec.tikhonovfit = spec.tikhonovfits[index]


# ── Background validation ────────────────────────────────────────────────────

@dataclass
class Validation:
    fit_points: np.ndarray      # background start time of each trial
    rmsd: np.ndarray            # reconstruction RMSD of each trial
    min_index: int              # RMSD minimum within the well
    well: tuple[int, int]       # default range-slider endpoints (inclusive)


def _min_in_range(rmsd, lo: int, hi: int) -> int:
    # NOTE(notebook): the notebook took argmin over a slice and used it as an
    # absolute index, which is wrong whenever the slice doesn't start at 0.
    hi = max(hi, lo + 1)
    return lo + int(rmsd[lo:hi].argmin())


def validate_background(spec, end_offset: int = 100) -> Validation:
    spec.validate_background(end_offset=end_offset)
    wr = spec.waveform_range

    amplitudes = [(1.0 - b[wr][0]) / b[wr][0] for b in spec.backgrounds]
    reconstructed = []
    for i, b in enumerate(spec.backgrounds):
        fit = np.dot(spec.kernel, spec.validation_solutions[i])
        reconstructed.append(fit / fit.max() * amplitudes[i] * b[wr] + b[wr])
    rmsd = np.sqrt(np.mean((spec.real[wr] - np.array(reconstructed)) ** 2, axis=1))

    try:
        maxima = int(argrelextrema(rmsd, np.greater)[0][0])
    except IndexError:
        maxima = len(rmsd) - 1
    try:
        well_l = int(np.where(rmsd[:maxima] >= rmsd[maxima])[0][-1])
    except IndexError:
        well_l = 0
    min_index = _min_in_range(rmsd, well_l, maxima)

    # DOF is approximate; it only needs to find the edges of the well.
    dof = len(spec.waveform) - 3
    stat = np.round(ftest.cdf(rmsd / rmsd[min_index], dof, dof), decimals=2)
    inside = np.where(stat <= 0.6)[0]
    gaps = np.where(np.diff(inside) > 1)[0]
    if gaps.size:
        inside = inside[:gaps[0] + 1]

    return Validation(
        fit_points=np.asarray(spec.background_fit_points),
        rmsd=rmsd,
        min_index=min_index,
        well=(int(inside[0]), int(inside[-1])),
    )


def validation_range_min(validation: Validation, lo: int, hi: int) -> int:
    """RMSD minimum within a user-chosen range (hi exclusive, as the notebook)."""
    return _min_in_range(validation.rmsd, lo, hi)


def select_validated_background(spec, index: int, range_=None) -> None:
    spec.validation_index = int(index)
    if range_ is not None:
        spec.validation_range = tuple(int(v) for v in range_)
    spec.background_start = spec.background_fit_points[index]
    spec.waveform = spec.waveforms[index]
    spec.tikhonovfit = spec.validation_tikhonovfits[index]
    spec.solution = spec.validation_solutions[index]
    # NOTE(notebook): the notebook left `background` at the pre-validation
    # fit, so the exported background didn't match the exported waveform.
    spec.background = spec.backgrounds[index]


def validation_envelope(spec, lo: int, hi: int) -> tuple[np.ndarray, np.ndarray]:
    sols = spec.validation_solutions[range(lo, max(hi, lo + 1))]
    return sols.min(axis=0), sols.max(axis=0)


# ── Stability (background supplementation) ───────────────────────────────────

def stability_limits(spec) -> tuple[float, float, int]:
    """Width and distance reliability limits for this trace length, plus the
    default instability index (midway between them).

    Jeschke, G. (2012). DEER distance measurements on proteins.
    Annu. Rev. Phys. Chem. 63, 419-446.
    """
    width_cutoff = 4 * (spec.t[-1] / 2) ** (1 / 3)
    distance_cutoff = 5 * (spec.t[-1] / 2) ** (1 / 3)
    idx = nearest_index(spec.r, width_cutoff + (distance_cutoff - width_cutoff) / 2)
    return width_cutoff, distance_cutoff, idx


def apply_stability_cutoff(spec, index: int) -> np.ndarray:
    """Remove the unreliable long-distance tail and re-regularize.

    Returns the unfiltered solution normalized over the stable region, for
    plotting against `spec.stable_solution`.
    """
    spec.instability_index = int(index)
    mask = np.zeros_like(spec.r)
    mask[index:] = 1
    spec.unstable_waveform = np.dot(spec.kernel, spec.solution * mask)
    spec.stable_waveform = (spec.waveform - spec.unstable_waveform) / (
        1.0 - spec.unstable_waveform.max())
    spec.stable_solution = spec.waveform_regularize(spec.stable_waveform)
    spec.stable_solution = spec.stable_solution / np.trapezoid(spec.stable_solution)
    spec.stable_fit = np.dot(spec.kernel, spec.stable_solution)
    spec.stable_fit = spec.stable_fit / spec.stable_fit[0]
    return spec.solution / np.trapezoid(spec.solution - spec.solution * mask)


# ── Gaussian modeling ────────────────────────────────────────────────────────

def gaussians(x, coeffs):
    """Sum of Gaussians; coeffs is (n, 3) with columns r, s, a."""
    peaks = (coeffs[:, 2] * np.exp(
        -((x[:, np.newaxis] - coeffs[:, 0]) ** 2 / (2 * coeffs[:, 1] ** 2)))).T
    return np.sum(peaks, axis=0)


def _curve_fit_gaussians(x, *coeffs):
    return gaussians(x, np.reshape(coeffs, (len(coeffs) // 3, 3)))


def bic(data, fit, n_model_params):
    k = n_model_params + 1
    n = len(data)
    return n * np.log(np.linalg.norm(data - fit) ** 2 / n) + k * np.log(n)


class GaussianModels:
    """Fitting state for the 1..n Gaussian models (the notebook's flag arrays)."""

    def __init__(self, spec, ngauss: int = 5):
        self.spec = spec
        self.ngauss = ngauss
        spec.gaussian_model_init(ngauss)
        self.fitted = np.zeros(ngauss, dtype=bool)
        self.optimized = np.zeros(ngauss, dtype=bool)
        self.bics = np.full(ngauss, np.inf)
        self.rmsds = np.zeros(ngauss)

    def model(self, m: int) -> np.ndarray:
        return self.spec.gaussian_model(m)

    def set_center(self, m: int, i: int, value: float) -> None:
        self.model(m)[i, 0] = value

    def fit_to_distribution(self, m: int) -> None:
        """'Initialize Gaussian(s)': fit model m to the stable distance distribution."""
        spec = self.spec
        initial = self.model(m).ravel().copy()
        # Reset widths and amplitudes so a previous bad fit can't trap this one.
        initial[1::3] = 0.1
        initial[2::3] = 0.1
        low = np.array([[-1.2 * r, -np.inf, -np.inf] for r in self.model(m)[:, 0]]).ravel()
        fit = optimize.curve_fit(_curve_fit_gaussians, spec.r, spec.stable_solution,
                                 p0=initial, bounds=(low, -low), method="trf")[0]
        self.model(m)[:] = np.abs(fit.reshape(m + 1, 3))
        self.fitted[m] = True
        self.optimized[m] = False
        self._score(m)

    def optimize_pending(self) -> list[int]:
        """'Fit Waveform': gradient-descent every fitted, not-yet-optimized model."""
        spec = self.spec
        done = []
        for m in np.where(self.fitted & ~self.optimized)[0]:
            self.model(m)[:] = gradient_descent(
                spec.stable_waveform, lambda x: spec.model_waveform(gaussians, x), self.model(m))
            self.optimized[m] = True
            self._score(m)
            done.append(int(m))
        return done

    def _score(self, m: int) -> None:
        wf = self.waveform(m)
        self.bics[m] = bic(self.spec.stable_waveform, wf, (m + 1) * 3)
        self.rmsds[m] = np.sqrt(np.mean((self.spec.stable_waveform - wf) ** 2)) * 1e2
        self.spec.best_model_id = self.best

    @property
    def best(self) -> int:
        return int(np.argmin(self.bics))

    def distribution(self, m: int) -> np.ndarray:
        curve = gaussians(self.spec.r, self.model(m))
        return curve / np.trapezoid(curve)

    def waveform(self, m: int) -> np.ndarray:
        return self.spec.model_waveform(gaussians, self.model(m))


# ── Monte Carlo error ────────────────────────────────────────────────────────

def default_rs_weights(spec, m: int) -> np.ndarray:
    """Larger r/s step sizes for small populations."""
    amps = spec.gaussian_model(m)[:, 2]
    return np.round(amps.max() / amps)


def step_weights(spec, rs_weights) -> np.ndarray:
    rs_weights = np.asarray(rs_weights, dtype=float)
    deltar = spec.r[1] - spec.r[0]
    return np.hstack((np.tile(rs_weights[:, None], (1, 2)) * deltar,
                      np.ones((len(rs_weights), 1)) * 0.01))


def run_monte_carlo(spec, m: int, rs_weights, iterations: int = 5000) -> None:
    spec.mc_model_id = int(m)
    spec.monte_carlo_error(data=spec.stable_waveform, origin=spec.gaussian_model(m),
                           weights=step_weights(spec, rs_weights), iterations=iterations)
    spec.norm_constants = np.array(
        [np.trapezoid(gaussians(spec.r, c), spec.r) for c in spec.model_landscape])


def mole_fraction_landscape(spec) -> np.ndarray:
    """Model landscape with amplitudes converted to normalized mole fractions."""
    scaled = spec.model_landscape.copy()
    scaled[:, :, 2] = (np.sqrt(2 * np.pi) * scaled[:, :, 1] * scaled[:, :, 2]
                       / spec.norm_constants[:, None])
    return scaled


# ── Results ──────────────────────────────────────────────────────────────────

def best_monte_carlo_model(spec, confidence: float = 0.95):
    """Best Monte Carlo model, its fit, and the min/max waveform envelope of
    every model under the confidence cutoff."""
    spec.best_model = spec.model_landscape[np.argmin(spec.landscape_statistics)]
    spec.best_model_fit = spec.model_waveform(gaussians, spec.best_model)
    solution = gaussians(spec.r, spec.best_model)
    spec.best_model_norm = np.trapezoid(solution, spec.r)
    spec.best_model_solution = solution / spec.best_model_norm
    variability = np.array([spec.model_waveform(gaussians, model) for model in
                            spec.model_landscape[spec.landscape_statistics < confidence]])
    return variability.min(axis=0), variability.max(axis=0)


def population_curves(spec) -> list[np.ndarray]:
    return [gaussians(spec.r, spec.best_model[p, None]) / spec.best_model_norm
            for p in range(spec.best_model.shape[0])]
