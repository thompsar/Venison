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

    Averaged from the slice mask rather than read from `spec.real`, which the
    phasing step overwrites with the phased trace. Same normalization as
    `spectrum.avg_slices`, without mutating the spectrum.
    """
    real = np.sum(spec.real2D[spec.slice_mask], axis=0)
    imag = np.sum(spec.imaginary2D[spec.slice_mask], axis=0)
    peak = np.max(real)
    return real / peak + 1j * imag / peak


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
    spec.fractal_fixed_d = None
    spec.fit_background(spec.bgmodel, spec.background_range)
    # The dimensionality validation can hold fixed (None unless Fractal fit).
    spec.bg_tab_fractal_d = (getattr(spec, "fractal_d", None)
                             if model_name == "Fractal" and not spec.bg_fell_back else None)
    return spec.background_correct()


@dataclass
class BackgroundSummary:
    model: str                  # model actually used ("3D" if a Fractal fit fell back)
    b: float
    c: float
    d: float                    # 3 for the 3D model
    rmsd: float                 # fit residual over the background range
    alt: "BackgroundSummary | None" = None   # 3D fit on the same range, for a Fractal fit
    alt_curve: np.ndarray | None = None


def summarize_background(spec, start, model_name) -> BackgroundSummary:
    """Parameters and fit quality of the background fit from `start` to the
    cutoff, plus (for Fractal) a 3D fit over the same range for comparison.

    Refits from scratch (milliseconds) rather than reading the spectrum, whose
    background may since have been replaced by a validation choice."""
    rng = slice(nearest_index(spec.raw_time, start), nearest_index(spec.raw_time, spec.cutoff) + 1)
    x, y = spec.raw_time[rng], spec.real[rng]

    def rmsd_of(curve):
        return float(np.sqrt(np.mean((y - curve[rng]) ** 2)))

    b3, c3 = optimize.curve_fit(bg3D, x, y)[0]
    curve3 = bg3D(spec.raw_time, b3, c3)
    three_d = BackgroundSummary("3D", b3, c3, 3.0, rmsd_of(curve3))
    if model_name != "Fractal":
        return three_d
    try:
        b, c, d = optimize.curve_fit(bgFractal, x, y, bounds=[[0, 0, 1], [np.inf, np.inf, 5]])[0]
    except RuntimeError:
        return three_d  # same fallback the analysis uses
    fractal = BackgroundSummary("Fractal", b, c, d, rmsd_of(bgFractal(spec.raw_time, b, c, d)))
    fractal.alt, fractal.alt_curve = three_d, curve3
    return fractal


# ── Tikhonov / L-curve ───────────────────────────────────────────────────────

@dataclass
class LCurve:
    residual_norm: np.ndarray
    solution_norm: np.ndarray
    loocv_index: int


def compute_lcurve(spec) -> LCurve:
    spec.generate_lcurve()
    # Validation later swaps spec.waveform; keep the one this L-curve fit.
    spec.lcurve_waveform = np.array(spec.waveform)
    spec.lcurve_residual = np.linalg.norm(spec.tikhonovfits - spec.waveform, axis=1) ** 2
    spec.lcurve_solution = np.linalg.norm(np.dot(spec.Lmatrix, spec.solutions.T).T, axis=1) ** 2
    # Carried over from the notebook, which flags this LOOCV estimate as
    # unverified.
    spec.lcurve_loocv = int(
        (spec.lcurve_residual / (np.sum(1 - spec.hMats, axis=1)) ** 2).argmin())
    return lcurve_of(spec)


def lcurve_of(spec) -> LCurve:
    """The L-curve last computed for this spectrum."""
    return LCurve(spec.lcurve_residual, spec.lcurve_solution, int(spec.lcurve_loocv))


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


def background_tab_fractal_d(spec) -> float | None:
    """Fractal dimension of the Background tab's fit, or None if that fit fell
    back to 3D. Analyses saved before it was recorded get it refit, without
    touching the spectrum's current background."""
    if not hasattr(spec, "bg_tab_fractal_d"):
        start = getattr(spec, "bg_tab_start", spec.background_start)
        fit_range = range(nearest_index(spec.raw_time, start),
                          nearest_index(spec.raw_time, spec.cutoff) + 1)
        try:
            spec.bg_tab_fractal_d = float(optimize.curve_fit(
                bgFractal, spec.raw_time[fit_range], spec.real[fit_range],
                bounds=[[0, 0, 1], [np.inf, np.inf, 5]])[0][2])
        except RuntimeError:
            spec.bg_tab_fractal_d = None
    return spec.bg_tab_fractal_d


VALIDATION_DEFAULTS = {"earliest_start": 0.2, "shortest_fit": 0.8, "n_trials": 100,
                       "hold_fractal_d": True}


def validate_background(spec, earliest_start: float = 0.2, shortest_fit: float = 0.8,
                        n_trials: int = 100, hold_fractal_d: bool = True) -> Validation:
    """Sweep the background start and score each trial by how well it rebuilds
    the raw echo.

    Trial starts run from `earliest_start` (µs) to the point that still
    leaves `shortest_fit` µs of trace for the background fit. With a Fractal
    background and `hold_fractal_d`, every trial keeps the dimensionality
    fitted on the Background tab; left free, it swings from trial to trial
    and makes the RMSD curve jump.
    """
    dt = float(spec.raw_time[1] - spec.raw_time[0])
    end_offset = max(1, int(round(shortest_fit / dt)))
    held_d = (background_tab_fractal_d(spec)
              if hold_fractal_d and spec.bgmodel_name == "Fractal" else None)
    spec.fractal_fixed_d = held_d
    try:
        spec.validate_background(background_start=earliest_start, n_backgrounds=n_trials,
                                 end_offset=end_offset)
    finally:
        spec.fractal_fixed_d = None
    spec.validation_settings = {"earliest_start": earliest_start, "shortest_fit": shortest_fit,
                                "n_trials": n_trials, "hold_fractal_d": hold_fractal_d,
                                "held_d": held_d,
                                "fallbacks": int(np.sum(spec.validation_fell_back))}
    wr = spec.waveform_range

    amplitudes = [(1.0 - b[wr][0]) / b[wr][0] for b in spec.backgrounds]
    reconstructed = []
    for i, b in enumerate(spec.backgrounds):
        fit = np.dot(spec.kernel, spec.validation_solutions[i])
        reconstructed.append(fit / fit.max() * amplitudes[i] * b[wr] + b[wr])
    rmsd = np.sqrt(np.mean((spec.real[wr] - np.array(reconstructed)) ** 2, axis=1))
    # DOF is approximate; it only needs to find the edges of the well.
    min_index, (lo, hi) = suggest_range(rmsd, dof=len(spec.waveform) - 3)

    spec.validation_rmsd = rmsd
    spec.validation_min_index = min_index
    spec.validation_well = np.array([lo, hi])
    return validation_of(spec)


def suggest_range(rmsd, dof: int, stat_cutoff: float = 0.6) -> tuple[int, tuple[int, int]]:
    """Suggested background start (trial index) and range (inclusive indices)."""
    # The suggested background is the bottom of the first RMSD well: earlier
    # starts fold dipolar signal into the background, while very late starts
    # score well only because Tikhonov soaks up the background error as
    # long-distance P(r).
    # NOTE(notebook): the notebook always used the first well. A minimum on the
    # very first trial is the edge of the sweep, not a well (the curve rises
    # from the start), so move on to the next well in that case.
    peaks = [int(p) for p in argrelextrema(rmsd, np.greater)[0]] + [len(rmsd) - 1]
    start = 0
    for peak in peaks:
        left = np.where(rmsd[start:peak] >= rmsd[peak])[0]
        well_l = start + (int(left[-1]) if left.size else 0)
        min_index = _min_in_range(rmsd, well_l, peak)
        if min_index > 0:
            break
        start = peak

    # The range: neighbours statistically indistinguishable from the minimum.
    ok = np.round(ftest.cdf(rmsd / rmsd[min_index], dof, dof), decimals=2) <= stat_cutoff
    # NOTE(notebook): the notebook took the first run of acceptable trials,
    # which need not contain the minimum; grow outward from the minimum instead.
    lo = hi = min_index
    while lo > 0 and ok[lo - 1]:
        lo -= 1
    while hi < len(rmsd) - 1 and ok[hi + 1]:
        hi += 1
    return min_index, (lo, hi)


def validation_of(spec) -> Validation:
    """The background validation last run for this spectrum."""
    return Validation(
        fit_points=np.asarray(spec.background_fit_points),
        rmsd=spec.validation_rmsd,
        min_index=int(spec.validation_min_index),
        well=tuple(int(v) for v in spec.validation_well),
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


def prepare_regularizer(spec):
    """QR factors of [K; αL] for the spectrum's current kernel and alpha,
    computed once (~0.3 s) and cached on the spectrum."""
    key = (float(spec.alpha), spec.kernel.shape)
    cache = getattr(spec, "_regularize_qr", None)
    if cache is None or cache[0] != key:
        from scipy.linalg import qr
        q, r = qr(np.concatenate([spec.kernel, spec.alpha * spec.Lmatrix]), mode="economic")
        cache = spec._regularize_qr = (key, q, r)
    return cache[1], cache[2]


def regularize_fast(spec, waveform, warm=None) -> np.ndarray:
    """Same result as `spec.waveform_regularize`, about 3× faster when called
    repeatedly at one alpha.

    The system [K; αL] is fixed for a given kernel and alpha, so it is QR
    factored once: ‖Ax − b‖² = ‖Rx − Qᵀb‖² + const, and the NNLS runs on the
    square R instead of the tall A. The factorization is cached on the
    spectrum (underscore attributes aren't saved). `warm`, a solution to a
    nearby problem, often lets it skip the NNLS iterations entirely.
    """
    q, r = prepare_regularizer(spec)
    rhs = q.T @ np.concatenate([waveform, np.zeros(spec.Lmatrix.shape[0])])
    if warm is not None and warm.shape == (r.shape[1],) and np.any(warm > 0):
        # Warm start: guess the support from a nearby solution, solve plain
        # least squares on it, and accept only if it satisfies the NNLS
        # optimality (KKT) conditions, in which case it IS the NNLS solution.
        support = warm > 0
        xs = np.linalg.lstsq(r[:, support], rhs, rcond=None)[0]
        if np.all(xs > 0):
            x = np.zeros(r.shape[1])
            x[support] = xs
            gradient = r.T @ (rhs - r @ x)
            if gradient[~support].max(initial=0.0) <= 1e-10 * np.abs(r.T @ rhs).max():
                return x
    return optimize.nnls(r, rhs)[0]


def stability_cutoff(spec, index: int) -> dict:
    """Remove the unreliable long-distance tail and re-regularize, without
    changing the spectrum (safe to run on a worker thread). Returns the
    attributes `apply_stability_cutoff` sets."""
    mask = np.zeros_like(spec.r)
    mask[index:] = 1
    unstable = np.dot(spec.kernel, spec.solution * mask)
    stable_waveform = (spec.waveform - unstable) / (1.0 - unstable.max())
    # Neighbouring cutoffs have nearly the same solution: warm-start from it.
    solution = regularize_fast(spec, stable_waveform, warm=getattr(spec, "stable_solution", None))
    solution = solution / np.trapezoid(solution)
    fit = np.dot(spec.kernel, solution)
    return {"instability_index": int(index), "unstable_waveform": unstable,
            "stable_waveform": stable_waveform, "stable_solution": solution,
            "stable_fit": fit / fit[0]}


def apply_stability_cutoff(spec, index: int) -> np.ndarray:
    """Remove the unreliable long-distance tail and re-regularize.

    Returns the unfiltered solution normalized over the stable region, for
    plotting against `spec.stable_solution`.
    """
    for key, value in stability_cutoff(spec, index).items():
        setattr(spec, key, value)
    return renormalized_solution(spec, index)


def renormalized_solution(spec, index: int) -> np.ndarray:
    """The unfiltered solution normalized over the stable region (r < r[index])."""
    kept = np.where(np.arange(len(spec.r)) < index, spec.solution, 0.0)
    return spec.solution / np.trapezoid(kept)


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

    # Fitting state is stored on the spectrum so it is saved with it.
    def __init__(self, spec, ngauss: int = 5, reset: bool = True):
        self.spec = spec
        self.ngauss = ngauss
        if reset or getattr(spec, "gauss_fitted", None) is None:
            spec.gaussian_model_init(ngauss)
            spec.gauss_fitted = np.zeros(ngauss, dtype=bool)
            spec.gauss_optimized = np.zeros(ngauss, dtype=bool)
            spec.gauss_bics = np.full(ngauss, np.inf)
            spec.gauss_rmsds = np.zeros(ngauss)
        else:
            self.ngauss = len(spec.gauss_fitted)

    fitted = property(lambda self: self.spec.gauss_fitted)
    optimized = property(lambda self: self.spec.gauss_optimized)
    bics = property(lambda self: self.spec.gauss_bics)
    rmsds = property(lambda self: self.spec.gauss_rmsds)

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
