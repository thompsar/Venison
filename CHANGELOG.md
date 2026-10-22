# Changelog

What changed for people using Venison, newest first. Versions are dates
(`YYYY.M.PATCH`); see RELEASING.md.

## 2026.10.0

The first version of Venison as a desktop app. The analysis is the same as
the notebook's; the notebook itself remains available as the
[`notebook-final`](https://github.com/thompsar/Venison/releases/tag/notebook-final)
release.

**New**
- A desktop app (`venison`) with one tab per step: Raw Signal, Background,
  Tikhonov, Validation, Stability, Gaussians, Monte Carlo and Results.
- Projects: drop `.DSC`/`.DTA` pairs into a project, label each spectrum with
  a sample and condition, and pick which one to analyze from the top bar.
  Every step is saved as you go and restored when you reopen a spectrum.
- Compare tab: overlay the project's spectra, waveforms with their fits and
  distance distributions, from either the Tikhonov or the Gaussian analysis.
- The Background tab shows the fitted parameters and, for the Fractal model,
  a 3D fit over the same range for comparison.
- Validation explains what it does, takes its settings in µs, can hold the
  Fractal dimension fixed across trials (smooth RMSD curves instead of jumps
  between bands), and warns when the Fractal fit falls back to 3D or the RMSD
  curve has no well.
- Sliders update the plots while you drag, as in the notebook; the Stability
  cutoff now does too, fast enough to follow the slider.
- Plot colors chosen to stay distinguishable with common color blindness.

**Fixed (compared with the notebook)**
- The validation minimum could be taken from the wrong position in the sweep.
- Choosing a validated background left the exported background unchanged.
- The best Gaussian model could be one that was never fitted.
- Monte Carlo failed when the fitted models didn't start at one Gaussian.
- Spectra in folders whose names contain a dot failed to load.

**Changed**
- Exports are the same CSV tables and `.npz` file as before, written into the
  project's `results/` and `processed/` folders.
- Needs Python 3.11 or newer (see INSTALL.md for a separate environment that
  leaves an existing notebook setup alone).
