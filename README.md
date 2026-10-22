# Venison: How the DEER sausage is made

A Python implementation of DEER waveform analysis (Tikhonov regularization and
Monte Carlo sum of Gaussians), with a desktop app for the full workflow.

## Getting started

From the Venison folder, with conda:

```bash
conda env create -f environment.yml    # once
conda activate venison
venison
```

or, without conda, with [uv](https://docs.astral.sh/uv/) (it brings its own Python):

```bash
uv tool install --python 3.12 .        # once
venison
```

Step-by-step instructions for a Mac, updating and troubleshooting are in
[INSTALL.md](INSTALL.md). `venison --help` lists the options (for example
`--no-webview` to use your browser instead of a window).

The last notebook-only version is the
[`notebook-final`](https://github.com/thompsar/Venison/releases/tag/notebook-final)
release. What's changed since is in [CHANGELOG.md](CHANGELOG.md).

## Workflow

Each tab on the left rail is one step. A green dot marks a step that is done
for the current spectrum. Changing an earlier step clears the steps after it.

| Tab | What it does |
|---|---|
| Project | Open or create a project, add `.DSC`/`.DTA` pairs, label and pick spectra |
| Raw Signal | Drop bad slices, phase the signal, choose the cutoff, build the kernel |
| Background | Fit zero time and V(t=0), choose the background model and start |
| Tikhonov | Compute the L-curve and choose the regularization parameter |
| Validation | Refit the background from many start points and pick the best |
| Stability | Remove the unreliable long-distance tail |
| Gaussians | Fit 1–5 Gaussian models and compare them by BIC |
| Monte Carlo | Map the error surface of the chosen model |
| Results | Best fit, confidence bounds, and export |
| Compare | Overlay the project's spectra: waveforms, fits and distance distributions |
| Console | Log output |

## Projects

Projects live in `~/venison_projects` by default. You can change the folder on
the Project tab, or set `VENISON_PROJECTS_ROOT`. A project holds any number of
spectra, analyzed one at a time:

```
<project>/
  raw_data/            ingested .DSC/.DTA files and spectra.json (sample, condition, notes)
  workspace/           saved analysis state, so a spectrum reopens where you left it
  results/             exported CSV tables
  processed/           exported .npz spectra (reload with venison.core.spectrum.spectrum)
  analysis_log.json    a record of every analysis step and its settings
```

## Code layout

```
venison/core/   analysis code with no UI dependencies (spectrum, analysis steps, export, projects)
venison/ui/     one module per tab
venison/app.py  main window; venison/launcher.py starts it
launch_venison.py   same as the `venison` command, from a source checkout
tests/          pytest suite (runs the full workflow on data/raw/example_2dDEER)
```

`Venison.ipynb` and `src/` are legacy. `src/` only re-exports `venison.core`
so the notebook keeps running, and both will be removed.
