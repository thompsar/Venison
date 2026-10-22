# Venison: How the DEER sausage is made

A Python implementation of DEER waveform analysis (Tikhonov regularization and
Monte Carlo sum of Gaussians), with a desktop app for the full workflow.

## Getting started

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
python launch_venison.py
```

This opens Venison in its own window. Use `--no-webview` to open it in your
browser instead, `--no-browser` to only start the server, and `--port` to
change the port (default 5007).

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
venison/app.py  main window
launch_venison.py
tests/          pytest suite (runs the full workflow on data/raw/example_2dDEER)
```

`Venison.ipynb` and `src/` are legacy. `src/` only re-exports `venison.core`
so the notebook keeps running, and both will be removed.
