# Installing Venison

These steps are for a Mac with Anaconda or Miniconda, which you likely
already have if you've run the Venison notebook. Venison gets its own conda
environment, so your existing setup, including the notebook, is left alone.

## 1. Get Venison

**Either** download it: open the
[latest release](https://github.com/thompsar/Venison/releases/latest), download
**Source code (zip)** under *Assets*, and unzip it somewhere you'll keep it
(for example `~/Applications/Venison`).

**Or**, if you use git:

```bash
git clone https://github.com/thompsar/Venison.git
```

## 2. Install (once)

Open Terminal, go to the Venison folder (type `cd `, drag the folder onto the
Terminal window, press Return), then run:

```bash
conda env create -f environment.yml
```

This creates an environment called `venison` with Python 3.12 and everything
Venison needs. It takes a few minutes and downloads a few hundred MB.

## 3. Run

```bash
conda activate venison
venison
```

Venison opens in its own window. To quit, close the window (or press
Control-C in Terminal).

Next time, just open Terminal and run those two lines again; you don't need
to be in the Venison folder.

## Your first analysis

1. On the **Project** tab, click **Create project** and give it a name.
   Projects are kept in `~/venison_projects` (you can choose another folder in
   the sidebar).
2. Drag your spectrum's `.DSC` and `.DTA` files into **Add spectra**, add an
   optional sample and condition, and click **Add**.
3. Work through the tabs on the left from top to bottom. A green dot marks
   each finished step, and your progress is saved as you go.
4. On **Results**, click **Export results**. Tables are written to the
   project's `results/` folder and a reloadable `.npz` to `processed/`.

To try it first, use the example spectrum in the Venison folder:
`data/raw/example_2dDEER.DSC` and `.DTA`.

## Updating

Download and unzip the new release (or `git pull` in your clone), then, from
the Venison folder:

```bash
conda env update -f environment.yml --prune
```

Your projects aren't affected: they live in their own folder, not inside
Venison's.

## If something goes wrong

- **`venison: command not found`**: run `conda activate venison` first.
- **`conda: command not found`**: open a new Terminal window. If it's still
  missing, install [Miniconda](https://docs.anaconda.com/miniconda/) and try
  again.
- **The window stays blank or doesn't open**: run `venison --no-webview` to
  use your web browser instead.
- **"Port 5007 is already in use"**: Venison is already running (look for its
  window), or run `venison --port 5008`.
- **Anything else**: open the **Console** tab (bottom of the left rail), copy
  what it shows, and include it with the Venison version (shown next to the
  name in the top left) when you
  [report the problem](https://github.com/thompsar/Venison/issues/new/choose).

## Without conda

With Python 3.11 or newer (e.g. from python.org), from the Venison folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
venison
```
