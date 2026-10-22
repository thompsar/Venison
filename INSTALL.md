# Installing Venison

These steps are for a Mac. Venison runs in its own Python environment, so
whatever you already use (for example to run the Venison notebook) is left
alone.

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

Open Terminal and go to the Venison folder: type `cd ` (with a space), drag
the folder onto the Terminal window, and press Return.

Then check whether you have conda (it comes with Anaconda, Miniconda and
Miniforge):

```bash
conda --version
```

### If that printed a version: use conda

```bash
conda env create -f environment.yml
```

This creates an environment called `venison` with Python 3.12 and everything
Venison needs (a few minutes, a few hundred MB). To run Venison:

```bash
conda activate venison
venison
```

### If it said `command not found`: use uv

[uv](https://docs.astral.sh/uv/) is a small installer that brings its own
Python, so it doesn't matter which Python you have, or whether you have one.
Install uv:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

Open a **new** Terminal window (so it finds `uv`), go to the Venison folder
again, and install Venison:

```bash
uv tool install --python 3.12 .
```

To run Venison:

```bash
venison
```

## 3. Run

Venison opens in its own window. To quit, close the window (or press
Control-C in Terminal). Next time, open Terminal and run the same command(s)
as above; you don't need to be in the Venison folder.

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

Download and unzip the new release (or `git pull` in your clone), go to the
Venison folder in Terminal, and run the line for your install:

```bash
conda env update -f environment.yml --prune      # conda
uv tool install --python 3.12 --reinstall .      # uv
```

Your projects aren't affected: they live in their own folder, not inside
Venison's.

## If something goes wrong

- **`venison: command not found`**: with conda, run `conda activate venison`
  first. With uv, open a new Terminal window; if it's still missing, run
  `uv tool update-shell` and open another one.
- **The window stays blank or doesn't open**: run `venison --no-webview` to
  use your web browser instead.
- **"Port 5007 is already in use"**: Venison is already running (look for its
  window), or run `venison --port 5008`.
- **Anything else**: open the **Console** tab (bottom of the left rail), copy
  what it shows, and include it with the Venison version (shown next to the
  name in the top left) when you
  [report the problem](https://github.com/thompsar/Venison/issues/new/choose).

## Other ways to install

Any Python 3.11 or newer works with a virtual environment, from the Venison
folder:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install .
venison
```

(Then run `source .venv/bin/activate` from the Venison folder before
`venison` each time.)
