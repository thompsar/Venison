# Installing Venison

These steps are for a Mac (Linux works the same way). Venison runs in its own
Python environment, so whatever you already have, including a setup for the
Venison notebook, is left alone.

## Install

Open **Terminal** (in Applications → Utilities, or search for it with
Spotlight), paste this line, and press Return:

```bash
curl -LsSf https://raw.githubusercontent.com/thompsar/Venison/main/install.sh | sh
```

It installs the latest release of Venison, together with
[uv](https://docs.astral.sh/uv/), a small tool that gives Venison its own copy
of Python. It takes a minute or two the first time.

## Run

Open a **new** Terminal window and type:

```bash
venison
```

Venison opens in its own window. To quit, close the window (or press
Control-C in Terminal). Next time, just open Terminal and type `venison` again.

## Update

When a new version is out, Venison shows a notice at the top of its window.
Quit Venison and run the **same install line** again:

```bash
curl -LsSf https://raw.githubusercontent.com/thompsar/Venison/main/install.sh | sh
```

Your projects aren't affected: they live in their own folder
(`~/venison_projects` unless you chose another).

## Your first analysis

1. On the **Project** tab, click **Create project** and give it a name.
2. Drag your spectrum's `.DSC` and `.DTA` files into **Add spectra**, add an
   optional sample and condition, and click **Add**.
3. Work through the tabs on the left from top to bottom. A green dot marks
   each finished step, and your progress is saved as you go.
4. On **Results**, click **Export results**. Tables are written to the
   project's `results/` folder and a reloadable `.npz` to `processed/`.

To try it first, download the example spectrum's
[`.DSC`](https://github.com/thompsar/Venison/raw/main/data/raw/example_2dDEER.DSC)
and [`.DTA`](https://github.com/thompsar/Venison/raw/main/data/raw/example_2dDEER.DTA)
files.

## If something goes wrong

- **`venison: command not found`**: open a new Terminal window. If it's still
  missing, run `~/.local/bin/uv tool update-shell` and open another one.
- **The window stays blank or doesn't open**: run `venison --no-webview` to
  use your web browser instead.
- **"Port 5007 is already in use"**: Venison is already running (look for its
  window), or run `venison --port 5008`.
- **Anything else**: open the **Console** tab (bottom of the left rail), copy
  what it shows, and include it with the Venison version (shown next to the
  name in the top left) when you
  [report the problem](https://github.com/thompsar/Venison/issues/new/choose).

To remove Venison: `uv tool uninstall venison` (your projects stay).

## Other ways to install

**With conda** (Anaconda, Miniconda, Miniforge): download **Source code
(zip)** from the [latest release](https://github.com/thompsar/Venison/releases/latest),
unzip it, and from that folder in Terminal:

```bash
conda env create -f environment.yml     # once; to update: conda env update -f environment.yml --prune
conda activate venison
venison
```

**From a source folder** (a release zip or a git clone), with any Python 3.11
or newer:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install .
venison
```

The update notice can be turned off by setting `VENISON_NO_UPDATE_CHECK=1`.
