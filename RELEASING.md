# Releasing Venison

## Model

- `main` is the trunk; PRs merge into it without touching the version.
- `__version__` in `venison/__init__.py` is the last **released** version and
  changes only when cutting a release. `pyproject.toml` reads it, so that is
  the one edit.
- A release is an annotated tag `vYYYY.M.PATCH` on `main` plus a GitHub
  Release whose notes are that version's section of `CHANGELOG.md`.
  Collaborators install from these releases (see INSTALL.md).

## Version scheme: CalVer `YYYY.M.PATCH`

- Month not zero-padded (`2026.10.0`, `2027.1.0`), since PEP 440 normalizes
  leading zeros and source and installed metadata should match.
- New month: roll `YYYY.M` forward and reset PATCH to 0. Same month: PATCH + 1.

## Cutting a release

1. In a PR: set `__version__`, add a `## <version>` section at the top of
   `CHANGELOG.md` written for users (what's new, what's fixed, what changed in
   how they work), and merge it.
2. On an up-to-date `main`:

   ```bash
   python scripts/release.py --check    # preflight only
   python scripts/release.py            # tag, push the tag, publish the release
   ```

   The script refuses to run off `main`, with uncommitted changes, out of
   sync with `origin`, without a changelog section, or if the tag exists. It
   warns if the PATCH number looks off.

## Other tags

- `notebook-final`: the last notebook-only version, before the desktop app.

## Before tagging

The repository is public, and release downloads contain every tracked file.
Raw data stays out via `.gitignore` (only `data/raw/example_2dDEER*` is
tracked); keep unpublished sample names, data and personal notes out of
commits, the changelog and release notes.
