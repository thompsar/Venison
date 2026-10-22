#!/bin/sh
# Install or update Venison to the latest release (macOS or Linux):
#
#   curl -LsSf https://raw.githubusercontent.com/thompsar/Venison/main/install.sh | sh
#
# Run the same line again to update. Venison gets its own Python 3.12
# environment via uv (installed first if missing); nothing else on the system
# changes. Then type `venison` in a new Terminal window.
#
# VENISON_REF overrides what gets installed, e.g. VENISON_REF=tags/v2026.10.0
# or VENISON_REF=heads/main.
set -eu

REPO="thompsar/Venison"

say() { printf '%s\n' "$*"; }
fail() { say "Error: $*" >&2; exit 1; }

command -v curl >/dev/null 2>&1 || fail "curl is required"

# uv: use it if present, otherwise install it (it goes to ~/.local/bin).
PATH="$HOME/.local/bin:$HOME/.cargo/bin:$PATH"
if ! command -v uv >/dev/null 2>&1; then
    say "Installing uv (a small tool that manages Venison's Python)..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    command -v uv >/dev/null 2>&1 || fail "uv didn't install; see https://docs.astral.sh/uv/"
fi

# Which version: the newest vYYYY.M.PATCH release unless overridden. (Not
# GitHub's "latest" label, which can point at other tags like notebook-final.)
if [ -n "${VENISON_REF:-}" ]; then
    REF="$VENISON_REF"
else
    RELEASES=$(curl -fsSL "https://api.github.com/repos/$REPO/releases?per_page=100") \
        || fail "couldn't reach GitHub to find the latest release"
    TAG=$(printf '%s\n' "$RELEASES" \
          | sed -n 's/.*"tag_name"[[:space:]]*:[[:space:]]*"\(v[0-9]*\.[0-9]*\.[0-9]*\)".*/\1/p' \
          | sort -t. -k1.2,1n -k2,2n -k3,3n | tail -n 1)
    [ -n "$TAG" ] || fail "no Venison release found on GitHub yet"
    REF="tags/$TAG"
fi

# A release tarball, so git isn't needed (on a fresh Mac, git prompts to
# install Apple's developer tools).
URL="https://github.com/$REPO/archive/refs/$REF.tar.gz"
say "Installing Venison from $REF..."
uv tool install --quiet --python 3.12 --reinstall "venison @ $URL"

# Make sure ~/.local/bin is on PATH in future Terminal windows.
uv tool update-shell >/dev/null 2>&1 || true

VERSION=$("$(uv tool dir --bin)/venison" --version 2>/dev/null || true)
say ""
say "${VERSION:-Venison} is installed."
say "Open a new Terminal window and type:  venison"
say "To update later, run this same install line again."
