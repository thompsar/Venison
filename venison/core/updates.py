"""Check GitHub for a newer Venison release."""

from __future__ import annotations

import json
import os
import re
import urllib.request

RELEASES_URL = "https://api.github.com/repos/thompsar/Venison/releases?per_page=100"
INSTALL_LINE = ("curl -LsSf https://raw.githubusercontent.com/thompsar/Venison/master/"
                "install.sh | sh")
_CALVER = re.compile(r"v?(\d+)\.(\d+)\.(\d+)")
_cache: dict[str, str | None] = {}  # one check per process, shared by all windows


def parse_version(text: str) -> tuple[int, int, int] | None:
    match = _CALVER.fullmatch(text.strip())
    return tuple(int(p) for p in match.groups()) if match else None


def fetch_latest_tag(timeout: float = 3.0) -> str | None:
    """Newest vYYYY.M.PATCH release tag. Not GitHub's "latest" label, which
    can point at other tags such as notebook-final."""
    request = urllib.request.Request(RELEASES_URL,
                                     headers={"Accept": "application/vnd.github+json"})
    with urllib.request.urlopen(request, timeout=timeout) as response:
        releases = json.load(response)
    tags = [r.get("tag_name", "") for r in releases if not r.get("draft") and not r.get("prerelease")]
    versions = [(parse_version(t), t) for t in tags if parse_version(t)]
    return max(versions)[1] if versions else None


def newer_release(current: str, fetch=fetch_latest_tag) -> str | None:
    """The latest release's version if newer than `current`, else None.

    Never raises: offline, rate-limited or odd tags (like `notebook-final`)
    just mean no notice. Set VENISON_NO_UPDATE_CHECK=1 to skip the request.
    """
    if os.environ.get("VENISON_NO_UPDATE_CHECK"):
        return None
    if current not in _cache:
        try:
            tag = fetch()
        except Exception:
            tag = None
        latest, mine = parse_version(tag or ""), parse_version(current)
        _cache[current] = (tag.lstrip("v") if latest and mine and latest > mine else None)
    return _cache[current]
