#!/usr/bin/env python3
"""Tag and publish a Venison release, with preflight checks.

Everything the tag and release need is derived from `venison/__init__.py`,
CHANGELOG.md and git state, so there is nothing to mistype. The release notes
are the version's CHANGELOG section (they are public, on GitHub).

    python scripts/release.py --check     # preflight only, changes nothing
    python scripts/release.py             # preflight, show plan, confirm, release
    python scripts/release.py --yes       # skip the confirmation prompt

Run it from a clean, up-to-date `master` after the release PR has merged.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = "thompsar/Venison"
BRANCH = "master"
ROOT = Path(__file__).resolve().parent.parent
VERSION_FILE = ROOT / "venison" / "__init__.py"
CHANGELOG = ROOT / "CHANGELOG.md"
# CalVer YYYY.M.PATCH; a zero-padded month would normalize away under PEP 440.
CALVER = re.compile(r"^(20\d\d)\.([1-9]|1[0-2])\.(0|[1-9]\d*)$")


def run(*args: str, check: bool = True) -> str:
    result = subprocess.run(args, cwd=ROOT, capture_output=True, text=True)
    if check and result.returncode != 0:
        fail(f"`{' '.join(args)}` failed:\n{result.stderr.strip()}")
    return result.stdout.strip()


def ok(msg: str) -> None:
    print(f"  ok    {msg}")


def warn(msg: str) -> None:
    print(f"  warn  {msg}")


def fail(msg: str) -> None:
    print(f"  FAIL  {msg}")
    sys.exit(1)


def read_version() -> str:
    match = re.search(r'^__version__\s*=\s*"([^"]+)"', VERSION_FILE.read_text(), re.M)
    if not match:
        fail(f"no __version__ in {VERSION_FILE.relative_to(ROOT)}")
    version = match.group(1)
    if not CALVER.match(version):
        fail(f"__version__ {version!r} is not CalVer YYYY.M.PATCH with an unpadded month")
    ok(f"version {version}")
    return version


def check_git() -> None:
    branch = run("git", "rev-parse", "--abbrev-ref", "HEAD")
    if branch != BRANCH:
        fail(f"on {branch!r}; releases are tagged on {BRANCH!r}")
    if run("git", "status", "--porcelain", "--untracked-files=no"):
        fail("working tree has uncommitted changes to tracked files")
    run("git", "fetch", "--quiet", "origin", BRANCH, "--tags")
    local, remote = run("git", "rev-parse", "HEAD"), run("git", "rev-parse", f"origin/{BRANCH}")
    if local != remote:
        fail(f"{BRANCH} is not in sync with origin/{BRANCH}; pull or push first")
    ok(f"on {BRANCH}, clean, in sync with origin")


def check_tag_free(tag: str) -> None:
    if run("git", "tag", "-l", tag):
        fail(f"tag {tag} already exists locally")
    if run("git", "ls-remote", "--tags", "origin", f"refs/tags/{tag}"):
        fail(f"tag {tag} already exists on origin")
    ok(f"tag {tag} is free")


def changelog_section(version: str) -> str:
    text = CHANGELOG.read_text()
    match = re.search(rf"^## {re.escape(version)}\s*$(.*?)(?=^## |\Z)", text, re.M | re.S)
    if not match or not match.group(1).strip():
        fail(f"CHANGELOG.md has no section for {version}")
    ok("CHANGELOG.md has a section for this version")
    return match.group(1).strip()


def check_progression(version: str) -> None:
    year, month, patch = map(int, version.split("."))
    tags = [t for t in run("git", "tag", "-l", "v*").split() if CALVER.match(t[1:])]
    if not tags:
        return
    prev = max(tags, key=lambda t: tuple(map(int, t[1:].split("."))))
    py, pm, pp = map(int, prev[1:].split("."))
    if (year, month) == (py, pm) and patch != pp + 1:
        warn(f"previous is {prev}; expected PATCH {pp + 1}")
    elif (year, month) != (py, pm) and patch != 0:
        warn(f"new month after {prev}; PATCH would normally reset to 0")
    else:
        ok(f"follows {prev}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--check", action="store_true", help="preflight only")
    parser.add_argument("--yes", action="store_true", help="don't ask for confirmation")
    args = parser.parse_args()

    print("Preflight")
    version = read_version()
    tag = f"v{version}"
    check_git()
    check_tag_free(tag)
    notes = changelog_section(version)
    check_progression(version)
    if args.check:
        print("\nPreflight passed; nothing changed (--check).")
        return 0

    head = run("git", "log", "-1", "--format=%h %s")
    print(f"\nWill tag {tag} on {head}, push the tag, and publish a GitHub Release on {REPO}.")
    if not args.yes and input("Proceed? [y/N] ").strip().lower() != "y":
        print("Cancelled.")
        return 1

    run("git", "tag", "-a", tag, "-m", f"Venison {version}")
    run("git", "push", "origin", tag)
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False) as f:
        f.write(notes + "\n")
    run("gh", "release", "create", tag, "--repo", REPO, "--verify-tag",
        "--title", f"Venison {version}", "--notes-file", f.name)
    described = run("git", "describe", "--tags", "--exact-match", "HEAD", check=False)
    if described != tag:
        warn(f"`git describe` on HEAD gives {described!r}, expected {tag}")
    print(f"\nReleased {tag}: https://github.com/{REPO}/releases/tag/{tag}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
