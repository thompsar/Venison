"""Update notice: version comparison, failure handling, and the install script."""

import subprocess

import pytest

from venison.core import updates


@pytest.fixture(autouse=True)
def fresh(monkeypatch):
    monkeypatch.delenv("VENISON_NO_UPDATE_CHECK", raising=False)
    monkeypatch.setattr(updates, "_cache", {})


@pytest.mark.parametrize("current, tag, expected", [
    ("2026.10.0", "v2026.10.1", "2026.10.1"),
    ("2026.10.3", "v2026.11.0", "2026.11.0"),
    ("2026.12.4", "v2027.1.0", "2027.1.0"),
    ("2026.10.0", "v2026.10.0", None),        # same
    ("2026.11.0", "v2026.10.9", None),        # running something newer
    ("2026.10.0", "notebook-final", None),    # not a version
    ("2026.10.0", None, None),
])
def test_newer_release(current, tag, expected):
    assert updates.newer_release(current, fetch=lambda: tag) == expected


def test_offline_means_no_notice():
    def offline():
        raise OSError("no network")
    assert updates.newer_release("2026.10.0", fetch=offline) is None


def test_opt_out(monkeypatch):
    monkeypatch.setenv("VENISON_NO_UPDATE_CHECK", "1")
    assert updates.newer_release("2026.10.0", fetch=lambda: "v2099.1.0") is None


def test_install_script_is_valid_sh():
    subprocess.run(["sh", "-n", "install.sh"], check=True)
