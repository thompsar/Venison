"""Tee stdout/stderr into a ring buffer so the Console tab can show it."""

from __future__ import annotations

import sys
import threading
from collections import deque
from datetime import datetime

_lines: deque[str] = deque(maxlen=2000)
_lock = threading.Lock()
_partial = {"out": "", "err": ""}
version = 0  # bumped on every new line; the console polls it


class _Tee:
    def __init__(self, stream, key):
        self._stream = stream
        self._key = key

    def write(self, text):
        global version
        try:
            self._stream.write(text)
        except Exception:
            pass
        with _lock:
            buf = _partial[self._key] + text
            *complete, _partial[self._key] = buf.split("\n")
            stamp = datetime.now().strftime("%H:%M:%S")
            for line in complete:
                if line.strip():
                    _lines.append(f"{stamp}  {line}")
                    version += 1
        return len(text)

    def flush(self):
        try:
            self._stream.flush()
        except Exception:
            pass

    def __getattr__(self, name):
        return getattr(self._stream, name)


def install() -> None:
    if not isinstance(sys.stdout, _Tee):
        sys.stdout = _Tee(sys.stdout, "out")
    if not isinstance(sys.stderr, _Tee):
        sys.stderr = _Tee(sys.stderr, "err")


def lines() -> list[str]:
    with _lock:
        return list(_lines)
