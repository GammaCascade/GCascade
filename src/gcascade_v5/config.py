from __future__ import annotations

import os
import sys

try:  # Optional acceleration
    from numba import njit

    _NUMBA_AVAILABLE = True
except Exception:  # pragma: no cover - fallback when numba is unavailable
    _NUMBA_AVAILABLE = False

    def njit(*args, **kwargs):  # type: ignore[override]
        def _decorator(func):
            return func

        return _decorator


PROGRESS_ENABLED = os.getenv("GCASCADE_PROGRESS", "1") != "0"
NUMBA_ENABLED = _NUMBA_AVAILABLE and os.getenv("GCASCADE_NUMBA", "1") != "0"


def set_progress(enabled: bool) -> None:
    global PROGRESS_ENABLED
    PROGRESS_ENABLED = bool(enabled)


def set_numba(enabled: bool) -> None:
    global NUMBA_ENABLED
    if enabled and not _NUMBA_AVAILABLE:
        raise RuntimeError("numba is not installed; install it or disable acceleration.")
    NUMBA_ENABLED = bool(enabled)


def is_numba_available() -> bool:
    return bool(_NUMBA_AVAILABLE)


def get_numba_enabled() -> bool:
    return bool(NUMBA_ENABLED)


def status(message: str) -> None:
    if PROGRESS_ENABLED:
        print(message, flush=True)


def progress_marks(total: int, n_marks: int = 10) -> set[int]:
    if total <= 0:
        return set()
    if total == 1:
        return {1}
    marks = {max(1, int(round(total * i / n_marks))) for i in range(1, n_marks + 1)}
    marks.add(total)
    return marks


class ProgressBar:
    def __init__(self, label: str, total: int) -> None:
        self.label = label
        self.total = max(1, int(total))
        self.enabled = PROGRESS_ENABLED
        self._last_percent = -1
        self._finished = False
        if self.enabled:
            self.update(0)

    def update(self, current: int) -> None:
        if not self.enabled or self._finished:
            return
        current_clamped = max(0, min(int(current), self.total))
        percent = int((100 * current_clamped) / self.total)
        if percent == self._last_percent and current_clamped not in (0, self.total):
            return

        width = 30
        filled = (width * percent) // 100
        bar = "#" * filled + "-" * (width - filled)
        sys.stdout.write(f"\r{self.label}: [{bar}] {percent:3d}%")
        sys.stdout.flush()
        self._last_percent = percent

        if current_clamped >= self.total:
            self.close()

    def close(self) -> None:
        if not self.enabled or self._finished:
            return
        sys.stdout.write("\n")
        sys.stdout.flush()
        self._finished = True
