from __future__ import annotations

"""Runtime verbosity helpers used by long propagation calls."""

import os
import sys


PROGRESS_ENABLED = os.getenv("GCASCADE_PROGRESS", "1") != "0"


def set_progress(enabled: bool) -> None:
    """Globally enable or disable status messages and progress bars."""
    global PROGRESS_ENABLED
    PROGRESS_ENABLED = bool(enabled)


def status(message: str) -> None:
    """Print one status line when runtime progress output is enabled."""
    if PROGRESS_ENABLED:
        print(message, flush=True)


def progress_marks(total: int, n_marks: int = 10) -> set[int]:
    """Return a sparse set of progress checkpoints for long build operations."""
    if total <= 0:
        return set()
    if total == 1:
        return {1}
    marks = {max(1, int(round(total * i / n_marks))) for i in range(1, n_marks + 1)}
    marks.add(total)
    return marks


class ProgressBar:
    """Lightweight terminal progress bar for cascade runs."""

    def __init__(self, label: str, total: int) -> None:
        """Create a progress bar with a short label and a known number of steps."""
        self.label = label
        self.total = max(1, int(total))
        self.enabled = PROGRESS_ENABLED
        self._last_percent = -1
        self._finished = False
        if self.enabled:
            self.update(0)

    def update(self, current: int) -> None:
        """Advance the bar to the current completed-step count."""
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
        """Finish the bar and move to the next output line."""
        if not self.enabled or self._finished:
            return
        sys.stdout.write("\n")
        sys.stdout.flush()
        self._finished = True
