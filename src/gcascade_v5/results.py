from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any

import numpy as np


@dataclass(frozen=True)
class CascadeResult:
    """Two-species cascade output.

    The gamma spectrum is the public observable returned by the standard
    return path. The electron spectrum is an extragalactic boundary diagnostic,
    not a charged-particle flux at Earth.
    """

    gamma: np.ndarray
    electron: np.ndarray
    diagnostics: dict[str, float]
    metadata: dict[str, Any]

    def __array__(self, dtype: np.dtype | None = None) -> np.ndarray:
        """Allow NumPy to view the result as its gamma-ray spectrum."""
        return np.asarray(self.gamma, dtype=dtype)

    def with_scaled_spectra(self, factor: float) -> "CascadeResult":
        """Scale both spectra together, preserving the diagnostic bookkeeping."""
        scale = float(factor)
        return replace(
            self,
            gamma=np.asarray(self.gamma, dtype=np.float64) * scale,
            electron=np.asarray(self.electron, dtype=np.float64) * scale,
            metadata={**self.metadata, "spectrum_scale_factor": scale},
        )
