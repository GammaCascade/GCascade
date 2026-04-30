from __future__ import annotations

"""Utilities for shifting spectra between adjacent redshift windows."""

import numpy as np
from scipy.interpolate import interp1d

from .physics import diffuseDistances, energies


def diffuse_distances_index(x: float) -> int:
    """Return the precomputed redshift-window index closest to a requested z value."""
    return int(np.argmin(np.abs(diffuseDistances - float(x))))


def prepare_log_spectrum(spec: np.ndarray) -> np.ndarray:
    """Convert a spectrum to log10 space while safely masking zeros and negatives."""
    with np.errstate(divide="ignore", invalid="ignore"):
        log_spec = np.log10(spec)
    return np.where(np.isfinite(log_spec), log_spec, -200.0)


def redshift_cycle(
    spec: np.ndarray,
    *,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    """Redshift one spectrum by interpolating it onto the next lower-energy grid."""
    log_spec = prepare_log_spectrum(np.asarray(spec, dtype=np.float64))
    shifted = np.empty_like(log_spec)
    left = left_idx[:-1].astype(np.int64, copy=False)
    shifted[:-1] = (1.0 - weights[:-1]) * log_spec[left] + weights[:-1] * log_spec[left + 1]
    shifted[-1] = log_spec[-1]
    return np.where(shifted >= -199.0, np.power(10.0, shifted), 0.0)


def reference_redshift_cycle(inj_spectra: np.ndarray, z_array_local: np.ndarray) -> np.ndarray:
    """Apply a direct interpolation redshift step used for validation tests."""
    stretched_energies = energies * ((1.0 + z_array_local[0]) / (1.0 + z_array_local[-1]))
    logfunc = prepare_log_spectrum(np.asarray(inj_spectra, dtype=np.float64))
    interp = interp1d(
        energies,
        logfunc,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True,
    )
    stretched = np.empty_like(logfunc)
    stretched[:-1] = interp(stretched_energies[:-1])
    stretched[-1] = interp(energies[-1])
    return np.where(stretched >= -199.0, np.power(10.0, stretched), 0.0)

