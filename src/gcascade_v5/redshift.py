from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d

from . import legacy
from .config import njit


def diffuse_distances_index(x: float) -> int:
    return int(np.argmin(np.abs(legacy.diffuseDistances - float(x))))


def prepare_log_spectrum(spec: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        log_spec = np.log10(spec)
    return np.where(np.isfinite(log_spec), log_spec, -200.0)


@njit(cache=True)
def _apply_redshift_numba(
    log_spec: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    n_energies = log_spec.shape[0]
    shifted = np.empty(n_energies, dtype=np.float64)
    for out_idx in range(n_energies - 1):
        left = int(left_idx[out_idx])
        weight = weights[out_idx]
        shifted[out_idx] = (1.0 - weight) * log_spec[left] + weight * log_spec[left + 1]
    shifted[n_energies - 1] = log_spec[n_energies - 1]

    out = np.empty(n_energies, dtype=np.float64)
    for idx in range(n_energies):
        out[idx] = 10.0 ** shifted[idx] if shifted[idx] >= -199.0 else 0.0
    return out


def redshift_cycle(
    spec: np.ndarray,
    *,
    left_idx: np.ndarray,
    weights: np.ndarray,
    use_numba: bool,
) -> np.ndarray:
    log_spec = prepare_log_spectrum(np.asarray(spec, dtype=np.float64))
    if use_numba:
        return _apply_redshift_numba(log_spec, left_idx, weights)

    shifted = np.empty_like(log_spec)
    shifted[:-1] = (
        (1.0 - weights[:-1]) * log_spec[left_idx[:-1].astype(np.int64, copy=False)]
        + weights[:-1] * log_spec[left_idx[:-1].astype(np.int64, copy=False) + 1]
    )
    shifted[-1] = log_spec[-1]
    return np.where(shifted >= -199.0, np.power(10.0, shifted), 0.0)


def legacy_redshift_cycle_for_validation(inj_spectra: np.ndarray, z_array_local: np.ndarray) -> np.ndarray:
    stretched_energies = legacy.energies * ((1.0 + z_array_local[0]) / (1.0 + z_array_local[-1]))
    logfunc = prepare_log_spectrum(np.asarray(inj_spectra, dtype=np.float64))
    interp = interp1d(
        legacy.energies,
        logfunc,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True,
    )
    stretched = np.empty_like(logfunc)
    stretched[:-1] = interp(stretched_energies[:-1])
    stretched[-1] = interp(legacy.energies[-1])
    return np.where(stretched >= -199.0, np.power(10.0, stretched), 0.0)
