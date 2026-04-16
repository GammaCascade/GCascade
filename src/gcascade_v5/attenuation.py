from __future__ import annotations

import numpy as np

from . import legacy
from .config import njit


@njit(cache=True)
def _redshift_numba(spec: np.ndarray, left_idx: np.ndarray, weights: np.ndarray) -> np.ndarray:
    log_spec = np.empty(spec.shape[0], dtype=np.float64)
    for idx in range(spec.shape[0]):
        value = spec[idx]
        log_spec[idx] = np.log10(value) if value > 0.0 else -200.0

    shifted = np.empty(spec.shape[0], dtype=np.float64)
    for out_idx in range(spec.shape[0] - 1):
        left = int(left_idx[out_idx])
        weight = weights[out_idx]
        shifted[out_idx] = (1.0 - weight) * log_spec[left] + weight * log_spec[left + 1]
    shifted[spec.shape[0] - 1] = log_spec[spec.shape[0] - 1]

    out = np.empty(spec.shape[0], dtype=np.float64)
    for idx in range(spec.shape[0]):
        out[idx] = 10.0 ** shifted[idx] if shifted[idx] >= -199.0 else 0.0
    return out


@njit(cache=True)
def propagate_point_attenuation_numba(
    inj: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = inj.copy()
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_point_attenuation_window_numba(
    final_result: np.ndarray,
    attenuation_row: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    final_result = final_result * attenuation_row
    return _redshift_numba(final_result, left_idx, weights)


@njit(cache=True)
def propagate_diffuse_attenuation_numba(
    inj: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = np.zeros(inj.shape[0], dtype=np.float64)
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result + volume_norms[window_idx] * inj
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_diffuse_attenuation_window_numba(
    final_result: np.ndarray,
    source_term: np.ndarray,
    attenuation_row: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    final_result = final_result + source_term
    final_result = final_result * attenuation_row
    return _redshift_numba(final_result, left_idx, weights)


@njit(cache=True)
def propagate_evolving_attenuation_numba(
    inj2d: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = np.zeros(inj2d.shape[1], dtype=np.float64)
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result + volume_norms[window_idx] * inj2d[window_idx]
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_evolving_attenuation_window_numba(
    final_result: np.ndarray,
    source_term: np.ndarray,
    attenuation_row: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    final_result = final_result + source_term
    final_result = final_result * attenuation_row
    return _redshift_numba(final_result, left_idx, weights)


def propagate_point_attenuation(
    inj: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = np.asarray(inj, dtype=np.float64).copy()
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = legacy.RedshiftingCycle(
            final_result,
            legacy._build_z_windows(window_idx + 1)[-1],  # pragma: no cover - compatibility fallback
        )
    return final_result


def propagate_point_attenuation_fast(
    inj: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_cycle,
    progress=None,
) -> np.ndarray:
    final_result = np.asarray(inj, dtype=np.float64).copy()
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_cycle(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_diffuse_attenuation_fast(
    inj: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_cycle,
    progress=None,
) -> np.ndarray:
    final_result = np.zeros(len(legacy.energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_cycle(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_evolving_attenuation_fast(
    inj2d: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_cycle,
    progress=None,
) -> np.ndarray:
    final_result = np.zeros(len(legacy.energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj2d[window_idx]
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_cycle(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result
