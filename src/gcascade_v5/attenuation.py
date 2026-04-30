from __future__ import annotations

"""Photon-only transport helpers used when pair production is treated as absorption only."""

import numpy as np

from .physics import energies


def propagate_point_attenuation(
    inj: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_step,
) -> np.ndarray:
    """Propagate a point-source photon spectrum with absorption and cosmological redshift."""
    final_result = np.asarray(inj, dtype=np.float64).copy()
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_step(final_result, window_idx)
    return final_result


def propagate_point_attenuation_fast(
    inj: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_step,
    progress=None,
) -> np.ndarray:
    """Propagate a point-source attenuation run while updating the progress bar."""
    final_result = np.asarray(inj, dtype=np.float64).copy()
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_step(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_diffuse_attenuation_fast(
    inj: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_step,
    progress=None,
) -> np.ndarray:
    """Propagate a diffuse photon emissivity while folding in absorption and redshift."""
    final_result = np.zeros(len(energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_step(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_evolving_attenuation_fast(
    inj2d: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    attenuation_vectors: np.ndarray,
    redshift_step,
    progress=None,
) -> np.ndarray:
    """Propagate a redshift-evolving photon emissivity without secondary production."""
    final_result = np.zeros(len(energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj2d[window_idx]
        final_result = final_result * attenuation_vectors[window_idx]
        final_result = redshift_step(final_result, window_idx)
        if progress is not None:
            progress.update(display_idx)
    return final_result

