from __future__ import annotations

import numpy as np

from . import legacy
from .config import njit
from .redshift import redshift_cycle


_D_ENERGY = np.asarray(legacy.dEnergiesGamma, dtype=np.float64)


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
def _cascade_step_numba(
    final_result: np.ndarray,
    coeff_row: np.ndarray,
    cycle_row: np.ndarray,
    step_size: float,
    d_energy: np.ndarray,
) -> np.ndarray:
    n_energies = final_result.shape[0]
    attenuated = np.empty(n_energies, dtype=np.float64)
    row_coeffs = np.empty(n_energies, dtype=np.float64)
    result = np.empty(n_energies, dtype=np.float64)

    for idx in range(n_energies):
        attenuated[idx] = (coeff_row[idx] ** step_size) * final_result[idx]
        row_coeffs[idx] = final_result[idx] - attenuated[idx]
        result[idx] = attenuated[idx]

    row_coeffs[0] = row_coeffs[0] * d_energy[0]
    for idx in range(1, n_energies - 1):
        row_coeffs[idx] = row_coeffs[idx] * (d_energy[idx - 1] + d_energy[idx])
    row_coeffs[n_energies - 1] = row_coeffs[n_energies - 1] * d_energy[n_energies - 2]

    cursor = 0
    for row_idx in range(n_energies):
        contrib = row_coeffs[row_idx]
        for out_idx in range(row_idx + 1):
            result[out_idx] += contrib * cycle_row[cursor + out_idx]
        cursor += row_idx + 1
    return result


@njit(cache=True)
def propagate_point_cascade_numba(
    inj: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = inj.copy()
    for window_idx in range(z_max_index - 1, -1, -1):
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        for step_idx in range(start, stop):
            z_index = int(zreg_indices[step_idx])
            final_result = _cascade_step_numba(
                final_result,
                extinction_coeffs[z_index],
                cycle_packed[z_index],
                step_sizes[step_idx],
                _D_ENERGY,
            )
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_point_cascade_window_numba(
    final_result: np.ndarray,
    start: int,
    stop: int,
    step_sizes: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    for step_idx in range(start, stop):
        z_index = int(zreg_indices[step_idx])
        final_result = _cascade_step_numba(
            final_result,
            extinction_coeffs[z_index],
            cycle_packed[z_index],
            step_sizes[step_idx],
            _D_ENERGY,
        )
    return _redshift_numba(final_result, left_idx, weights)


@njit(cache=True)
def propagate_diffuse_cascade_numba(
    inj: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = np.zeros(inj.shape[0], dtype=np.float64)
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result + volume_norms[window_idx] * inj
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        for step_idx in range(start, stop):
            z_index = int(zreg_indices[step_idx])
            final_result = _cascade_step_numba(
                final_result,
                extinction_coeffs[z_index],
                cycle_packed[z_index],
                step_sizes[step_idx],
                _D_ENERGY,
            )
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_diffuse_cascade_window_numba(
    final_result: np.ndarray,
    source_term: np.ndarray,
    start: int,
    stop: int,
    step_sizes: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    final_result = final_result + source_term
    for step_idx in range(start, stop):
        z_index = int(zreg_indices[step_idx])
        final_result = _cascade_step_numba(
            final_result,
            extinction_coeffs[z_index],
            cycle_packed[z_index],
            step_sizes[step_idx],
            _D_ENERGY,
        )
    return _redshift_numba(final_result, left_idx, weights)


@njit(cache=True)
def propagate_evolving_cascade_numba(
    inj2d: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
) -> np.ndarray:
    final_result = np.zeros(inj2d.shape[1], dtype=np.float64)
    for window_idx in range(z_max_index - 1, -1, -1):
        final_result = final_result + volume_norms[window_idx] * inj2d[window_idx]
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        for step_idx in range(start, stop):
            z_index = int(zreg_indices[step_idx])
            final_result = _cascade_step_numba(
                final_result,
                extinction_coeffs[z_index],
                cycle_packed[z_index],
                step_sizes[step_idx],
                _D_ENERGY,
            )
        final_result = _redshift_numba(final_result, redshift_left_idx[window_idx], redshift_weights[window_idx])
    return final_result


@njit(cache=True)
def propagate_evolving_cascade_window_numba(
    final_result: np.ndarray,
    source_term: np.ndarray,
    start: int,
    stop: int,
    step_sizes: np.ndarray,
    zreg_indices: np.ndarray,
    extinction_coeffs: np.ndarray,
    cycle_packed: np.ndarray,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> np.ndarray:
    final_result = final_result + source_term
    for step_idx in range(start, stop):
        z_index = int(zreg_indices[step_idx])
        final_result = _cascade_step_numba(
            final_result,
            extinction_coeffs[z_index],
            cycle_packed[z_index],
            step_sizes[step_idx],
            _D_ENERGY,
        )
    return _redshift_numba(final_result, left_idx, weights)


def _cascade_step_python(
    final_result: np.ndarray,
    log_coeff_row: np.ndarray,
    weighted_cycle: np.ndarray,
    step_size: float,
) -> np.ndarray:
    attenuated = np.exp(log_coeff_row * step_size) * final_result
    delta = final_result - attenuated
    return attenuated + weighted_cycle @ delta


def propagate_point_cascade_fast(
    inj: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    log_extinction_coeffs: np.ndarray,
    get_weighted_cycle_slice,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    progress=None,
) -> np.ndarray:
    final_result = np.asarray(inj, dtype=np.float64).copy()
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        block_start = start
        while block_start < stop:
            z_index = int(zreg_indices[block_start])
            block_stop = block_start + 1
            while block_stop < stop and int(zreg_indices[block_stop]) == z_index:
                block_stop += 1
            weighted_cycle = get_weighted_cycle_slice(z_index)
            log_coeff_row = log_extinction_coeffs[z_index]
            for step_idx in range(block_start, block_stop):
                final_result = _cascade_step_python(
                    final_result,
                    log_coeff_row,
                    weighted_cycle,
                    float(step_sizes[step_idx]),
                )
            block_start = block_stop
        final_result = redshift_cycle(
            final_result,
            left_idx=redshift_left_idx[window_idx],
            weights=redshift_weights[window_idx],
            use_numba=False,
        )
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_diffuse_cascade_fast(
    inj: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    log_extinction_coeffs: np.ndarray,
    get_weighted_cycle_slice,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    progress=None,
) -> np.ndarray:
    final_result = np.zeros(len(legacy.energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        block_start = start
        while block_start < stop:
            z_index = int(zreg_indices[block_start])
            block_stop = block_start + 1
            while block_stop < stop and int(zreg_indices[block_stop]) == z_index:
                block_stop += 1
            weighted_cycle = get_weighted_cycle_slice(z_index)
            log_coeff_row = log_extinction_coeffs[z_index]
            for step_idx in range(block_start, block_stop):
                final_result = _cascade_step_python(
                    final_result,
                    log_coeff_row,
                    weighted_cycle,
                    float(step_sizes[step_idx]),
                )
            block_start = block_stop
        final_result = redshift_cycle(
            final_result,
            left_idx=redshift_left_idx[window_idx],
            weights=redshift_weights[window_idx],
            use_numba=False,
        )
        if progress is not None:
            progress.update(display_idx)
    return final_result


def propagate_evolving_cascade_fast(
    inj2d: np.ndarray,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    log_extinction_coeffs: np.ndarray,
    get_weighted_cycle_slice,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    progress=None,
) -> np.ndarray:
    final_result = np.zeros(len(legacy.energies), dtype=np.float64)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = final_result + volume_norms[window_idx] * inj2d[window_idx]
        start = int(row_ptr[window_idx])
        stop = int(row_ptr[window_idx + 1])
        block_start = start
        while block_start < stop:
            z_index = int(zreg_indices[block_start])
            block_stop = block_start + 1
            while block_stop < stop and int(zreg_indices[block_stop]) == z_index:
                block_stop += 1
            weighted_cycle = get_weighted_cycle_slice(z_index)
            log_coeff_row = log_extinction_coeffs[z_index]
            for step_idx in range(block_start, block_stop):
                final_result = _cascade_step_python(
                    final_result,
                    log_coeff_row,
                    weighted_cycle,
                    float(step_sizes[step_idx]),
                )
            block_start = block_stop
        final_result = redshift_cycle(
            final_result,
            left_idx=redshift_left_idx[window_idx],
            weights=redshift_weights[window_idx],
            use_numba=False,
        )
        if progress is not None:
            progress.update(display_idx)
    return final_result
