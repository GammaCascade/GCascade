from __future__ import annotations

"""Two-species gamma/electron transport operators for the modern cascade."""

from collections.abc import Callable
from typing import Any

import numpy as np

from .physics import Mpc, c, echarge, elecmass, energies, mu0, sigmaTe, zReg, dEnergiesGamma
from .redshift import redshift_cycle
from .results import CascadeResult


TRAPEZOID_WEIGHTS = np.empty(len(energies), dtype=np.float64)
TRAPEZOID_WEIGHTS[0] = dEnergiesGamma[0]
TRAPEZOID_WEIGHTS[-1] = dEnergiesGamma[-1]
TRAPEZOID_WEIGHTS[1:-1] = dEnergiesGamma[:-1] + dEnergiesGamma[1:]
SECONDS_PER_MPC = Mpc / (c * 1.0e5)
DEFAULT_MAX_ICS_TAU_PER_SUBSTEP = 0.05
DEFAULT_MAX_ICS_SUBSTEPS = 32
DEFAULT_CEL_LOG_BIN_THRESHOLD = 0.25


WeightedKernelGetter = Callable[[int], np.ndarray]
CelDataGetter = Callable[[int, float], tuple[np.ndarray, np.ndarray]]


def spectrum_energy(spec: np.ndarray) -> float:
    """Integrate a differential spectrum over the GCascade energy grid."""
    arr = np.asarray(spec, dtype=np.float64)
    return float(np.sum(arr * energies * TRAPEZOID_WEIGHTS))


def zero_diagnostics() -> dict[str, float]:
    """Create the bookkeeping dictionary used to track the cascade energy budget."""
    return {
        "initial_gamma_energy": 0.0,
        "initial_electron_energy": 0.0,
        "final_gamma_energy": 0.0,
        "final_electron_energy": 0.0,
        "pp_absorbed_energy": 0.0,
        "ics_scattered_electron_energy": 0.0,
        "synchrotron_energy_lost": 0.0,
        "below_grid_energy_lost": 0.0,
        "redshift_energy_lost": 0.0,
        "cel_scattered_electron_energy": 0.0,
        "energy_residual": 0.0,
        "relative_energy_residual": 0.0,
    }


def synchrotron_loss_rates(b_field_gauss: float, gamma: float, z_value: float) -> np.ndarray:
    """Return dE/dt for synchrotron cooling in eV s^-1 on the electron grid."""
    b_tesla = float(b_field_gauss) / 10000.0
    if b_tesla == 0.0:
        return np.zeros(len(energies), dtype=np.float64)
    z_factor = (1.0 + float(z_value)) ** float(gamma)
    synch_prefactor = ((b_tesla * z_factor) ** 2) / (2.0 * mu0)
    gamma_factor = np.sqrt((energies * 1.0e9) ** 2 - elecmass**2) / elecmass
    return (
        # The synchrotron power is first computed in J/s and then converted
        # to eV/s so it can be compared consistently to the stored ICS dE/dt tables.
        (1.0 / echarge)
        * (4.0 / 3.0)
        * sigmaTe
        * synch_prefactor
        * c
        * 1000.0
        * np.power(gamma_factor, 2.0)
    )


def _deposit_to_log_grid(amount: float, input_weight: float, target_energy: float, out: np.ndarray) -> float:
    """Deposit particles onto the logarithmic energy grid and return any sub-grid energy."""
    if amount <= 0.0:
        return 0.0
    particle_count = float(amount) * float(input_weight)
    if target_energy < energies[0]:
        return float(particle_count * target_energy)
    if target_energy >= energies[-1]:
        out[-1] += particle_count / TRAPEZOID_WEIGHTS[-1]
        return 0.0

    idx = int(np.searchsorted(energies, target_energy, side="right") - 1)
    idx = max(0, min(idx, len(energies) - 2))
    left = energies[idx]
    right = energies[idx + 1]
    weight = (np.log(target_energy) - np.log(left)) / (np.log(right) - np.log(left))
    weight = float(np.clip(weight, 0.0, 1.0))
    out[idx] += particle_count * (1.0 - weight) / TRAPEZOID_WEIGHTS[idx]
    out[idx + 1] += particle_count * weight / TRAPEZOID_WEIGHTS[idx + 1]
    return 0.0


def _apply_synchrotron_losses(
    electron: np.ndarray,
    *,
    d_edt_sync_row: np.ndarray,
    sub_step: float,
    diagnostics: dict[str, float],
) -> np.ndarray:
    """Cool the electron spectrum continuously by synchrotron emission over one sub-step."""
    if float(sub_step) <= 0.0 or not np.any(electron > 0.0) or not np.any(d_edt_sync_row > 0.0):
        return electron

    dt_seconds = float(sub_step) * SECONDS_PER_MPC
    if dt_seconds <= 0.0:
        return electron

    energy_losses = np.maximum(np.asarray(d_edt_sync_row, dtype=np.float64), 0.0) * dt_seconds / 1.0e9
    if not np.any(energy_losses > 0.0):
        return electron

    shifted = np.zeros_like(electron)
    synchrotron_lost = 0.0
    below_grid_lost = 0.0

    for idx, amount in enumerate(np.asarray(electron, dtype=np.float64)):
        if amount <= 0.0:
            continue
        initial_energy = float(energies[idx])
        target_energy = max(initial_energy - float(energy_losses[idx]), 0.0)
        particle_count = float(amount) * float(TRAPEZOID_WEIGHTS[idx])

        if target_energy <= 0.0:
            synchrotron_lost += particle_count * initial_energy
            continue

        synchrotron_lost += particle_count * max(initial_energy - target_energy, 0.0)

        if target_energy < energies[0]:
            below_grid_lost += particle_count * target_energy
            continue

        below_grid_lost += _deposit_to_log_grid(
            float(amount),
            float(TRAPEZOID_WEIGHTS[idx]),
            target_energy,
            shifted,
        )

    diagnostics["synchrotron_energy_lost"] += float(synchrotron_lost)
    diagnostics["below_grid_energy_lost"] += float(below_grid_lost)
    return shifted


def transport_step(
    gamma: np.ndarray,
    electron: np.ndarray,
    *,
    step_size: float,
    z_index: int,
    pp_log_extinction_row: np.ndarray,
    ics_log_extinction_row: np.ndarray,
    d_edt_ics_row: np.ndarray,
    pp_kernel: np.ndarray,
    pp_below_grid_row: np.ndarray,
    ics_gamma_kernel: np.ndarray,
    ics_electron_kernel: np.ndarray,
    ics_gamma_energy_row: np.ndarray,
    ics_below_grid_row: np.ndarray,
    cel_data: tuple[np.ndarray, np.ndarray] | None,
    diagnostics: dict[str, float],
    b_field_gauss: float = 0.0,
    b_field_gamma: float = 0.0,
    max_ics_tau_per_substep: float = DEFAULT_MAX_ICS_TAU_PER_SUBSTEP,
    max_ics_substeps: int = DEFAULT_MAX_ICS_SUBSTEPS,
) -> tuple[np.ndarray, np.ndarray]:
    """Advance the coupled gamma/electron cascade through one fixed-redshift block."""
    gamma = np.asarray(gamma, dtype=np.float64)
    electron = np.asarray(electron, dtype=np.float64)

    gamma_survived = np.exp(pp_log_extinction_row * float(step_size)) * gamma
    delta_gamma = gamma - gamma_survived
    diagnostics["pp_absorbed_energy"] += spectrum_energy(delta_gamma)
    diagnostics["below_grid_energy_lost"] += float(
        np.sum(delta_gamma * TRAPEZOID_WEIGHTS * pp_below_grid_row)
    )

    electron = electron + pp_kernel @ delta_gamma
    gamma = gamma_survived

    tau = np.maximum(-ics_log_extinction_row * float(step_size), 0.0)
    active_tau = tau[electron > 0.0]
    finite_tau = active_tau[np.isfinite(active_tau)]
    if active_tau.size and finite_tau.size != active_tau.size:
        n_substeps = int(max_ics_substeps)
    else:
        max_tau = float(np.max(finite_tau)) if finite_tau.size else 0.0
        n_substeps = max(1, int(np.ceil(max_tau / max(float(max_ics_tau_per_substep), 1.0e-12))))
        n_substeps = min(n_substeps, int(max_ics_substeps))
    sub_step = float(step_size) / n_substeps

    d_edt_sync_row = synchrotron_loss_rates(b_field_gauss, b_field_gamma, zReg[int(z_index)])
    cel_mask = None
    cel_targets = None
    if cel_data is not None:
        cel_mask, cel_targets = cel_data

    for _ in range(n_substeps):
        electron_survived = np.exp(ics_log_extinction_row * sub_step) * electron
        delta_electron = electron - electron_survived
        diagnostics["ics_scattered_electron_energy"] += spectrum_energy(delta_electron)

        if cel_mask is not None and np.any(cel_mask):
            cel_delta = np.where(cel_mask, delta_electron, 0.0)
            discrete_delta = delta_electron - cel_delta
        else:
            cel_delta = np.zeros_like(delta_electron)
            discrete_delta = delta_electron

        gamma = gamma + ics_gamma_kernel @ delta_electron
        electron_next = electron_survived + ics_electron_kernel @ discrete_delta

        if np.any(discrete_delta > 0.0):
            diagnostics["below_grid_energy_lost"] += float(
                np.sum(discrete_delta * TRAPEZOID_WEIGHTS * ics_below_grid_row)
            )

        if cel_mask is not None and np.any(cel_delta > 0.0):
            diagnostics["cel_scattered_electron_energy"] += spectrum_energy(cel_delta)
            for idx, amount in enumerate(cel_delta):
                if amount > 0.0:
                    diagnostics["below_grid_energy_lost"] += float(
                        amount
                        * TRAPEZOID_WEIGHTS[idx]
                        * max(energies[idx] - ics_gamma_energy_row[idx] - float(cel_targets[idx]), 0.0)
                    )
                    diagnostics["below_grid_energy_lost"] += _deposit_to_log_grid(
                        float(amount),
                        float(TRAPEZOID_WEIGHTS[idx]),
                        float(cel_targets[idx]),
                        electron_next,
                    )

        if b_field_gauss != 0.0:
            electron_next = _apply_synchrotron_losses(
                electron_next,
                d_edt_sync_row=d_edt_sync_row,
                sub_step=sub_step,
                diagnostics=diagnostics,
            )

        electron = np.maximum(electron_next, 0.0)
        gamma = np.maximum(gamma, 0.0)

    return gamma, electron


def _redshift_pair(
    gamma: np.ndarray,
    electron: np.ndarray,
    *,
    left_idx: np.ndarray,
    weights: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Redshift the gamma and electron spectra together across one window boundary."""
    return (
        redshift_cycle(gamma, left_idx=left_idx, weights=weights),
        redshift_cycle(electron, left_idx=left_idx, weights=weights),
    )


def _finalize_result(
    gamma: np.ndarray,
    electron: np.ndarray,
    diagnostics: dict[str, float],
    metadata: dict[str, Any],
) -> CascadeResult:
    """Package the final spectra and close the energy budget bookkeeping."""
    diagnostics = {key: float(value) for key, value in diagnostics.items()}
    diagnostics["final_gamma_energy"] = spectrum_energy(gamma)
    diagnostics["final_electron_energy"] = spectrum_energy(electron)
    initial = diagnostics["initial_gamma_energy"] + diagnostics["initial_electron_energy"]
    final_known = (
        diagnostics["final_gamma_energy"]
        + diagnostics["final_electron_energy"]
        + diagnostics["synchrotron_energy_lost"]
        + diagnostics["below_grid_energy_lost"]
        + diagnostics["redshift_energy_lost"]
    )
    diagnostics["energy_residual"] = float(initial - final_known)
    diagnostics["relative_energy_residual"] = (
        diagnostics["energy_residual"] / initial if initial > 0.0 else 0.0
    )
    return CascadeResult(
        gamma=np.asarray(gamma, dtype=np.float64),
        electron=np.asarray(electron, dtype=np.float64),
        diagnostics=diagnostics,
        metadata=metadata,
    )


def _run_transport_window(
    gamma: np.ndarray,
    electron: np.ndarray,
    *,
    window_idx: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    pp_log_extinction: np.ndarray,
    ics_log_extinction: np.ndarray,
    d_edt_ics: np.ndarray,
    get_pp_kernel: WeightedKernelGetter,
    get_pp_below_grid_energy: WeightedKernelGetter,
    get_ics_gamma_kernel: WeightedKernelGetter,
    get_ics_electron_kernel: WeightedKernelGetter,
    get_ics_gamma_row_energy: WeightedKernelGetter,
    get_ics_below_grid_energy: WeightedKernelGetter,
    get_cel_data: CelDataGetter | None,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    diagnostics: dict[str, float],
    b_field_gauss: float,
    b_field_gamma: float,
    cel_threshold_log_width: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Evolve the cascade through one redshift window and then redshift both spectra."""
    start = int(row_ptr[window_idx])
    stop = int(row_ptr[window_idx + 1])
    block_start = start
    while block_start < stop:
        z_index = int(zreg_indices[block_start])
        block_stop = block_start + 1
        while block_stop < stop and int(zreg_indices[block_stop]) == z_index:
            block_stop += 1
        block_step_size = float(np.sum(step_sizes[block_start:block_stop]))
        cel_data = None
        if get_cel_data is not None:
            cel_data = get_cel_data(z_index, cel_threshold_log_width)
        gamma, electron = transport_step(
            gamma,
            electron,
            step_size=block_step_size,
            z_index=z_index,
            pp_log_extinction_row=pp_log_extinction[z_index],
            ics_log_extinction_row=ics_log_extinction[z_index],
            d_edt_ics_row=d_edt_ics[z_index],
            pp_kernel=get_pp_kernel(z_index),
            pp_below_grid_row=get_pp_below_grid_energy(z_index),
            ics_gamma_kernel=get_ics_gamma_kernel(z_index),
            ics_electron_kernel=get_ics_electron_kernel(z_index),
            ics_gamma_energy_row=get_ics_gamma_row_energy(z_index),
            ics_below_grid_row=get_ics_below_grid_energy(z_index),
            cel_data=cel_data,
            diagnostics=diagnostics,
            b_field_gauss=b_field_gauss,
            b_field_gamma=b_field_gamma,
        )
        block_start = block_stop

    energy_before_redshift = spectrum_energy(gamma) + spectrum_energy(electron)
    gamma, electron = _redshift_pair(
        gamma,
        electron,
        left_idx=redshift_left_idx[window_idx],
        weights=redshift_weights[window_idx],
    )
    energy_after_redshift = spectrum_energy(gamma) + spectrum_energy(electron)
    diagnostics["redshift_energy_lost"] += max(energy_before_redshift - energy_after_redshift, 0.0)
    return gamma, electron


def propagate_point_transport(
    gamma_inj: np.ndarray,
    electron_inj: np.ndarray,
    *,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    pp_log_extinction: np.ndarray,
    ics_log_extinction: np.ndarray,
    d_edt_ics: np.ndarray,
    get_pp_kernel: WeightedKernelGetter,
    get_pp_below_grid_energy: WeightedKernelGetter,
    get_ics_gamma_kernel: WeightedKernelGetter,
    get_ics_electron_kernel: WeightedKernelGetter,
    get_ics_gamma_row_energy: WeightedKernelGetter,
    get_ics_below_grid_energy: WeightedKernelGetter,
    get_cel_data: CelDataGetter | None,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    b_field_gauss: float,
    b_field_gamma: float,
    cel_threshold_log_width: float = DEFAULT_CEL_LOG_BIN_THRESHOLD,
    progress=None,
) -> CascadeResult:
    """Propagate one source spectrum from its emission redshift down to z=0."""
    gamma = np.asarray(gamma_inj, dtype=np.float64).copy()
    electron = np.asarray(electron_inj, dtype=np.float64).copy()
    diagnostics = zero_diagnostics()
    diagnostics["initial_gamma_energy"] = spectrum_energy(gamma)
    diagnostics["initial_electron_energy"] = spectrum_energy(electron)

    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        gamma, electron = _run_transport_window(
            gamma,
            electron,
            window_idx=window_idx,
            step_sizes=step_sizes,
            row_ptr=row_ptr,
            zreg_indices=zreg_indices,
            pp_log_extinction=pp_log_extinction,
            ics_log_extinction=ics_log_extinction,
            d_edt_ics=d_edt_ics,
            get_pp_kernel=get_pp_kernel,
            get_pp_below_grid_energy=get_pp_below_grid_energy,
            get_ics_gamma_kernel=get_ics_gamma_kernel,
            get_ics_electron_kernel=get_ics_electron_kernel,
            get_ics_gamma_row_energy=get_ics_gamma_row_energy,
            get_ics_below_grid_energy=get_ics_below_grid_energy,
            get_cel_data=get_cel_data,
            redshift_left_idx=redshift_left_idx,
            redshift_weights=redshift_weights,
            diagnostics=diagnostics,
            b_field_gauss=b_field_gauss,
            b_field_gamma=b_field_gamma,
            cel_threshold_log_width=cel_threshold_log_width,
        )
        if progress is not None:
            progress.update(display_idx)

    return _finalize_result(gamma, electron, diagnostics, {"mode": "point", "z_max_index": int(z_max_index)})


def propagate_diffuse_transport(
    gamma_inj: np.ndarray,
    electron_inj: np.ndarray,
    *,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    pp_log_extinction: np.ndarray,
    ics_log_extinction: np.ndarray,
    d_edt_ics: np.ndarray,
    get_pp_kernel: WeightedKernelGetter,
    get_pp_below_grid_energy: WeightedKernelGetter,
    get_ics_gamma_kernel: WeightedKernelGetter,
    get_ics_electron_kernel: WeightedKernelGetter,
    get_ics_gamma_row_energy: WeightedKernelGetter,
    get_ics_below_grid_energy: WeightedKernelGetter,
    get_cel_data: CelDataGetter | None,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    b_field_gauss: float,
    b_field_gamma: float,
    cel_threshold_log_width: float = DEFAULT_CEL_LOG_BIN_THRESHOLD,
    progress=None,
) -> CascadeResult:
    """Propagate a diffuse source population with a single injected spectrum shape."""
    gamma = np.zeros(len(energies), dtype=np.float64)
    electron = np.zeros(len(energies), dtype=np.float64)
    diagnostics = zero_diagnostics()
    diagnostics["initial_gamma_energy"] = 0.0
    diagnostics["initial_electron_energy"] = 0.0

    gamma_template = np.asarray(gamma_inj, dtype=np.float64)
    electron_template = np.asarray(electron_inj, dtype=np.float64)

    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        gamma = gamma + volume_norms[window_idx] * gamma_template
        electron = electron + volume_norms[window_idx] * electron_template
        diagnostics["initial_gamma_energy"] += spectrum_energy(volume_norms[window_idx] * gamma_template)
        diagnostics["initial_electron_energy"] += spectrum_energy(volume_norms[window_idx] * electron_template)
        gamma, electron = _run_transport_window(
            gamma,
            electron,
            window_idx=window_idx,
            step_sizes=step_sizes,
            row_ptr=row_ptr,
            zreg_indices=zreg_indices,
            pp_log_extinction=pp_log_extinction,
            ics_log_extinction=ics_log_extinction,
            d_edt_ics=d_edt_ics,
            get_pp_kernel=get_pp_kernel,
            get_pp_below_grid_energy=get_pp_below_grid_energy,
            get_ics_gamma_kernel=get_ics_gamma_kernel,
            get_ics_electron_kernel=get_ics_electron_kernel,
            get_ics_gamma_row_energy=get_ics_gamma_row_energy,
            get_ics_below_grid_energy=get_ics_below_grid_energy,
            get_cel_data=get_cel_data,
            redshift_left_idx=redshift_left_idx,
            redshift_weights=redshift_weights,
            diagnostics=diagnostics,
            b_field_gauss=b_field_gauss,
            b_field_gamma=b_field_gamma,
            cel_threshold_log_width=cel_threshold_log_width,
        )
        if progress is not None:
            progress.update(display_idx)

    return _finalize_result(gamma, electron, diagnostics, {"mode": "diffuse", "z_max_index": int(z_max_index)})


def propagate_evolving_transport(
    gamma_inj: np.ndarray,
    electron_inj: np.ndarray,
    *,
    volume_norms: np.ndarray,
    z_max_index: int,
    step_sizes: np.ndarray,
    row_ptr: np.ndarray,
    zreg_indices: np.ndarray,
    pp_log_extinction: np.ndarray,
    ics_log_extinction: np.ndarray,
    d_edt_ics: np.ndarray,
    get_pp_kernel: WeightedKernelGetter,
    get_pp_below_grid_energy: WeightedKernelGetter,
    get_ics_gamma_kernel: WeightedKernelGetter,
    get_ics_electron_kernel: WeightedKernelGetter,
    get_ics_gamma_row_energy: WeightedKernelGetter,
    get_ics_below_grid_energy: WeightedKernelGetter,
    get_cel_data: CelDataGetter | None,
    redshift_left_idx: np.ndarray,
    redshift_weights: np.ndarray,
    b_field_gauss: float,
    b_field_gamma: float,
    cel_threshold_log_width: float = DEFAULT_CEL_LOG_BIN_THRESHOLD,
    progress=None,
) -> CascadeResult:
    """Propagate a diffuse population whose injection spectrum changes with redshift."""
    gamma = np.zeros(len(energies), dtype=np.float64)
    electron = np.zeros(len(energies), dtype=np.float64)
    diagnostics = zero_diagnostics()
    diagnostics["initial_gamma_energy"] = 0.0
    diagnostics["initial_electron_energy"] = 0.0

    gamma_grid = np.asarray(gamma_inj, dtype=np.float64)
    electron_grid = np.asarray(electron_inj, dtype=np.float64)

    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        gamma_source = volume_norms[window_idx] * gamma_grid[window_idx]
        electron_source = volume_norms[window_idx] * electron_grid[window_idx]
        gamma = gamma + gamma_source
        electron = electron + electron_source
        diagnostics["initial_gamma_energy"] += spectrum_energy(gamma_source)
        diagnostics["initial_electron_energy"] += spectrum_energy(electron_source)
        gamma, electron = _run_transport_window(
            gamma,
            electron,
            window_idx=window_idx,
            step_sizes=step_sizes,
            row_ptr=row_ptr,
            zreg_indices=zreg_indices,
            pp_log_extinction=pp_log_extinction,
            ics_log_extinction=ics_log_extinction,
            d_edt_ics=d_edt_ics,
            get_pp_kernel=get_pp_kernel,
            get_pp_below_grid_energy=get_pp_below_grid_energy,
            get_ics_gamma_kernel=get_ics_gamma_kernel,
            get_ics_electron_kernel=get_ics_electron_kernel,
            get_ics_gamma_row_energy=get_ics_gamma_row_energy,
            get_ics_below_grid_energy=get_ics_below_grid_energy,
            get_cel_data=get_cel_data,
            redshift_left_idx=redshift_left_idx,
            redshift_weights=redshift_weights,
            diagnostics=diagnostics,
            b_field_gauss=b_field_gauss,
            b_field_gamma=b_field_gamma,
            cel_threshold_log_width=cel_threshold_log_width,
        )
        if progress is not None:
            progress.update(display_idx)

    return _finalize_result(gamma, electron, diagnostics, {"mode": "evolving", "z_max_index": int(z_max_index)})
