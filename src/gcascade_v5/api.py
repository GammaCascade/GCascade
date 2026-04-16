from __future__ import annotations

import importlib
import os
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d

from . import bundle, builders, legacy
from .attenuation import (
    propagate_diffuse_attenuation_fast,
    propagate_diffuse_attenuation_numba,
    propagate_diffuse_attenuation_window_numba,
    propagate_evolving_attenuation_fast,
    propagate_evolving_attenuation_numba,
    propagate_evolving_attenuation_window_numba,
    propagate_point_attenuation_fast,
    propagate_point_attenuation_numba,
    propagate_point_attenuation_window_numba,
)
from .cascade import (
    propagate_diffuse_cascade_fast,
    propagate_diffuse_cascade_numba,
    propagate_diffuse_cascade_window_numba,
    propagate_evolving_cascade_fast,
    propagate_evolving_cascade_numba,
    propagate_evolving_cascade_window_numba,
    propagate_point_cascade_fast,
    propagate_point_cascade_numba,
    propagate_point_cascade_window_numba,
)
from .config import (
    ProgressBar,
    get_numba_enabled,
    is_numba_available,
    set_numba,
    set_progress,
    status,
)
from .redshift import diffuse_distances_index, redshift_cycle
from .state import (
    get_generated_library_path,
    get_library_path,
    reset_state,
    set_generated_library_path,
    set_library_path,
    state,
)


Mpc = legacy.Mpc
c = legacy.c
t = legacy.t
mu0 = legacy.mu0
elecmass = legacy.elecmass
echarge = legacy.echarge
sigmaTe = legacy.sigmaTe

energies = legacy.energies
dEnergiesGamma = legacy.dEnergiesGamma
diffuseSteps = legacy.diffuseSteps
diffuseDistances = legacy.diffuseDistances
zReg = legacy.zReg
EBL_NAME_MAP = legacy.EBL_NAME_MAP
EBL_DESCRIPTION_MAP = legacy.EBL_DESCRIPTION_MAP
LIB_PATH_ENV_VAR = legacy.LIB_PATH_ENV_VAR
GENERATED_LIB_PATH_ENV_VAR = legacy.GENERATED_LIB_PATH_ENV_VAR

hubble = legacy.hubble
cutoffPowerLaw = legacy.cutoffPowerLaw
specPlot = legacy.specPlot

_FALLBACK_RED_LEFT_IDX, _FALLBACK_RED_WEIGHTS, _ = bundle._build_redshift_tables()
_state_module = importlib.import_module(".state", __package__)


def _validate_z_start(z_start: float) -> float:
    z = float(z_start)
    if z < 0.0:
        raise ValueError("zStart must be non-negative")
    if z > 10.0:
        raise ValueError("GCascadeV5 does not allow sources with redshift above z=10")
    return z


def _validate_spectrum_1d(inj: np.ndarray | list[float], announce: bool = True) -> np.ndarray:
    arr = np.asarray(inj, dtype=np.float64)
    if arr.shape != energies.shape:
        raise ValueError(f"Injected spectrum must have shape {energies.shape}, got {arr.shape}")
    if announce:
        status("Injected spectrum properly formatted.")
    return arr


def _validate_spectrum_2d(inj: np.ndarray | list[list[float]], announce: bool = True) -> np.ndarray:
    arr = np.asarray(inj, dtype=np.float64)
    target = (len(diffuseDistances), len(energies))
    if arr.shape != target:
        raise ValueError(f"Injected evolving spectrum must have shape {target}, got {arr.shape}")
    if announce:
        status("Injected spectrum properly formatted.")
    return arr


def _validate_z_distribution(z_distrib: np.ndarray | list[float], announce: bool = True) -> np.ndarray:
    arr = np.asarray(z_distrib, dtype=np.float64)
    if arr.shape != diffuseDistances.shape:
        raise ValueError(
            f"Comoving density distribution must have shape {diffuseDistances.shape}, got {arr.shape}"
        )
    if announce:
        status("Comoving density distribution properly formatted.")
    return arr


def _redshift_cycle_callback(spec: np.ndarray, window_idx: int) -> np.ndarray:
    runtime = state()
    left_idx = _FALLBACK_RED_LEFT_IDX[window_idx]
    weights = _FALLBACK_RED_WEIGHTS[window_idx]
    if runtime.redshift_left_idx is not None and runtime.redshift_weights is not None:
        left_idx = runtime.redshift_left_idx[window_idx]
        weights = runtime.redshift_weights[window_idx]
    return redshift_cycle(spec, left_idx=left_idx, weights=weights, use_numba=get_numba_enabled())


def _volume_norms(z_start: float, z_distrib: np.ndarray) -> np.ndarray:
    z_max_index = diffuse_distances_index(z_start) + 1
    z_interp = interp1d(
        diffuseDistances,
        z_distrib,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True,
    )
    z_grid = diffuseDistances[:z_max_index]
    return (c * Mpc * z_interp(z_grid) * diffuseSteps[:z_max_index]) / hubble(z_grid)


def diffuseDistancesIndex(x: float) -> int:
    return diffuse_distances_index(x)


def RedshiftPoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)
    status(f"Running RedshiftPoint for zStart={z_start:.6g}")

    z_max_index = diffuse_distances_index(z_start) + 1
    final_result = inj_spectra.copy()
    progress = ProgressBar("RedshiftPoint", z_max_index)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = _redshift_cycle_callback(final_result, window_idx)
        progress.update(display_idx)

    d_l = legacy._luminosity_distance_mpc(z_start)
    status("RedshiftPoint completed.")
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def AttenuatePoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)
    status(f"Running AttenuatePoint for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    z_max_index = diffuse_distances_index(z_start) + 1
    progress = ProgressBar("AttenuatePoint", z_max_index)

    if get_numba_enabled() and runtime.attenuation_vectors is not None and runtime.redshift_left_idx is not None:
        final_result = inj_spectra.copy()
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_point_attenuation_window_numba(
                final_result,
                runtime.attenuation_vectors[window_idx],
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_point_attenuation_fast(
            inj_spectra,
            z_max_index,
            runtime.attenuation_vectors,
            _redshift_cycle_callback,
            progress=progress,
        )

    d_l = legacy._luminosity_distance_mpc(z_start)
    status("AttenuatePoint completed.")
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def CascadePoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)
    status(f"Running CascadePoint for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    z_max_index = diffuse_distances_index(z_start) + 1
    progress = ProgressBar("CascadePoint", z_max_index)

    if get_numba_enabled() and runtime.cycle_packed_array is not None:
        final_result = inj_spectra.copy()
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_point_cascade_window_numba(
                final_result,
                int(runtime.row_ptr[window_idx]),
                int(runtime.row_ptr[window_idx + 1]),
                runtime.step_sizes,
                runtime.zreg_indices,
                runtime.extinction_coeffs,
                runtime.cycle_packed_array,
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_point_cascade_fast(
            inj_spectra,
            z_max_index,
            runtime.step_sizes,
            runtime.row_ptr,
            runtime.zreg_indices,
            runtime.log_extinction_coeffs,
            runtime.get_weighted_cycle_slice,
            runtime.redshift_left_idx,
            runtime.redshift_weights,
            progress=progress,
        )

    d_l = legacy._luminosity_distance_mpc(z_start)
    status("CascadePoint completed.")
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def RedshiftDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running RedshiftDiffuse for zStart={z_start:.6g}")

    volume_norms = _volume_norms(z_start, z_distrib)
    final_result = np.zeros_like(inj_spectra)
    progress = ProgressBar("RedshiftDiffuse", len(volume_norms))
    for display_idx, window_idx in enumerate(range(len(volume_norms) - 1, -1, -1), start=1):
        final_result = _redshift_cycle_callback(final_result + volume_norms[window_idx] * inj_spectra, window_idx)
        progress.update(display_idx)

    status("RedshiftDiffuse completed.")
    return final_result / (4.0 * np.pi)


def AttenuateDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running AttenuateDiffuse for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("AttenuateDiffuse", z_max_index)

    if get_numba_enabled() and runtime.attenuation_vectors is not None:
        final_result = np.zeros(len(energies), dtype=np.float64)
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_diffuse_attenuation_window_numba(
                final_result,
                volume_norms[window_idx] * inj_spectra,
                runtime.attenuation_vectors[window_idx],
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_diffuse_attenuation_fast(
            inj_spectra,
            volume_norms,
            z_max_index,
            runtime.attenuation_vectors,
            _redshift_cycle_callback,
            progress=progress,
        )

    status("AttenuateDiffuse completed.")
    return final_result / (4.0 * np.pi)


def CascadeDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running CascadeDiffuse for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("CascadeDiffuse", z_max_index)

    if get_numba_enabled() and runtime.cycle_packed_array is not None:
        final_result = np.zeros(len(energies), dtype=np.float64)
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_diffuse_cascade_window_numba(
                final_result,
                volume_norms[window_idx] * inj_spectra,
                int(runtime.row_ptr[window_idx]),
                int(runtime.row_ptr[window_idx + 1]),
                runtime.step_sizes,
                runtime.zreg_indices,
                runtime.extinction_coeffs,
                runtime.cycle_packed_array,
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_diffuse_cascade_fast(
            inj_spectra,
            volume_norms,
            z_max_index,
            runtime.step_sizes,
            runtime.row_ptr,
            runtime.zreg_indices,
            runtime.log_extinction_coeffs,
            runtime.get_weighted_cycle_slice,
            runtime.redshift_left_idx,
            runtime.redshift_weights,
            progress=progress,
        )

    status("CascadeDiffuse completed.")
    return final_result / (4.0 * np.pi)


def RedshiftEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running RedshiftEvolving for zStart={z_start:.6g}")

    volume_norms = _volume_norms(z_start, z_distrib)
    final_result = np.zeros(len(energies), dtype=np.float64)
    progress = ProgressBar("RedshiftEvolving", len(volume_norms))
    for display_idx, window_idx in enumerate(range(len(volume_norms) - 1, -1, -1), start=1):
        final_result = _redshift_cycle_callback(final_result + volume_norms[window_idx] * inj_spectra[window_idx], window_idx)
        progress.update(display_idx)

    status("RedshiftEvolving completed.")
    return final_result / (4.0 * np.pi)


def AttenuateEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running AttenuateEvolving for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("AttenuateEvolving", z_max_index)

    if get_numba_enabled() and runtime.attenuation_vectors is not None:
        final_result = np.zeros(len(energies), dtype=np.float64)
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_evolving_attenuation_window_numba(
                final_result,
                volume_norms[window_idx] * inj_spectra[window_idx],
                runtime.attenuation_vectors[window_idx],
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_evolving_attenuation_fast(
            inj_spectra,
            volume_norms,
            z_max_index,
            runtime.attenuation_vectors,
            _redshift_cycle_callback,
            progress=progress,
        )

    status("AttenuateEvolving completed.")
    return final_result / (4.0 * np.pi)


def CascadeEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running CascadeEvolving for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("CascadeEvolving", z_max_index)

    if get_numba_enabled() and runtime.cycle_packed_array is not None:
        final_result = np.zeros(len(energies), dtype=np.float64)
        for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
            final_result = propagate_evolving_cascade_window_numba(
                final_result,
                volume_norms[window_idx] * inj_spectra[window_idx],
                int(runtime.row_ptr[window_idx]),
                int(runtime.row_ptr[window_idx + 1]),
                runtime.step_sizes,
                runtime.zreg_indices,
                runtime.extinction_coeffs,
                runtime.cycle_packed_array,
                runtime.redshift_left_idx[window_idx],
                runtime.redshift_weights[window_idx],
            )
            progress.update(display_idx)
    else:
        final_result = propagate_evolving_cascade_fast(
            inj_spectra,
            volume_norms,
            z_max_index,
            runtime.step_sizes,
            runtime.row_ptr,
            runtime.zreg_indices,
            runtime.log_extinction_coeffs,
            runtime.get_weighted_cycle_slice,
            runtime.redshift_left_idx,
            runtime.redshift_weights,
            progress=progress,
        )

    status("CascadeEvolving completed.")
    return final_result / (4.0 * np.pi)


def changeEBLModel(EBL: int) -> None:
    runtime = state()
    new_ebl = int(EBL)
    if new_ebl not in EBL_NAME_MAP:
        raise ValueError(f"Invalid EBL index: {new_ebl}")
    if new_ebl == runtime.ebl_index:
        raise ValueError(f"{EBL_DESCRIPTION_MAP[new_ebl]} is already the current EBL model")

    old_ebl = runtime.ebl_index
    runtime.set_ebl_index(new_ebl)
    runtime.ensure_ebl_loaded()
    _state_module.EBLindex = new_ebl
    status(
        "EBL model changed successfully from "
        f"{EBL_DESCRIPTION_MAP[old_ebl]} to {EBL_DESCRIPTION_MAP[new_ebl]}."
    )


def changeMagneticField(BField: float, gamma: float, EBL: int) -> None:
    runtime = state()
    runtime.ensure_bundle_loaded()

    ebl = int(EBL)
    if ebl not in EBL_NAME_MAP:
        raise ValueError(f"Invalid EBL index: {ebl}")

    status(
        "Changing magnetic field cycle tables for EBL index "
        f"{ebl} with B(z)={float(BField):.6g}*(1+z)^{float(gamma):.6g} Gauss."
    )
    target = builders.generate_magnetic_field_variant(
        runtime.library_path,
        runtime.generated_library_path,
        ebl_index=ebl,
        b_field=float(BField),
        gamma=float(gamma),
    )

    runtime.set_ebl_index(ebl)
    runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    _state_module.EBLindex = ebl
    status(f"Magnetic field update completed. {target.name} was written, activated, and EBL model {EBL_DESCRIPTION_MAP[ebl]} is now active.")


def convert_legacy_library(
    source_path: str | os.PathLike[str],
    target_path: str | os.PathLike[str],
    *,
    overwrite: bool = False,
) -> Path:
    return bundle.convert_legacy_library(source_path, target_path, overwrite=overwrite)


def get_bundle_info() -> dict[str, object]:
    return bundle.bundle_info(get_library_path())


def list_generated_variants(ebl_index: int) -> list[dict[str, object]]:
    return bundle.generated_variants(get_library_path(), int(ebl_index))


def get_active_cycle_path(ebl_index: int | None = None) -> Path:
    runtime = state()
    target_ebl = runtime.ebl_index if ebl_index is None else int(ebl_index)
    return bundle.resolve_active_cycle_path(get_library_path(), target_ebl)


def set_active_cycle_path(
    path: str | os.PathLike[str],
    ebl_index: int | None = None,
    *,
    switch_ebl: bool = True,
) -> Path:
    target = bundle.activate_cycle_path(get_library_path(), path, ebl_index=ebl_index)
    target_ebl = int(ebl_index) if ebl_index is not None else bundle.infer_cycle_ebl_index(target)
    runtime = state()
    if switch_ebl:
        runtime.set_ebl_index(target_ebl)
        _state_module.EBLindex = target_ebl
    if runtime.ebl_index == target_ebl:
        runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    return target


def reset_factory_settings() -> None:
    runtime = state()
    runtime.ensure_bundle_loaded()
    bundle.clear_active_generated_variant(runtime.library_path)
    runtime.set_ebl_index(1)
    runtime.ensure_ebl_loaded(refresh_active_cycle=True)
    _state_module.EBLindex = 1
    status("Factory settings restored: EBL model set to Saldana-Lopez et al. (2021) and cycle tables reset to defaults.")


def export_active_cycle_to_legacy_mat(ebl_index: int, target_path: str | os.PathLike[str]) -> Path:
    return bundle.export_active_cycle_to_legacy_mat(get_library_path(), int(ebl_index), target_path)


redshift_point = RedshiftPoint
attenuate_point = AttenuatePoint
cascade_point = CascadePoint
redshift_diffuse = RedshiftDiffuse
attenuate_diffuse = AttenuateDiffuse
cascade_diffuse = CascadeDiffuse
redshift_evolving = RedshiftEvolving
attenuate_evolving = AttenuateEvolving
cascade_evolving = CascadeEvolving
change_ebl_model = changeEBLModel
change_magnetic_field = changeMagneticField


__all__ = [
    "Mpc",
    "c",
    "t",
    "mu0",
    "elecmass",
    "echarge",
    "sigmaTe",
    "energies",
    "dEnergiesGamma",
    "diffuseSteps",
    "diffuseDistances",
    "zReg",
    "EBL_NAME_MAP",
    "EBL_DESCRIPTION_MAP",
    "LIB_PATH_ENV_VAR",
    "GENERATED_LIB_PATH_ENV_VAR",
    "hubble",
    "cutoffPowerLaw",
    "specPlot",
    "set_progress",
    "set_numba",
    "is_numba_available",
    "get_numba_enabled",
    "reset_state",
    "set_library_path",
    "set_generated_library_path",
    "get_library_path",
    "get_generated_library_path",
    "diffuseDistancesIndex",
    "RedshiftPoint",
    "AttenuatePoint",
    "CascadePoint",
    "RedshiftDiffuse",
    "AttenuateDiffuse",
    "CascadeDiffuse",
    "RedshiftEvolving",
    "AttenuateEvolving",
    "CascadeEvolving",
    "changeEBLModel",
    "changeMagneticField",
    "convert_legacy_library",
    "get_bundle_info",
    "list_generated_variants",
    "get_active_cycle_path",
    "set_active_cycle_path",
    "reset_factory_settings",
    "export_active_cycle_to_legacy_mat",
    "redshift_point",
    "attenuate_point",
    "cascade_point",
    "redshift_diffuse",
    "attenuate_diffuse",
    "cascade_diffuse",
    "redshift_evolving",
    "attenuate_evolving",
    "cascade_evolving",
    "change_ebl_model",
    "change_magnetic_field",
]
