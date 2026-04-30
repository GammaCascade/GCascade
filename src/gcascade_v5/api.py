from __future__ import annotations

"""Public GCascadeV5 interface for photon and electron cascade calculations."""

import importlib

import numpy as np
from scipy.interpolate import interp1d

from . import bundle
from .attenuation import (
    propagate_diffuse_attenuation_fast,
    propagate_evolving_attenuation_fast,
    propagate_point_attenuation_fast,
)
from .config import ProgressBar, set_progress, status
from .physics import (
    EBL_DESCRIPTION_MAP,
    EBL_NAME_MAP,
    Mpc,
    c,
    cutoffPowerLaw,
    dEnergiesGamma,
    diffuseDistances,
    diffuseSteps,
    echarge,
    elecmass,
    energies,
    hubble,
    luminosity_distance_mpc,
    mu0,
    sigmaTe,
    specPlot,
    t,
    zReg,
)
from .redshift import diffuse_distances_index, redshift_cycle
from .results import CascadeResult
from .state import LIB_PATH_ENV_VAR, get_library_path, reset_state, set_library_path, state
from .transport import (
    DEFAULT_CEL_LOG_BIN_THRESHOLD,
    propagate_diffuse_transport,
    propagate_evolving_transport,
    propagate_point_transport,
)


_FALLBACK_RED_LEFT_IDX, _FALLBACK_RED_WEIGHTS, _ = bundle._build_redshift_tables()
_state_module = importlib.import_module(".state", __package__)


def _validate_z_start(z_start: float) -> float:
    """Check that the source redshift lies inside the tabulated propagation range."""
    z = float(z_start)
    if z < 0.0:
        raise ValueError("zStart must be non-negative")
    if z > 10.0:
        raise ValueError("GCascadeV5 does not allow sources with redshift above z=10")
    return z


def _validate_spectrum_1d(inj: np.ndarray | list[float], announce: bool = True) -> np.ndarray:
    """Check that a one-dimensional injected spectrum matches the main energy grid."""
    arr = np.asarray(inj, dtype=np.float64)
    if arr.shape != energies.shape:
        raise ValueError(f"Injected spectrum must have shape {energies.shape}, got {arr.shape}")
    if announce:
        status("Injected spectrum properly formatted.")
    return arr


def _validate_optional_spectrum_1d(
    inj: np.ndarray | list[float] | None,
    *,
    label: str,
    announce: bool = True,
) -> np.ndarray:
    """Validate an optional one-dimensional spectrum, defaulting to zero injection."""
    if inj is None:
        return np.zeros_like(energies, dtype=np.float64)
    arr = np.asarray(inj, dtype=np.float64)
    if arr.shape != energies.shape:
        raise ValueError(f"{label} spectrum must have shape {energies.shape}, got {arr.shape}")
    if announce:
        status(f"{label} spectrum properly formatted.")
    return arr


def _validate_spectrum_2d(inj: np.ndarray | list[list[float]], announce: bool = True) -> np.ndarray:
    """Check that an evolving injection history covers all redshift windows and energies."""
    arr = np.asarray(inj, dtype=np.float64)
    target = (len(diffuseDistances), len(energies))
    if arr.shape != target:
        raise ValueError(f"Injected evolving spectrum must have shape {target}, got {arr.shape}")
    if announce:
        status("Injected spectrum properly formatted.")
    return arr


def _validate_optional_spectrum_2d(
    inj: np.ndarray | list[list[float]] | None,
    *,
    label: str,
    announce: bool = True,
) -> np.ndarray:
    """Validate an optional evolving spectrum, defaulting to no electron injection."""
    target = (len(diffuseDistances), len(energies))
    if inj is None:
        return np.zeros(target, dtype=np.float64)
    arr = np.asarray(inj, dtype=np.float64)
    if arr.shape != target:
        raise ValueError(f"{label} evolving spectrum must have shape {target}, got {arr.shape}")
    if announce:
        status(f"{label} spectrum properly formatted.")
    return arr


def _validate_z_distribution(z_distrib: np.ndarray | list[float], announce: bool = True) -> np.ndarray:
    """Check that the comoving source-density history is tabulated on the bundle redshift grid."""
    arr = np.asarray(z_distrib, dtype=np.float64)
    if arr.shape != diffuseDistances.shape:
        raise ValueError(
            f"Comoving density distribution must have shape {diffuseDistances.shape}, got {arr.shape}"
        )
    if announce:
        status("Comoving density distribution properly formatted.")
    return arr


def _redshift_cycle_callback(spec: np.ndarray, window_idx: int) -> np.ndarray:
    """Apply one precomputed redshift step using the active bundle tables."""
    runtime = state()
    left_idx = _FALLBACK_RED_LEFT_IDX[window_idx]
    weights = _FALLBACK_RED_WEIGHTS[window_idx]
    if runtime.redshift_left_idx is not None and runtime.redshift_weights is not None:
        left_idx = runtime.redshift_left_idx[window_idx]
        weights = runtime.redshift_weights[window_idx]
    return redshift_cycle(spec, left_idx=left_idx, weights=weights)


def _volume_norms(z_start: float, z_distrib: np.ndarray) -> np.ndarray:
    """Convert a comoving emissivity history into per-window source power weights."""
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
    """Return the closest tabulated redshift-window index for a given z value."""
    return diffuse_distances_index(x)


def RedshiftPoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    """Redshift a point-source spectrum to the observer without interactions."""
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)
    status(f"Running RedshiftPoint for zStart={z_start:.6g}")

    z_max_index = diffuse_distances_index(z_start) + 1
    final_result = inj_spectra.copy()
    progress = ProgressBar("RedshiftPoint", z_max_index)
    for display_idx, window_idx in enumerate(range(z_max_index - 1, -1, -1), start=1):
        final_result = _redshift_cycle_callback(final_result, window_idx)
        progress.update(display_idx)

    d_l = luminosity_distance_mpc(z_start)
    status("RedshiftPoint completed.")
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def AttenuatePoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    """Propagate a point-source photon spectrum with absorption but no secondaries."""
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)
    status(f"Running AttenuatePoint for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    z_max_index = diffuse_distances_index(z_start) + 1
    progress = ProgressBar("AttenuatePoint", z_max_index)
    final_result = propagate_point_attenuation_fast(
        inj_spectra,
        z_max_index,
        runtime.attenuation_vectors,
        _redshift_cycle_callback,
        progress=progress,
    )

    d_l = luminosity_distance_mpc(z_start)
    status("AttenuatePoint completed.")
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def CascadePoint(
    gammaSpectraPre: np.ndarray | list[float],
    zStart: float,
    electronSpectraPre: np.ndarray | list[float] | None = None,
    *,
    return_state: bool = False,
) -> np.ndarray | CascadeResult:
    """Run the full gamma/electron cascade for a single source redshift."""
    inj_spectra = _validate_spectrum_1d(gammaSpectraPre)
    electron_spectra = _validate_optional_spectrum_1d(
        electronSpectraPre,
        label="Injected electron",
        announce=electronSpectraPre is not None,
    )
    z_start = _validate_z_start(zStart)
    status(f"Running CascadePoint for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_transport_loaded()
    z_max_index = diffuse_distances_index(z_start) + 1
    progress = ProgressBar("CascadePoint", z_max_index)

    result = propagate_point_transport(
        inj_spectra,
        electron_spectra,
        z_max_index=z_max_index,
        step_sizes=runtime.step_sizes,
        row_ptr=runtime.row_ptr,
        zreg_indices=runtime.zreg_indices,
        pp_log_extinction=runtime.log_extinction_coeffs,
        ics_log_extinction=runtime.log_ics_extinction_coeffs,
        d_edt_ics=runtime.d_edt_ics,
        get_pp_kernel=runtime.get_weighted_pp_slice,
        get_pp_below_grid_energy=runtime.get_pp_below_grid_energy,
        get_ics_gamma_kernel=runtime.get_weighted_ics_gamma_slice,
        get_ics_electron_kernel=runtime.get_weighted_ics_electron_slice,
        get_ics_gamma_row_energy=runtime.get_ics_gamma_row_energy,
        get_ics_below_grid_energy=runtime.get_ics_below_grid_energy,
        get_cel_data=runtime.get_ics_cel_data,
        redshift_left_idx=runtime.redshift_left_idx,
        redshift_weights=runtime.redshift_weights,
        b_field_gauss=runtime.b_field_gauss,
        b_field_gamma=runtime.b_field_gamma,
        cel_threshold_log_width=DEFAULT_CEL_LOG_BIN_THRESHOLD,
        progress=progress,
    )

    d_l = luminosity_distance_mpc(z_start)
    scaled = result.with_scaled_spectra(((1.0 + z_start) ** 2) / (4.0 * np.pi * (d_l * Mpc) ** 2))
    status("CascadePoint completed.")
    return scaled if return_state else scaled.gamma


def RedshiftDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    """Redshift a diffuse emissivity history without any interactions."""
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
    """Propagate a diffuse photon emissivity with absorption but no secondaries."""
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running AttenuateDiffuse for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("AttenuateDiffuse", z_max_index)
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
    gammaSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
    electronSpectra: np.ndarray | list[float] | None = None,
    *,
    return_state: bool = False,
) -> np.ndarray | CascadeResult:
    """Run the full cascade for a diffuse source class with one spectral template."""
    inj_spectra = _validate_spectrum_1d(gammaSpectra)
    electron_spectra = _validate_optional_spectrum_1d(
        electronSpectra,
        label="Injected electron",
        announce=electronSpectra is not None,
    )
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running CascadeDiffuse for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_transport_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("CascadeDiffuse", z_max_index)

    result = propagate_diffuse_transport(
        inj_spectra,
        electron_spectra,
        volume_norms=volume_norms,
        z_max_index=z_max_index,
        step_sizes=runtime.step_sizes,
        row_ptr=runtime.row_ptr,
        zreg_indices=runtime.zreg_indices,
        pp_log_extinction=runtime.log_extinction_coeffs,
        ics_log_extinction=runtime.log_ics_extinction_coeffs,
        d_edt_ics=runtime.d_edt_ics,
        get_pp_kernel=runtime.get_weighted_pp_slice,
        get_pp_below_grid_energy=runtime.get_pp_below_grid_energy,
        get_ics_gamma_kernel=runtime.get_weighted_ics_gamma_slice,
        get_ics_electron_kernel=runtime.get_weighted_ics_electron_slice,
        get_ics_gamma_row_energy=runtime.get_ics_gamma_row_energy,
        get_ics_below_grid_energy=runtime.get_ics_below_grid_energy,
        get_cel_data=runtime.get_ics_cel_data,
        redshift_left_idx=runtime.redshift_left_idx,
        redshift_weights=runtime.redshift_weights,
        b_field_gauss=runtime.b_field_gauss,
        b_field_gamma=runtime.b_field_gamma,
        cel_threshold_log_width=DEFAULT_CEL_LOG_BIN_THRESHOLD,
        progress=progress,
    ).with_scaled_spectra(1.0 / (4.0 * np.pi))

    status("CascadeDiffuse completed.")
    return result if return_state else result.gamma


def RedshiftEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    """Redshift a source history whose injected spectrum itself evolves with redshift."""
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
    """Propagate an evolving photon emissivity with absorption but no secondary particles."""
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running AttenuateEvolving for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_ebl_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("AttenuateEvolving", z_max_index)
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
    gammaSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
    electronSpectra: np.ndarray | list[list[float]] | None = None,
    *,
    return_state: bool = False,
) -> np.ndarray | CascadeResult:
    """Run the full cascade for a source history with redshift-dependent injection."""
    inj_spectra = _validate_spectrum_2d(gammaSpectra)
    electron_spectra = _validate_optional_spectrum_2d(
        electronSpectra,
        label="Injected electron",
        announce=electronSpectra is not None,
    )
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)
    status(f"Running CascadeEvolving for zStart={z_start:.6g}")

    runtime = state()
    runtime.ensure_transport_loaded()
    volume_norms = _volume_norms(z_start, z_distrib)
    z_max_index = len(volume_norms)
    progress = ProgressBar("CascadeEvolving", z_max_index)

    result = propagate_evolving_transport(
        inj_spectra,
        electron_spectra,
        volume_norms=volume_norms,
        z_max_index=z_max_index,
        step_sizes=runtime.step_sizes,
        row_ptr=runtime.row_ptr,
        zreg_indices=runtime.zreg_indices,
        pp_log_extinction=runtime.log_extinction_coeffs,
        ics_log_extinction=runtime.log_ics_extinction_coeffs,
        d_edt_ics=runtime.d_edt_ics,
        get_pp_kernel=runtime.get_weighted_pp_slice,
        get_pp_below_grid_energy=runtime.get_pp_below_grid_energy,
        get_ics_gamma_kernel=runtime.get_weighted_ics_gamma_slice,
        get_ics_electron_kernel=runtime.get_weighted_ics_electron_slice,
        get_ics_gamma_row_energy=runtime.get_ics_gamma_row_energy,
        get_ics_below_grid_energy=runtime.get_ics_below_grid_energy,
        get_cel_data=runtime.get_ics_cel_data,
        redshift_left_idx=runtime.redshift_left_idx,
        redshift_weights=runtime.redshift_weights,
        b_field_gauss=runtime.b_field_gauss,
        b_field_gamma=runtime.b_field_gamma,
        cel_threshold_log_width=DEFAULT_CEL_LOG_BIN_THRESHOLD,
        progress=progress,
    ).with_scaled_spectra(1.0 / (4.0 * np.pi))

    status("CascadeEvolving completed.")
    return result if return_state else result.gamma


def changeEBLModel(EBL: int) -> None:
    """Switch the cascade to a different EBL background model."""
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
    """Set the magnetic-field law that competes with ICS through synchrotron cooling."""
    runtime = state()
    runtime.ensure_bundle_loaded()

    ebl = int(EBL)
    if ebl not in EBL_NAME_MAP:
        raise ValueError(f"Invalid EBL index: {ebl}")

    status(
        "Setting electron-transport magnetic field for EBL index "
        f"{ebl} with B(z)={float(BField):.6g}*(1+z)^{float(gamma):.6g} Gauss."
    )
    runtime.set_ebl_index(ebl)
    runtime.set_magnetic_field(float(BField), float(gamma))
    runtime.ensure_transport_loaded()
    _state_module.EBLindex = ebl
    status(
        "Magnetic field update completed. "
        "Electron transport will apply synchrotron-loss competition with "
        f"EBL model {EBL_DESCRIPTION_MAP[ebl]} active."
    )


def get_bundle_info() -> dict[str, object]:
    """Return a summary of the active runtime bundle."""
    return bundle.bundle_info(get_library_path())


def reset_factory_settings() -> None:
    """Restore the default EBL model and turn off synchrotron competition."""
    runtime = state()
    runtime.ensure_bundle_loaded()
    runtime.set_ebl_index(1)
    runtime.set_magnetic_field(0.0, 0.0)
    runtime.ensure_ebl_loaded()
    _state_module.EBLindex = 1
    status("Factory settings restored: EBL model set to Saldana-Lopez et al. (2021) and magnetic field reset.")


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
    "CascadeResult",
    "energies",
    "dEnergiesGamma",
    "diffuseSteps",
    "diffuseDistances",
    "zReg",
    "EBL_NAME_MAP",
    "EBL_DESCRIPTION_MAP",
    "LIB_PATH_ENV_VAR",
    "hubble",
    "cutoffPowerLaw",
    "specPlot",
    "set_progress",
    "reset_state",
    "set_library_path",
    "get_library_path",
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
    "get_bundle_info",
    "reset_factory_settings",
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
