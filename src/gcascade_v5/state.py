from __future__ import annotations

"""Runtime state for the active GCascadeV5 library bundle and EBL selection."""

from collections import OrderedDict
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from . import bundle
from .physics import EBL_NAME_MAP, dEnergiesGamma, energies


LIB_PATH_ENV_VAR = "GCASCADE_LIB_PATH"
_TRAPEZOID_WEIGHTS = np.empty(len(energies), dtype=np.float64)
_TRAPEZOID_WEIGHTS[0] = dEnergiesGamma[0]
_TRAPEZOID_WEIGHTS[-1] = dEnergiesGamma[-1]
_TRAPEZOID_WEIGHTS[1:-1] = dEnergiesGamma[:-1] + dEnergiesGamma[1:]


def _library_path_config_message() -> str:
    """Explain how to point GCascade at a different runtime bundle."""
    return (
        f"Configure with {LIB_PATH_ENV_VAR} (before import) or call "
        "gcascade_v5.set_library_path('/path/to/LibrariesV5')."
    )


def _discover_default_library_path() -> Path:
    """Search standard locations for the default GCascadeV5 runtime bundle."""
    env_path = os.getenv(LIB_PATH_ENV_VAR)
    if env_path:
        return Path(env_path).expanduser().resolve()

    candidates = [
        Path.cwd() / "LibrariesV5",
        Path.cwd() / "GCascade" / "LibrariesV5",
        Path.home() / "GCascade" / "LibrariesV5",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate.resolve()
    return (Path.cwd() / "LibrariesV5").resolve()


@dataclass
class RuntimeBundleState:
    """In-memory handles and cached tables for the active cascade bundle."""

    library_path: Path
    ebl_index: int = 1

    manifest: dict[str, Any] | None = None
    step_sizes: np.ndarray | None = None
    row_ptr: np.ndarray | None = None
    zreg_indices: np.ndarray | None = None
    redshift_left_idx: np.ndarray | None = None
    redshift_weights: np.ndarray | None = None
    redshift_scales: np.ndarray | None = None
    packed_row_ptr: np.ndarray | None = None

    runtime_ebl_file: h5py.File | None = None
    imfp: np.ndarray | None = None
    extinction_coeffs: np.ndarray | None = None
    log_extinction_coeffs: np.ndarray | None = None
    attenuation_vectors: np.ndarray | None = None
    ics_imfp: np.ndarray | None = None
    ics_extinction_coeffs: np.ndarray | None = None
    log_ics_extinction_coeffs: np.ndarray | None = None
    d_edt_ics: np.ndarray | None = None
    _weighted_pp_cache: OrderedDict[int, np.ndarray] | None = None
    _pp_below_grid_cache: OrderedDict[int, np.ndarray] | None = None
    _weighted_ics_gamma_cache: OrderedDict[int, np.ndarray] | None = None
    _weighted_ics_electron_cache: OrderedDict[int, np.ndarray] | None = None
    _ics_gamma_energy_cache: OrderedDict[int, np.ndarray] | None = None
    _ics_below_grid_cache: OrderedDict[int, np.ndarray] | None = None
    _ics_cel_cache: OrderedDict[tuple[int, float], tuple[np.ndarray, np.ndarray]] | None = None
    b_field_gauss: float = 0.0
    b_field_gamma: float = 0.0

    def reset(self) -> None:
        """Drop open files and cached tables so the next call reloads the bundle."""
        if self.runtime_ebl_file is not None:
            self.runtime_ebl_file.close()
        self.runtime_ebl_file = None
        self.manifest = None
        self.step_sizes = None
        self.row_ptr = None
        self.zreg_indices = None
        self.redshift_left_idx = None
        self.redshift_weights = None
        self.redshift_scales = None
        self.packed_row_ptr = None
        self.imfp = None
        self.extinction_coeffs = None
        self.log_extinction_coeffs = None
        self.attenuation_vectors = None
        self.ics_imfp = None
        self.ics_extinction_coeffs = None
        self.log_ics_extinction_coeffs = None
        self.d_edt_ics = None
        self._weighted_pp_cache = None
        self._pp_below_grid_cache = None
        self._weighted_ics_gamma_cache = None
        self._weighted_ics_electron_cache = None
        self._ics_gamma_energy_cache = None
        self._ics_below_grid_cache = None
        self._ics_cel_cache = None

    def ensure_bundle_loaded(self) -> None:
        """Load the manifest once and verify that the selected path is a valid bundle."""
        if self.manifest is not None:
            return
        if not bundle.is_bundle_root(self.library_path):
            raise FileNotFoundError(
                f"Could not find a GCascade HDF5 bundle at {self.library_path}. "
                f"{_library_path_config_message()}"
            )
        self.manifest = bundle.read_manifest(self.library_path)

    def ensure_common_loaded(self) -> None:
        """Load the shared redshift-window tables used by every EBL model."""
        self.ensure_bundle_loaded()
        if self.step_sizes is not None:
            return
        common_path = bundle.bundle_paths(self.library_path).runtime_common
        with h5py.File(common_path, "r") as handle:
            self.step_sizes = np.asarray(handle["step_sizes"], dtype=np.float64)
            self.row_ptr = np.asarray(handle["step_row_ptr"], dtype=np.int32)
            self.zreg_indices = np.asarray(handle["zreg_indices"], dtype=np.uint16)
            self.redshift_left_idx = np.asarray(handle["redshift_left_idx"], dtype=np.uint16)
            self.redshift_weights = np.asarray(handle["redshift_weights"], dtype=np.float64)
            self.redshift_scales = np.asarray(handle["redshift_scales"], dtype=np.float64)
            self.packed_row_ptr = np.asarray(handle["packed_row_ptr"], dtype=np.int32)

    def ensure_ebl_loaded(self) -> None:
        """Load the photon attenuation tables for the currently selected EBL model."""
        self.ensure_common_loaded()
        if self.imfp is not None and self.extinction_coeffs is not None and self.attenuation_vectors is not None:
            return

        runtime_path = bundle.resolve_runtime_ebl_path(self.library_path, self.ebl_index)
        self.runtime_ebl_file = h5py.File(runtime_path, "r")
        self.imfp = np.asarray(self.runtime_ebl_file["imfp"], dtype=np.float64)
        self.extinction_coeffs = np.asarray(self.runtime_ebl_file["extinction_coeffs"], dtype=np.float64)
        with np.errstate(divide="ignore"):
            self.log_extinction_coeffs = np.log(self.extinction_coeffs)
        self.attenuation_vectors = np.asarray(self.runtime_ebl_file["attenuation_vectors"], dtype=np.float64)

    def ensure_transport_loaded(self) -> None:
        """Load the pair-production and inverse-Compton kernels for electron tracking."""
        self.ensure_ebl_loaded()
        assert self.runtime_ebl_file is not None
        required = {
            "pp_packed",
            "ics_imfp",
            "ics_extinction_coeffs",
            "ics_gamma_packed",
            "ics_electron_packed",
            "dEdt_ics",
        }
        missing = sorted(required - set(self.runtime_ebl_file.keys()))
        if missing:
            raise RuntimeError(
                "The active GCascade library is missing electron-tracking tables. "
                f"Missing datasets in {self.runtime_ebl_file.filename}: {missing}"
            )
        if self.ics_imfp is None:
            self.ics_imfp = np.asarray(self.runtime_ebl_file["ics_imfp"], dtype=np.float64)
            self.ics_extinction_coeffs = np.asarray(self.runtime_ebl_file["ics_extinction_coeffs"], dtype=np.float64)
            with np.errstate(divide="ignore"):
                self.log_ics_extinction_coeffs = np.log(self.ics_extinction_coeffs)
            self.d_edt_ics = np.asarray(self.runtime_ebl_file["dEdt_ics"], dtype=np.float64)
            self._weighted_pp_cache = OrderedDict()
            self._pp_below_grid_cache = OrderedDict()
            self._weighted_ics_gamma_cache = OrderedDict()
            self._weighted_ics_electron_cache = OrderedDict()
            self._ics_gamma_energy_cache = OrderedDict()
            self._ics_below_grid_cache = OrderedDict()
            self._ics_cel_cache = OrderedDict()

    def set_ebl_index(self, ebl_index: int) -> None:
        """Switch to a different EBL model and clear all EBL-specific caches."""
        new_ebl = int(ebl_index)
        if new_ebl == self.ebl_index:
            return
        self.ebl_index = new_ebl
        if self.runtime_ebl_file is not None:
            self.runtime_ebl_file.close()
            self.runtime_ebl_file = None
        self.imfp = None
        self.extinction_coeffs = None
        self.log_extinction_coeffs = None
        self.attenuation_vectors = None
        self.ics_imfp = None
        self.ics_extinction_coeffs = None
        self.log_ics_extinction_coeffs = None
        self.d_edt_ics = None
        self._weighted_pp_cache = None
        self._pp_below_grid_cache = None
        self._weighted_ics_gamma_cache = None
        self._weighted_ics_electron_cache = None
        self._ics_gamma_energy_cache = None
        self._ics_below_grid_cache = None
        self._ics_cel_cache = None

    def _cache_transport_vector(self, cache: OrderedDict[int, np.ndarray], key: int, value: np.ndarray) -> np.ndarray:
        """Store one recently used transport slice in a tiny least-recently-used cache."""
        cache[key] = value
        while len(cache) > 4:
            cache.popitem(last=False)
        return value

    def _load_weighted_pp_slice(self, z_index: int) -> tuple[np.ndarray, np.ndarray]:
        """Load the pair-production electron yield kernel for one redshift layer."""
        self.ensure_transport_loaded()
        assert self.runtime_ebl_file is not None
        dense = bundle._unpack_lower_slice(np.asarray(self.runtime_ebl_file["pp_packed"][z_index], dtype=np.float64))
        row_energy = dense @ (energies * _TRAPEZOID_WEIGHTS)
        valid = np.logical_and(energies > 0.0, row_energy > 0.0)
        if np.any(valid):
            scale = 1.0 / float(np.median(row_energy[valid] / energies[valid]))
            dense *= scale
            row_energy *= scale
        weighted = np.ascontiguousarray(dense.T * _TRAPEZOID_WEIGHTS)
        below_grid = np.maximum(energies - row_energy, 0.0)
        return weighted, below_grid.astype(np.float64, copy=False)

    def _load_weighted_ics_slices(
        self,
        z_index: int,
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Load the photon and electron inverse-Compton kernels for one redshift layer."""
        self.ensure_transport_loaded()
        assert self.runtime_ebl_file is not None
        gamma_dense = bundle._unpack_lower_slice(
            np.asarray(self.runtime_ebl_file["ics_gamma_packed"][z_index], dtype=np.float64)
        )
        electron_dense = bundle._unpack_lower_slice(
            np.asarray(self.runtime_ebl_file["ics_electron_packed"][z_index], dtype=np.float64)
        )
        gamma_energy = gamma_dense @ (energies * _TRAPEZOID_WEIGHTS)
        electron_energy = electron_dense @ (energies * _TRAPEZOID_WEIGHTS)
        total_energy = gamma_energy + electron_energy
        valid = np.logical_and(energies > 0.0, total_energy > 0.0)
        if np.any(valid):
            scale = 1.0 / float(np.median(total_energy[valid] / energies[valid]))
            gamma_dense *= scale
            electron_dense *= scale
            gamma_energy *= scale
            electron_energy *= scale
            total_energy *= scale
        below_grid = np.maximum(energies - total_energy, 0.0)
        weighted_gamma = np.ascontiguousarray(gamma_dense.T * _TRAPEZOID_WEIGHTS)
        weighted_electron = np.ascontiguousarray(electron_dense.T * _TRAPEZOID_WEIGHTS)
        return (
            weighted_gamma,
            weighted_electron,
            gamma_energy.astype(np.float64, copy=False),
            below_grid.astype(np.float64, copy=False),
        )

    def get_weighted_pp_slice(self, z_index: int) -> np.ndarray:
        """Return the pair-production electron kernel weighted for fast matrix products."""
        if self._weighted_pp_cache is None:
            self._weighted_pp_cache = OrderedDict()
        if self._pp_below_grid_cache is None:
            self._pp_below_grid_cache = OrderedDict()
        if z_index not in self._weighted_pp_cache:
            weighted, below_grid = self._load_weighted_pp_slice(int(z_index))
            self._cache_transport_vector(self._weighted_pp_cache, int(z_index), weighted)
            self._cache_transport_vector(self._pp_below_grid_cache, int(z_index), below_grid)
        cached = self._weighted_pp_cache.pop(int(z_index))
        self._weighted_pp_cache[int(z_index)] = cached
        return cached

    def get_pp_below_grid_energy(self, z_index: int) -> np.ndarray:
        """Return the pair-production energy that falls below the tracked grid."""
        self.get_weighted_pp_slice(int(z_index))
        assert self._pp_below_grid_cache is not None
        cached = self._pp_below_grid_cache.pop(int(z_index))
        self._pp_below_grid_cache[int(z_index)] = cached
        return cached

    def get_weighted_ics_gamma_slice(self, z_index: int) -> np.ndarray:
        """Return the inverse-Compton photon yield kernel weighted for fast matrix products."""
        if self._weighted_ics_gamma_cache is None:
            self._weighted_ics_gamma_cache = OrderedDict()
        if self._weighted_ics_electron_cache is None:
            self._weighted_ics_electron_cache = OrderedDict()
        if self._ics_gamma_energy_cache is None:
            self._ics_gamma_energy_cache = OrderedDict()
        if self._ics_below_grid_cache is None:
            self._ics_below_grid_cache = OrderedDict()
        if z_index not in self._weighted_ics_gamma_cache:
            weighted_gamma, weighted_electron, gamma_energy, below_grid = self._load_weighted_ics_slices(int(z_index))
            self._cache_transport_vector(self._weighted_ics_gamma_cache, int(z_index), weighted_gamma)
            self._cache_transport_vector(self._weighted_ics_electron_cache, int(z_index), weighted_electron)
            self._cache_transport_vector(self._ics_gamma_energy_cache, int(z_index), gamma_energy)
            self._cache_transport_vector(self._ics_below_grid_cache, int(z_index), below_grid)
        cached = self._weighted_ics_gamma_cache.pop(int(z_index))
        self._weighted_ics_gamma_cache[int(z_index)] = cached
        return cached

    def get_weighted_ics_electron_slice(self, z_index: int) -> np.ndarray:
        """Return the inverse-Compton electron redistribution kernel."""
        self.get_weighted_ics_gamma_slice(int(z_index))
        assert self._weighted_ics_electron_cache is not None
        cached = self._weighted_ics_electron_cache.pop(int(z_index))
        self._weighted_ics_electron_cache[int(z_index)] = cached
        return cached

    def get_ics_gamma_row_energy(self, z_index: int) -> np.ndarray:
        """Return the mean photon energy emitted by one ICS interaction from each bin."""
        self.get_weighted_ics_gamma_slice(int(z_index))
        assert self._ics_gamma_energy_cache is not None
        cached = self._ics_gamma_energy_cache.pop(int(z_index))
        self._ics_gamma_energy_cache[int(z_index)] = cached
        return cached

    def get_ics_below_grid_energy(self, z_index: int) -> np.ndarray:
        """Return the ICS energy that leaves the tracked grid entirely."""
        self.get_weighted_ics_gamma_slice(int(z_index))
        assert self._ics_below_grid_cache is not None
        cached = self._ics_below_grid_cache.pop(int(z_index))
        self._ics_below_grid_cache[int(z_index)] = cached
        return cached

    def get_ics_cel_data(self, z_index: int, threshold_log_width: float) -> tuple[np.ndarray, np.ndarray]:
        """Mark the electron bins where ICS is treated as a continuous energy loss."""
        self.ensure_transport_loaded()
        key = (int(z_index), float(threshold_log_width))
        if self._ics_cel_cache is None:
            self._ics_cel_cache = OrderedDict()
        if key in self._ics_cel_cache:
            cached = self._ics_cel_cache.pop(key)
            self._ics_cel_cache[key] = cached
            return cached

        assert self.runtime_ebl_file is not None
        dense = bundle._unpack_lower_slice(
            np.asarray(self.runtime_ebl_file["ics_electron_packed"][int(z_index)], dtype=np.float64)
        )
        output_weights = dense * _TRAPEZOID_WEIGHTS[None, :]
        norm = np.sum(output_weights, axis=1)
        mean_energy = np.divide(
            output_weights @ energies,
            norm,
            out=energies.copy(),
            where=norm > 0.0,
        )
        mean_energy = np.minimum(mean_energy, energies)
        log_width = np.empty_like(energies)
        log_width[:-1] = np.diff(np.log(energies))
        log_width[-1] = log_width[-2]
        with np.errstate(divide="ignore", invalid="ignore"):
            loss_widths = np.log(energies / np.maximum(mean_energy, energies[0])) / log_width
        cel_mask = np.logical_and(norm > 0.0, loss_widths < float(threshold_log_width))
        result = (cel_mask.astype(bool), mean_energy.astype(np.float64, copy=False))
        self._ics_cel_cache[key] = result
        while len(self._ics_cel_cache) > 4:
            self._ics_cel_cache.popitem(last=False)
        return result

    def set_magnetic_field(self, b_field_gauss: float, gamma: float) -> None:
        """Store the magnetic-field law B(z)=B0(1+z)^gamma used by synchrotron cooling."""
        self.b_field_gauss = float(b_field_gauss)
        self.b_field_gamma = float(gamma)


_STATE: RuntimeBundleState | None = None
EBLindex = 1


def state() -> RuntimeBundleState:
    """Return the singleton runtime state used by the public API."""
    global _STATE, EBLindex
    if _STATE is None:
        _STATE = RuntimeBundleState(
            library_path=_discover_default_library_path(),
            ebl_index=EBLindex,
        )
    return _STATE


def reset_state() -> None:
    """Reset the global runtime state to its default configuration."""
    global _STATE, EBLindex
    if _STATE is not None:
        _STATE.reset()
    _STATE = None
    EBLindex = 1


def set_library_path(path: str | os.PathLike[str]) -> None:
    """Point the runtime at a different GCascadeV5 bundle directory."""
    global _STATE
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Library path does not exist: {resolved}\n{_library_path_config_message()}")
    if not resolved.is_dir():
        raise NotADirectoryError(f"Library path is not a directory: {resolved}")
    if _STATE is None:
        _STATE = RuntimeBundleState(resolved)
    else:
        _STATE.library_path = resolved
        _STATE.reset()


def get_library_path() -> Path:
    """Return the active bundle directory."""
    return state().library_path

