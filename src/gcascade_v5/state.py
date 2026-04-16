from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass
import os
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from . import config
from . import bundle, legacy


LIB_PATH_ENV_VAR = legacy.LIB_PATH_ENV_VAR
GENERATED_LIB_PATH_ENV_VAR = legacy.GENERATED_LIB_PATH_ENV_VAR
_TRAPEZOID_WEIGHTS = np.empty(len(legacy.energies), dtype=np.float64)
_TRAPEZOID_WEIGHTS[0] = legacy.dEnergiesGamma[0]
_TRAPEZOID_WEIGHTS[-1] = legacy.dEnergiesGamma[-1]
_TRAPEZOID_WEIGHTS[1:-1] = legacy.dEnergiesGamma[:-1] + legacy.dEnergiesGamma[1:]


def _library_path_config_message() -> str:
    return (
        f"Configure with {LIB_PATH_ENV_VAR} (before import) or call "
        "gcascade_v5.set_library_path('/path/to/gcascade_bundle')."
    )


def _generated_library_path_config_message() -> str:
    return (
        f"Configure with {GENERATED_LIB_PATH_ENV_VAR} (before import) or call "
        "gcascade_v5.set_generated_library_path('/path/to/generated')."
    )


def _discover_default_library_path() -> Path:
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


def _discover_default_generated_library_path(library_path: Path | None = None) -> Path:
    env_path = os.getenv(GENERATED_LIB_PATH_ENV_VAR)
    if env_path:
        return Path(env_path).expanduser().resolve()
    if library_path is not None and bundle.is_bundle_root(library_path):
        return (Path(library_path).expanduser().resolve() / bundle.GENERATED_DIRNAME).resolve()
    return (Path.cwd() / "generated_libraries").resolve()


@dataclass
class RuntimeBundleState:
    library_path: Path
    generated_library_path: Path
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
    cycle_file: h5py.File | None = None
    active_cycle_path: Path | None = None
    imfp: np.ndarray | None = None
    extinction_coeffs: np.ndarray | None = None
    log_extinction_coeffs: np.ndarray | None = None
    attenuation_vectors: np.ndarray | None = None
    cycle_packed_array: np.ndarray | None = None
    _cycle_cache: OrderedDict[int, np.ndarray] | None = None
    _dense_cycle_cache: OrderedDict[int, np.ndarray] | None = None
    _weighted_cycle_cache: OrderedDict[int, np.ndarray] | None = None

    def reset(self) -> None:
        if self.runtime_ebl_file is not None:
            self.runtime_ebl_file.close()
        if self.cycle_file is not None and self.cycle_file is not self.runtime_ebl_file:
            self.cycle_file.close()
        self.runtime_ebl_file = None
        self.cycle_file = None
        self.active_cycle_path = None
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
        self.cycle_packed_array = None
        self._cycle_cache = None
        self._dense_cycle_cache = None
        self._weighted_cycle_cache = None

    def ensure_bundle_loaded(self) -> None:
        if self.manifest is not None:
            return
        if not bundle.is_bundle_root(self.library_path):
            if bundle.is_legacy_root(self.library_path):
                raise RuntimeError(
                    "The active library path points to a legacy MAT/CSV library. "
                    "Convert it first with gcascade_v5.convert_legacy_library(source_path, target_path). "
                    f"Active path: {self.library_path}"
                )
            raise FileNotFoundError(
                f"Could not find a GCascade HDF5 bundle at {self.library_path}. "
                f"{_library_path_config_message()}"
            )
        self.manifest = bundle.read_manifest(self.library_path)

    def ensure_common_loaded(self) -> None:
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

    def ensure_ebl_loaded(self, *, refresh_active_cycle: bool = False) -> None:
        self.ensure_common_loaded()
        if self.imfp is not None and self.extinction_coeffs is not None and self.attenuation_vectors is not None:
            if refresh_active_cycle:
                active_path = bundle.resolve_active_cycle_path(self.library_path, self.ebl_index)
                if self.active_cycle_path != active_path:
                    self._reload_cycle_source(active_path)
            return

        runtime_path = bundle.resolve_runtime_ebl_path(self.library_path, self.ebl_index)
        self.runtime_ebl_file = h5py.File(runtime_path, "r")
        self.imfp = np.asarray(self.runtime_ebl_file["imfp"], dtype=np.float64)
        self.extinction_coeffs = np.asarray(self.runtime_ebl_file["extinction_coeffs"], dtype=np.float64)
        with np.errstate(divide="ignore"):
            self.log_extinction_coeffs = np.log(self.extinction_coeffs)
        self.attenuation_vectors = np.asarray(self.runtime_ebl_file["attenuation_vectors"], dtype=np.float64)
        self._reload_cycle_source(bundle.resolve_active_cycle_path(self.library_path, self.ebl_index))

    def _reload_cycle_source(self, cycle_path: Path) -> None:
        if self.cycle_file is not None and self.cycle_file is not self.runtime_ebl_file:
            self.cycle_file.close()
        if self.runtime_ebl_file is not None and cycle_path.resolve() == Path(self.runtime_ebl_file.filename).resolve():
            self.cycle_file = self.runtime_ebl_file
        else:
            self.cycle_file = h5py.File(cycle_path, "r")
        self.active_cycle_path = cycle_path
        self._cycle_cache = OrderedDict()
        self._dense_cycle_cache = OrderedDict()
        self._weighted_cycle_cache = OrderedDict()
        self.cycle_packed_array = None
        if config.get_numba_enabled():
            self.cycle_packed_array = np.ascontiguousarray(self.cycle_file["cycle_packed"], dtype=np.float64)

    def set_ebl_index(self, ebl_index: int) -> None:
        new_ebl = int(ebl_index)
        if new_ebl == self.ebl_index:
            return
        self.ebl_index = new_ebl
        if self.runtime_ebl_file is not None:
            self.runtime_ebl_file.close()
            self.runtime_ebl_file = None
        if self.cycle_file is not None and self.cycle_file is not self.runtime_ebl_file:
            self.cycle_file.close()
            self.cycle_file = None
        self.imfp = None
        self.extinction_coeffs = None
        self.log_extinction_coeffs = None
        self.attenuation_vectors = None
        self.cycle_packed_array = None
        self._cycle_cache = None
        self._dense_cycle_cache = None
        self._weighted_cycle_cache = None
        self.active_cycle_path = None

    def get_row_bounds(self, window_idx: int) -> tuple[int, int]:
        self.ensure_common_loaded()
        assert self.row_ptr is not None
        return int(self.row_ptr[window_idx]), int(self.row_ptr[window_idx + 1])

    def get_cycle_slice(self, z_index: int) -> np.ndarray:
        if self.cycle_file is None:
            self.ensure_ebl_loaded()
        if self.cycle_packed_array is not None:
            return self.cycle_packed_array[z_index]
        assert self.cycle_file is not None
        assert self._cycle_cache is not None
        if z_index in self._cycle_cache:
            cached = self._cycle_cache.pop(z_index)
            self._cycle_cache[z_index] = cached
            return cached
        loaded = np.asarray(self.cycle_file["cycle_packed"][z_index], dtype=np.float64)
        self._cycle_cache[z_index] = loaded
        while len(self._cycle_cache) > 4:
            self._cycle_cache.popitem(last=False)
        return loaded

    def get_dense_cycle_slice(self, z_index: int) -> np.ndarray:
        if self.cycle_file is None:
            self.ensure_ebl_loaded()
        if self._dense_cycle_cache is None:
            self._dense_cycle_cache = OrderedDict()
        if z_index in self._dense_cycle_cache:
            cached = self._dense_cycle_cache.pop(z_index)
            self._dense_cycle_cache[z_index] = cached
            return cached
        dense = bundle._unpack_lower_slice(self.get_cycle_slice(z_index))
        self._dense_cycle_cache[z_index] = dense
        while len(self._dense_cycle_cache) > 4:
            self._dense_cycle_cache.popitem(last=False)
        return dense

    def get_weighted_cycle_slice(self, z_index: int) -> np.ndarray:
        if self.cycle_file is None:
            self.ensure_ebl_loaded()
        if self._weighted_cycle_cache is None:
            self._weighted_cycle_cache = OrderedDict()
        if z_index in self._weighted_cycle_cache:
            cached = self._weighted_cycle_cache.pop(z_index)
            self._weighted_cycle_cache[z_index] = cached
            return cached
        dense = bundle._unpack_lower_slice(self.get_cycle_slice(z_index))
        weighted = np.ascontiguousarray(dense.T * _TRAPEZOID_WEIGHTS)
        self._weighted_cycle_cache[z_index] = weighted
        while len(self._weighted_cycle_cache) > 4:
            self._weighted_cycle_cache.popitem(last=False)
        return weighted


_STATE: RuntimeBundleState | None = None
EBLindex = 1


def state() -> RuntimeBundleState:
    global _STATE, EBLindex
    if _STATE is None:
        library_path = _discover_default_library_path()
        _STATE = RuntimeBundleState(
            library_path=library_path,
            generated_library_path=_discover_default_generated_library_path(library_path),
            ebl_index=EBLindex,
        )
    return _STATE


def reset_state() -> None:
    global _STATE, EBLindex
    if _STATE is not None:
        _STATE.reset()
    _STATE = None
    EBLindex = 1


def set_library_path(path: str | os.PathLike[str]) -> None:
    global _STATE
    resolved = Path(path).expanduser().resolve()
    if not resolved.exists():
        raise FileNotFoundError(f"Library path does not exist: {resolved}\n{_library_path_config_message()}")
    if not resolved.is_dir():
        raise NotADirectoryError(f"Library path is not a directory: {resolved}")
    if _STATE is None:
        _STATE = RuntimeBundleState(resolved, _discover_default_generated_library_path(resolved))
    else:
        old_library = _STATE.library_path
        old_generated_default = _discover_default_generated_library_path(old_library)
        old_generated_path = _STATE.generated_library_path
        _STATE.library_path = resolved
        if old_generated_path == old_generated_default:
            _STATE.generated_library_path = _discover_default_generated_library_path(resolved)
        _STATE.reset()


def set_generated_library_path(path: str | os.PathLike[str]) -> None:
    global _STATE
    resolved = Path(path).expanduser().resolve()
    if resolved.exists() and not resolved.is_dir():
        raise NotADirectoryError(f"Generated library path is not a directory: {resolved}")
    if _STATE is None:
        library_path = _discover_default_library_path()
        _STATE = RuntimeBundleState(library_path, resolved)
    else:
        _STATE.generated_library_path = resolved


def get_library_path() -> Path:
    return state().library_path


def get_generated_library_path() -> Path:
    return state().generated_library_path
