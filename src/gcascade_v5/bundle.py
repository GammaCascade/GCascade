from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

import h5py
import numpy as np
from scipy.io import loadmat, savemat

from . import legacy


BUNDLE_SCHEMA_VERSION = 1
MANIFEST_FILENAME = "bundle_manifest.json"
RUNTIME_DIRNAME = "runtime"
BUILDER_DIRNAME = "builder"
GENERATED_DIRNAME = "generated"
COMMON_RUNTIME_FILENAME = "common.h5"
PACKED_TRIANGULAR_SIZE = (len(legacy.energies) * (len(legacy.energies) + 1)) // 2
PACKED_ROW_PTR = np.cumsum(
    np.concatenate([np.array([0], dtype=np.int32), np.arange(1, len(legacy.energies) + 1, dtype=np.int32)])
)
PACKAGE_VERSION = "5.0"


@dataclass(frozen=True)
class BundlePaths:
    root: Path
    manifest: Path
    runtime_common: Path
    runtime_dir: Path
    builder_dir: Path
    generated_dir: Path


def bundle_paths(root: str | Path) -> BundlePaths:
    bundle_root = Path(root).expanduser().resolve()
    return BundlePaths(
        root=bundle_root,
        manifest=bundle_root / MANIFEST_FILENAME,
        runtime_common=bundle_root / RUNTIME_DIRNAME / COMMON_RUNTIME_FILENAME,
        runtime_dir=bundle_root / RUNTIME_DIRNAME,
        builder_dir=bundle_root / BUILDER_DIRNAME,
        generated_dir=bundle_root / GENERATED_DIRNAME,
    )


def is_bundle_root(path: str | Path) -> bool:
    return bundle_paths(path).manifest.exists()


def is_legacy_root(path: str | Path) -> bool:
    root = Path(path).expanduser().resolve()
    return (
        (root / "stepSizeArray.mat").exists()
        and (root / "zRegIndexArray.mat").exists()
        and (root / "cycle-spec").is_dir()
        and (root / "pp-IMFPs").is_dir()
    )


def read_manifest(root: str | Path) -> dict[str, Any]:
    manifest_path = bundle_paths(root).manifest
    if not manifest_path.exists():
        raise FileNotFoundError(
            f"Bundle manifest not found at {manifest_path}. "
            "Convert a legacy library first with gcascade_v5.convert_legacy_library(...)."
        )
    with manifest_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_manifest(root: str | Path, manifest: dict[str, Any]) -> None:
    manifest_path = bundle_paths(root).manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = manifest_path.with_suffix(".tmp")
    with temp_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temp_path.replace(manifest_path)


def bundle_info(root: str | Path) -> dict[str, Any]:
    manifest = read_manifest(root)
    return {
        "root": str(Path(root).expanduser().resolve()),
        "schema_version": manifest["schema_version"],
        "code_version": manifest["code_version"],
        "available_ebl_indices": [int(x) for x in manifest["available_ebl_indices"]],
        "active_generated_variants": manifest["generated"]["active"],
        "files": manifest["files"],
        "provenance": manifest["provenance"],
    }


def generated_variants(root: str | Path, ebl_index: int) -> list[dict[str, Any]]:
    manifest = read_manifest(root)
    return list(manifest["generated"]["variants"].get(str(int(ebl_index)), []))


def resolve_runtime_ebl_path(root: str | Path, ebl_index: int) -> Path:
    manifest = read_manifest(root)
    return (Path(root).expanduser().resolve() / manifest["files"]["runtime"][str(int(ebl_index))]).resolve()


def resolve_builder_ebl_path(root: str | Path, ebl_index: int) -> Path:
    manifest = read_manifest(root)
    return (Path(root).expanduser().resolve() / manifest["files"]["builder"][str(int(ebl_index))]).resolve()


def resolve_active_cycle_path(root: str | Path, ebl_index: int) -> Path:
    manifest = read_manifest(root)
    active = manifest["generated"]["active"].get(str(int(ebl_index)))
    if active is not None:
        return (Path(root).expanduser().resolve() / active["path"]).resolve()
    return resolve_runtime_ebl_path(root, ebl_index)


def _checksum_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _relative_to_root(path: Path, root: Path) -> str:
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _load_legacy_csv(path: Path) -> np.ndarray:
    try:
        return np.loadtxt(path, delimiter=",", dtype=np.float64)
    except ValueError:
        return np.loadtxt(path, dtype=np.float64)


def _load_legacy_mat(path: Path) -> np.ndarray:
    payload = loadmat(path)
    keys = [key for key in payload if not key.startswith("__")]
    if not keys:
        raise ValueError(f"No non-metadata arrays found in MAT file: {path}")
    return np.asarray(payload[keys[0]])


def _ensure_cycle_shape(arr: np.ndarray) -> np.ndarray:
    target = (len(legacy.zReg), len(legacy.energies), len(legacy.energies))
    if arr.shape == target:
        return arr
    legacy_shape = (len(legacy.energies), len(legacy.energies), len(legacy.zReg))
    if arr.shape == legacy_shape:
        return np.transpose(arr, (2, 0, 1))
    raise ValueError(f"Unexpected cycle table shape: {arr.shape}")


def _ensure_e_loss_shape(arr: np.ndarray) -> np.ndarray:
    target = (len(legacy.zReg), len(legacy.energies))
    if arr.shape == target:
        return arr
    if arr.shape == (len(legacy.energies), len(legacy.zReg)):
        return arr.T
    raise ValueError(f"Unexpected energy-loss table shape: {arr.shape}")


def _flatten_rows(
    step_array: np.ndarray,
    zreg_array: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if step_array.shape == (10000, len(legacy.diffuseDistances)):
        step_array = step_array.T
    if zreg_array.shape == (10000, len(legacy.diffuseDistances)):
        zreg_array = zreg_array.T

    if step_array.shape != (len(legacy.diffuseDistances), 10000):
        raise ValueError(f"Unexpected stepSizeArray shape: {step_array.shape}")
    if zreg_array.shape != (len(legacy.diffuseDistances), 10000):
        raise ValueError(f"Unexpected zRegIndexArray shape: {zreg_array.shape}")

    row_ptr = np.zeros(len(legacy.diffuseDistances) + 1, dtype=np.int32)
    step_rows: list[np.ndarray] = []
    zreg_rows: list[np.ndarray] = []
    cursor = 0
    for row_idx in range(len(legacy.diffuseDistances)):
        valid = step_array[row_idx] > 0.0
        step_row = np.asarray(step_array[row_idx, valid], dtype=np.float64)
        zreg_row = np.asarray(zreg_array[row_idx, valid], dtype=np.uint16) - 1
        step_rows.append(step_row)
        zreg_rows.append(zreg_row)
        cursor += int(step_row.size)
        row_ptr[row_idx + 1] = cursor

    flat_steps = np.concatenate(step_rows, dtype=np.float64)
    flat_zreg = np.concatenate(zreg_rows, dtype=np.uint16)
    return flat_steps, flat_zreg, row_ptr


def _build_redshift_tables() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n_windows = len(legacy.diffuseDistances)
    n_energies = len(legacy.energies)
    left_idx = np.empty((n_windows, n_energies), dtype=np.uint16)
    weights = np.empty((n_windows, n_energies), dtype=np.float64)
    scale = np.empty(n_windows, dtype=np.float64)

    energies = legacy.energies
    for window_idx in range(n_windows):
        z_hi = float(legacy.diffuseDistances[window_idx])
        z_lo = 0.0 if window_idx == 0 else float(legacy.diffuseDistances[window_idx - 1])
        stretch = (1.0 + z_hi) / (1.0 + z_lo)
        scale[window_idx] = stretch

        stretched = energies * stretch
        indices = np.searchsorted(energies, stretched, side="right") - 1
        indices = np.clip(indices, 0, len(energies) - 2)
        denom = energies[indices + 1] - energies[indices]
        local_weights = (stretched - energies[indices]) / denom

        left_idx[window_idx, :-1] = indices[:-1].astype(np.uint16)
        weights[window_idx, :-1] = local_weights[:-1]

        left_idx[window_idx, -1] = np.uint16(len(energies) - 2)
        weights[window_idx, -1] = 1.0

    return left_idx, weights, scale


def _pack_lower_cube(arr: np.ndarray) -> np.ndarray:
    if arr.shape[1:] != (len(legacy.energies), len(legacy.energies)):
        raise ValueError(f"Unexpected cube shape: {arr.shape}")
    packed = np.empty((arr.shape[0], PACKED_TRIANGULAR_SIZE), dtype=np.float64)
    for idx in range(arr.shape[0]):
        cursor = 0
        for row in range(len(legacy.energies)):
            size = row + 1
            packed[idx, cursor : cursor + size] = arr[idx, row, :size]
            cursor += size
    return packed


def _unpack_lower_slice(packed: np.ndarray) -> np.ndarray:
    out = np.zeros((len(legacy.energies), len(legacy.energies)), dtype=np.float64)
    cursor = 0
    for row in range(len(legacy.energies)):
        size = row + 1
        out[row, :size] = packed[cursor : cursor + size]
        cursor += size
    return out


def _write_h5_dataset(handle: h5py.File, name: str, array: np.ndarray, *, compression: str | None) -> None:
    handle.create_dataset(name, data=array, compression=compression, shuffle=True)


def _compute_attenuation_vectors(
    extinction_coeffs: np.ndarray,
    flat_steps: np.ndarray,
    flat_zreg: np.ndarray,
    row_ptr: np.ndarray,
) -> np.ndarray:
    with np.errstate(divide="ignore"):
        log_extinction = np.log(extinction_coeffs)
    attenuation = np.empty((len(legacy.diffuseDistances), len(legacy.energies)), dtype=np.float64)
    for row_idx in range(len(legacy.diffuseDistances)):
        start = int(row_ptr[row_idx])
        stop = int(row_ptr[row_idx + 1])
        if start == stop:
            attenuation[row_idx] = np.ones(len(legacy.energies), dtype=np.float64)
            continue
        step_sizes = flat_steps[start:stop]
        z_indices = flat_zreg[start:stop].astype(np.int64, copy=False)
        attenuation[row_idx] = np.exp(np.sum(log_extinction[z_indices] * step_sizes[:, None], axis=0))
    return attenuation


def _build_runtime_ebl_file(
    source_root: Path,
    target_root: Path,
    ebl_index: int,
    flat_steps: np.ndarray,
    flat_zreg: np.ndarray,
    row_ptr: np.ndarray,
) -> Path:
    ebl_name = legacy.EBL_NAME_MAP[ebl_index]
    imfp_cmb = _load_legacy_csv(source_root / "pp-IMFPs" / "IMFPcmb.csv")
    if ebl_index == 0:
        imfp_ebl = np.zeros_like(imfp_cmb, dtype=np.float64)
    else:
        imfp_ebl = _load_legacy_csv(source_root / "pp-IMFPs" / f"IMFPebl{ebl_name}.csv")

    imfp = np.asarray(imfp_cmb + imfp_ebl, dtype=np.float64)
    extinction = np.asarray(np.exp(-legacy.Mpc * imfp), dtype=np.float64)
    attenuation_vectors = _compute_attenuation_vectors(extinction, flat_steps, flat_zreg, row_ptr)

    cycle_raw = _ensure_cycle_shape(
        np.asarray(
            _load_legacy_mat(source_root / "cycle-spec" / f"cyclespec{ebl_name}.mat"),
            dtype=np.float64,
        )
    )
    runtime_cycle = np.ascontiguousarray(cycle_raw * 1.0e9, dtype=np.float64)
    packed_cycle = _pack_lower_cube(runtime_cycle)

    runtime_path = target_root / RUNTIME_DIRNAME / f"ebl_{ebl_name}.h5"
    runtime_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(runtime_path, "w") as handle:
        _write_h5_dataset(handle, "imfp", imfp, compression="lzf")
        _write_h5_dataset(handle, "extinction_coeffs", extinction, compression="lzf")
        _write_h5_dataset(handle, "attenuation_vectors", attenuation_vectors, compression="lzf")
        _write_h5_dataset(handle, "cycle_packed", packed_cycle, compression="lzf")
        handle.attrs["ebl_index"] = int(ebl_index)
        handle.attrs["ebl_name"] = ebl_name
        handle.attrs["packed_triangular_size"] = PACKED_TRIANGULAR_SIZE
    return runtime_path


def _build_builder_ebl_file(source_root: Path, target_root: Path, ebl_index: int) -> Path:
    ebl_name = legacy.EBL_NAME_MAP[ebl_index]
    pp_raw = _ensure_cycle_shape(
        np.asarray(
            _load_legacy_mat(source_root / "pp-spec" / f"normalizedPPspec{ebl_name}.mat"),
            dtype=np.float64,
        )
    )
    ots_raw = _ensure_cycle_shape(
        np.asarray(
            _load_legacy_mat(source_root / "on-the-spot-ics-spec" / f"onthespotICSspec{ebl_name}.mat"),
            dtype=np.float64,
        )
    )
    d_edt = _ensure_e_loss_shape(
        np.asarray(
            _load_legacy_mat(source_root / "ics-E-loss-rates" / f"dEdt{ebl_name}.mat"),
            dtype=np.float64,
        )
    )

    builder_path = target_root / BUILDER_DIRNAME / f"ebl_{ebl_name}.h5"
    builder_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(builder_path, "w") as handle:
        _write_h5_dataset(handle, "pp_packed", _pack_lower_cube(pp_raw), compression="gzip")
        _write_h5_dataset(handle, "ots_packed", _pack_lower_cube(ots_raw), compression="gzip")
        _write_h5_dataset(handle, "dEdt_ics", d_edt, compression="gzip")
        handle.attrs["ebl_index"] = int(ebl_index)
        handle.attrs["ebl_name"] = ebl_name
        handle.attrs["packed_triangular_size"] = PACKED_TRIANGULAR_SIZE
    return builder_path


def _import_extra_cycle_variants(
    *,
    bundle_root: Path,
    source_dir: Path,
    indices: list[int],
    active_standard_variants: set[Path],
) -> None:
    cycle_dir = source_dir / "cycle-spec"
    if not cycle_dir.is_dir():
        return
    for mat_file in sorted(cycle_dir.glob("*.mat")):
        matched_ebl: int | None = None
        matched_name = ""
        for ebl_index in indices:
            ebl_name = legacy.EBL_NAME_MAP[ebl_index]
            if ebl_name in mat_file.stem:
                matched_ebl = ebl_index
                matched_name = ebl_name
                break
        if matched_ebl is None:
            continue
        standard_name = f"cyclespec{matched_name}.mat"
        if source_dir.name.endswith("-legacy") and mat_file.name == standard_name:
            continue
        if mat_file in active_standard_variants:
            continue
        import_cycle_variant_mat(
            bundle_root,
            mat_file,
            variant_id=f"imported_{mat_file.stem}",
            active=False,
        )


def convert_legacy_library(
    source_path: str | Path,
    target_path: str | Path,
    *,
    overwrite: bool = False,
    ebl_indices: list[int] | tuple[int, ...] | None = None,
    include_builder: bool = True,
    generated_source_path: str | Path | None = None,
) -> Path:
    source_root = Path(source_path).expanduser().resolve()
    target_root = Path(target_path).expanduser().resolve()

    if not is_legacy_root(source_root):
        raise FileNotFoundError(
            f"Legacy GCascade library not found or incomplete at {source_root}. "
            "Expected MAT/CSV tables such as stepSizeArray.mat and cycle-spec/."
        )

    if target_root.exists():
        if not overwrite:
            raise FileExistsError(f"Target bundle path already exists: {target_root}")
        shutil.rmtree(target_root)

    paths = bundle_paths(target_root)
    paths.runtime_dir.mkdir(parents=True, exist_ok=True)
    if include_builder:
        paths.builder_dir.mkdir(parents=True, exist_ok=True)
    paths.generated_dir.mkdir(parents=True, exist_ok=True)

    step_array = np.asarray(_load_legacy_mat(source_root / "stepSizeArray.mat"), dtype=np.float64)
    zreg_array = np.asarray(_load_legacy_mat(source_root / "zRegIndexArray.mat"), dtype=np.uint16)
    flat_steps, flat_zreg, row_ptr = _flatten_rows(step_array, zreg_array)
    left_idx, weights, scales = _build_redshift_tables()

    with h5py.File(paths.runtime_common, "w") as handle:
        _write_h5_dataset(handle, "energies", np.asarray(legacy.energies, dtype=np.float64), compression=None)
        _write_h5_dataset(handle, "diffuse_steps", np.asarray(legacy.diffuseSteps, dtype=np.float64), compression=None)
        _write_h5_dataset(
            handle,
            "diffuse_distances",
            np.asarray(legacy.diffuseDistances, dtype=np.float64),
            compression=None,
        )
        _write_h5_dataset(handle, "z_reg", np.asarray(legacy.zReg, dtype=np.float64), compression=None)
        _write_h5_dataset(handle, "step_sizes", flat_steps, compression="lzf")
        _write_h5_dataset(handle, "step_row_ptr", row_ptr, compression=None)
        _write_h5_dataset(handle, "zreg_indices", flat_zreg, compression="lzf")
        _write_h5_dataset(handle, "zreg_row_ptr", row_ptr, compression=None)
        _write_h5_dataset(handle, "redshift_left_idx", left_idx, compression="lzf")
        _write_h5_dataset(handle, "redshift_weights", weights, compression="lzf")
        _write_h5_dataset(handle, "redshift_scales", scales, compression="lzf")
        _write_h5_dataset(handle, "packed_row_ptr", PACKED_ROW_PTR, compression=None)

    indices = [int(idx) for idx in (ebl_indices if ebl_indices is not None else sorted(legacy.EBL_NAME_MAP))]
    runtime_files: dict[str, str] = {}
    builder_files: dict[str, str] = {}
    checksum_files: dict[str, str] = {}

    checksum_files[_relative_to_root(paths.runtime_common, target_root)] = _checksum_file(paths.runtime_common)
    for ebl_index in indices:
        runtime_path = _build_runtime_ebl_file(source_root, target_root, ebl_index, flat_steps, flat_zreg, row_ptr)
        runtime_files[str(ebl_index)] = _relative_to_root(runtime_path, target_root)
        checksum_files[runtime_files[str(ebl_index)]] = _checksum_file(runtime_path)

        if include_builder:
            builder_path = _build_builder_ebl_file(source_root, target_root, ebl_index)
            builder_files[str(ebl_index)] = _relative_to_root(builder_path, target_root)
            checksum_files[builder_files[str(ebl_index)]] = _checksum_file(builder_path)

    manifest = {
        "schema_version": BUNDLE_SCHEMA_VERSION,
        "code_version": PACKAGE_VERSION,
        "available_ebl_indices": indices,
        "ebl_name_map": {str(key): value for key, value in legacy.EBL_NAME_MAP.items()},
        "ebl_description_map": {str(key): value for key, value in legacy.EBL_DESCRIPTION_MAP.items()},
        "files": {
            "runtime_common": _relative_to_root(paths.runtime_common, target_root),
            "runtime": runtime_files,
            "builder": builder_files,
        },
        "generated": {"active": {}, "variants": {}},
        "checksums": checksum_files,
        "provenance": [
            {
                "kind": "legacy-conversion",
                "source_path": str(source_root),
                "created_at": _now_iso(),
                "include_builder": bool(include_builder),
                "ebl_indices": indices,
            }
        ],
    }
    write_manifest(target_root, manifest)

    generated_source_root: Path | None = None
    if generated_source_path is not None:
        generated_source_root = Path(generated_source_path).expanduser().resolve()
    else:
        default_generated = legacy._discover_default_generated_library_path()
        if default_generated.exists():
            generated_source_root = default_generated

    imported_active_sources: set[Path] = set()
    if generated_source_root is not None and generated_source_root.exists():
        generated_cycle_dir = generated_source_root / "cycle-spec"
        if generated_cycle_dir.is_dir():
            for ebl_index in indices:
                ebl_name = legacy.EBL_NAME_MAP[ebl_index]
                legacy_generated_path = generated_cycle_dir / f"cyclespec{ebl_name}.mat"
                if not legacy_generated_path.exists():
                    continue
                raw = _ensure_cycle_shape(np.asarray(_load_legacy_mat(legacy_generated_path), dtype=np.float64))
                runtime_cycle = np.ascontiguousarray(raw * 1.0e9, dtype=np.float64)
                target_generated = paths.generated_dir / f"imported_legacy_{ebl_name}.h5"
                with h5py.File(target_generated, "w") as handle:
                    _write_h5_dataset(handle, "cycle_packed", _pack_lower_cube(runtime_cycle), compression="lzf")
                    handle.attrs["ebl_index"] = int(ebl_index)
                    handle.attrs["ebl_name"] = ebl_name
                    handle.attrs["variant_id"] = f"imported_legacy_{ebl_name}"
                    handle.attrs["source_path"] = str(legacy_generated_path)
                    handle.attrs["created_at"] = _now_iso()
                register_generated_variant(
                    target_root,
                    ebl_index=ebl_index,
                    file_path=target_generated,
                    parameters={"imported_from": str(legacy_generated_path)},
                    variant_id=f"imported_legacy_{ebl_name}",
                )
                imported_active_sources.add(legacy_generated_path)

    _import_extra_cycle_variants(
        bundle_root=target_root,
        source_dir=source_root,
        indices=indices,
        active_standard_variants=imported_active_sources,
    )
    if generated_source_root is not None and generated_source_root.exists():
        _import_extra_cycle_variants(
            bundle_root=target_root,
            source_dir=generated_source_root,
            indices=indices,
            active_standard_variants=imported_active_sources,
        )
    return target_root


def register_generated_variant(
    bundle_root: str | Path,
    *,
    ebl_index: int,
    file_path: Path,
    parameters: dict[str, Any],
    variant_id: str,
    active: bool = True,
) -> None:
    root = Path(bundle_root).expanduser().resolve()
    manifest = read_manifest(root)
    key = str(int(ebl_index))
    rel_path = _relative_to_root(file_path, root)
    record = {
        "variant_id": variant_id,
        "path": rel_path,
        "parameters": parameters,
        "created_at": _now_iso(),
    }
    variants = list(manifest["generated"]["variants"].get(key, []))
    variants = [item for item in variants if item["variant_id"] != variant_id]
    variants.append(record)
    manifest["generated"]["variants"][key] = variants
    if active:
        manifest["generated"]["active"][key] = record
    manifest["checksums"][rel_path] = _checksum_file(file_path)
    manifest["provenance"].append(
        {
            "kind": "generated-cycle",
            "ebl_index": int(ebl_index),
            "variant_id": variant_id,
            "parameters": parameters,
            "path": rel_path,
            "created_at": _now_iso(),
        }
    )
    write_manifest(root, manifest)


def clear_active_generated_variant(bundle_root: str | Path, ebl_index: int | None = None) -> None:
    root = Path(bundle_root).expanduser().resolve()
    manifest = read_manifest(root)
    if ebl_index is None:
        manifest["generated"]["active"] = {}
    else:
        manifest["generated"]["active"].pop(str(int(ebl_index)), None)
    write_manifest(root, manifest)


def activate_cycle_path(
    bundle_root: str | Path,
    file_path: str | Path,
    *,
    ebl_index: int | None = None,
) -> Path:
    root = Path(bundle_root).expanduser().resolve()
    source = Path(file_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Cycle file does not exist: {source}")

    with h5py.File(source, "r") as handle:
        if "cycle_packed" not in handle:
            raise ValueError(f"Cycle file does not contain a 'cycle_packed' dataset: {source}")
        file_ebl_index = handle.attrs.get("ebl_index")
        file_variant_id = handle.attrs.get("variant_id", source.stem)

    matched_ebl = int(ebl_index) if ebl_index is not None else None
    if file_ebl_index is not None:
        detected = int(file_ebl_index)
        if matched_ebl is None:
            matched_ebl = detected
        elif matched_ebl != detected:
            raise ValueError(
                f"Requested ebl_index={matched_ebl} does not match file attribute ebl_index={detected}: {source}"
            )
    if matched_ebl is None:
        raise ValueError("Could not determine EBL index for cycle file. Provide ebl_index explicitly.")

    register_generated_variant(
        root,
        ebl_index=matched_ebl,
        file_path=source,
        parameters={"activated_from": str(source)},
        variant_id=str(file_variant_id),
        active=True,
    )
    return source


def infer_cycle_ebl_index(file_path: str | Path) -> int:
    source = Path(file_path).expanduser().resolve()
    if not source.exists():
        raise FileNotFoundError(f"Cycle file does not exist: {source}")
    with h5py.File(source, "r") as handle:
        if "cycle_packed" not in handle:
            raise ValueError(f"Cycle file does not contain a 'cycle_packed' dataset: {source}")
        file_ebl_index = handle.attrs.get("ebl_index")
    if file_ebl_index is None:
        raise ValueError(f"Cycle file does not declare an ebl_index attribute: {source}")
    return int(file_ebl_index)


def export_active_cycle_to_legacy_mat(
    bundle_root: str | Path,
    ebl_index: int,
    target_path: str | Path,
) -> Path:
    root = Path(bundle_root).expanduser().resolve()
    target = Path(target_path).expanduser().resolve()
    cycle_path = resolve_active_cycle_path(root, ebl_index)
    with h5py.File(cycle_path, "r") as handle:
        packed = np.asarray(handle["cycle_packed"], dtype=np.float64)
    dense = np.empty((packed.shape[0], len(legacy.energies), len(legacy.energies)), dtype=np.float64)
    for idx in range(packed.shape[0]):
        dense[idx] = _unpack_lower_slice(packed[idx])
    savemat(target, {"Expression1": dense / 1.0e9}, do_compression=True)
    return target


def import_cycle_variant_mat(
    bundle_root: str | Path,
    mat_path: str | Path,
    *,
    variant_id: str | None = None,
    active: bool = False,
) -> Path:
    root = Path(bundle_root).expanduser().resolve()
    source = Path(mat_path).expanduser().resolve()
    raw = _ensure_cycle_shape(np.asarray(_load_legacy_mat(source), dtype=np.float64))
    runtime_cycle = np.ascontiguousarray(raw * 1.0e9, dtype=np.float64)

    stem = source.stem
    matched_ebl: int | None = None
    matched_name = ""
    for ebl_index, ebl_name in legacy.EBL_NAME_MAP.items():
        if ebl_name in stem:
            matched_ebl = int(ebl_index)
            matched_name = ebl_name
            break
    if matched_ebl is None:
        raise ValueError(f"Could not infer EBL model from cycle file name: {source.name}")

    safe_variant = variant_id or stem.replace(" ", "_")
    target = bundle_paths(root).generated_dir / f"{safe_variant}.h5"
    target.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(target, "w") as handle:
        _write_h5_dataset(handle, "cycle_packed", _pack_lower_cube(runtime_cycle), compression="lzf")
        handle.attrs["ebl_index"] = matched_ebl
        handle.attrs["ebl_name"] = matched_name
        handle.attrs["variant_id"] = safe_variant
        handle.attrs["source_path"] = str(source)
        handle.attrs["created_at"] = _now_iso()

    register_generated_variant(
        root,
        ebl_index=matched_ebl,
        file_path=target,
        parameters={"imported_from": str(source)},
        variant_id=safe_variant,
        active=active,
    )
    return target
