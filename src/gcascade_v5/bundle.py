from __future__ import annotations

"""Helpers for reading the modern GCascadeV5 runtime bundle."""

from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np

from .physics import EBL_NAME_MAP, diffuseDistances, energies


BUNDLE_SCHEMA_VERSION = 2
MANIFEST_FILENAME = "bundle_manifest.json"
RUNTIME_DIRNAME = "runtime"
COMMON_RUNTIME_FILENAME = "common.h5"
PACKED_TRIANGULAR_SIZE = (len(energies) * (len(energies) + 1)) // 2
PACKED_ROW_PTR = np.cumsum(
    np.concatenate([np.array([0], dtype=np.int32), np.arange(1, len(energies) + 1, dtype=np.int32)])
)
PACKAGE_VERSION = "5.1"


@dataclass(frozen=True)
class BundlePaths:
    """Canonical file locations inside a modern GCascadeV5 library bundle."""

    root: Path
    manifest: Path
    runtime_common: Path
    runtime_dir: Path


def bundle_paths(root: str | Path) -> BundlePaths:
    """Return the standard manifest and runtime paths for a bundle root."""
    bundle_root = Path(root).expanduser().resolve()
    return BundlePaths(
        root=bundle_root,
        manifest=bundle_root / MANIFEST_FILENAME,
        runtime_common=bundle_root / RUNTIME_DIRNAME / COMMON_RUNTIME_FILENAME,
        runtime_dir=bundle_root / RUNTIME_DIRNAME,
    )


def is_bundle_root(path: str | Path) -> bool:
    """Check whether a directory looks like a GCascadeV5 runtime bundle."""
    return bundle_paths(path).manifest.exists()


def read_manifest(root: str | Path) -> dict[str, Any]:
    """Load the bundle manifest JSON from disk."""
    manifest_path = bundle_paths(root).manifest
    if not manifest_path.exists():
        raise FileNotFoundError(f"Bundle manifest not found at {manifest_path}.")
    with manifest_path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_manifest(root: str | Path, manifest: dict[str, Any]) -> None:
    """Write a manifest JSON file atomically."""
    manifest_path = bundle_paths(root).manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = manifest_path.with_suffix(".tmp")
    with temp_path.open("w", encoding="utf-8") as handle:
        json.dump(manifest, handle, indent=2, sort_keys=True)
        handle.write("\n")
    temp_path.replace(manifest_path)


def _checksum_file(path: Path) -> str:
    """Return a SHA256 checksum for a bundle file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def _relative_to_root(path: Path, root: Path) -> str:
    """Express a path relative to the bundle root when possible."""
    try:
        return path.resolve().relative_to(root.resolve()).as_posix()
    except ValueError:
        return path.resolve().as_posix()


def bundle_info(root: str | Path) -> dict[str, Any]:
    """Summarize the current runtime bundle in a user-facing dictionary."""
    manifest = read_manifest(root)
    return {
        "root": str(Path(root).expanduser().resolve()),
        "schema_version": int(manifest["schema_version"]),
        "code_version": manifest["code_version"],
        "available_ebl_indices": [int(x) for x in manifest["available_ebl_indices"]],
        "files": {
            "runtime_common": manifest["files"]["runtime_common"],
            "runtime": dict(manifest["files"]["runtime"]),
        },
        "transport_tables": manifest.get("transport_tables"),
        "provenance": manifest.get("provenance", []),
    }


def resolve_runtime_ebl_path(root: str | Path, ebl_index: int) -> Path:
    """Return the HDF5 runtime file used for one EBL background model."""
    manifest = read_manifest(root)
    return (Path(root).expanduser().resolve() / manifest["files"]["runtime"][str(int(ebl_index))]).resolve()


def _build_redshift_tables() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Precompute the interpolation indices and weights for one-window redshifting."""
    n_windows = len(diffuseDistances)
    n_energies = len(energies)
    left_idx = np.empty((n_windows, n_energies), dtype=np.uint16)
    weights = np.empty((n_windows, n_energies), dtype=np.float64)
    scale = np.empty(n_windows, dtype=np.float64)

    for window_idx in range(n_windows):
        z_hi = float(diffuseDistances[window_idx])
        z_lo = 0.0 if window_idx == 0 else float(diffuseDistances[window_idx - 1])
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
    """Pack lower-triangular transport kernels into the compact HDF5 layout."""
    if arr.shape[1:] != (len(energies), len(energies)):
        raise ValueError(f"Unexpected cube shape: {arr.shape}")
    packed = np.empty((arr.shape[0], PACKED_TRIANGULAR_SIZE), dtype=np.float64)
    for idx in range(arr.shape[0]):
        cursor = 0
        for row in range(len(energies)):
            size = row + 1
            packed[idx, cursor : cursor + size] = arr[idx, row, :size]
            cursor += size
    return packed


def _unpack_lower_slice(packed: np.ndarray) -> np.ndarray:
    """Reconstruct one lower-triangular transport kernel from its packed form."""
    out = np.zeros((len(energies), len(energies)), dtype=np.float64)
    cursor = 0
    for row in range(len(energies)):
        size = row + 1
        out[row, :size] = packed[cursor : cursor + size]
        cursor += size
    return out


def prune_bundle_artifacts(root: str | Path, *, delete_files: bool = True) -> Path:
    """Remove unused non-runtime metadata and optional files from a bundle."""
    bundle_root = Path(root).expanduser().resolve()
    manifest = read_manifest(bundle_root)

    files = dict(manifest.get("files", {}))
    files.pop("builder", None)
    manifest["files"] = files

    checksums = {
        key: value
        for key, value in dict(manifest.get("checksums", {})).items()
        if not key.startswith("builder/") and not key.startswith("generated/")
    }
    manifest["checksums"] = checksums

    keep_kinds = {"auxiliary-v4-promotion", "transport-table-refresh"}
    provenance = [item for item in manifest.get("provenance", []) if item.get("kind") in keep_kinds]
    manifest["provenance"] = provenance

    manifest.pop("generated", None)
    write_manifest(bundle_root, manifest)

    if delete_files:
        for dirname in ("builder", "generated"):
            target = bundle_root / dirname
            if target.exists():
                shutil.rmtree(target)
    return bundle_root
