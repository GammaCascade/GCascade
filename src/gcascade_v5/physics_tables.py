from __future__ import annotations

"""Helpers for checking that a GCascade bundle contains the transport tables it needs."""

from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

from .physics import EBL_NAME_MAP


REQUIRED_RUNTIME_DATASETS = (
    "imfp",
    "extinction_coeffs",
    "attenuation_vectors",
    "pp_packed",
    "ics_imfp",
    "ics_extinction_coeffs",
    "ics_gamma_packed",
    "ics_electron_packed",
    "dEdt_ics",
)


@dataclass(frozen=True)
class TransportSourcePaths:
    """Filesystem layout of the runtime bundle currently used by the code."""

    root: Path
    manifest: Path
    runtime_dir: Path


def transport_source_paths(root: str | Path) -> TransportSourcePaths:
    """Return the standard locations of the manifest and runtime files inside a bundle."""
    source_root = Path(root).expanduser().resolve()
    return TransportSourcePaths(
        root=source_root,
        manifest=source_root / "bundle_manifest.json",
        runtime_dir=source_root / "runtime",
    )


def required_transport_files(root: str | Path, ebl_index: int) -> list[Path]:
    """List the manifest and runtime file required for a given EBL model."""
    source = transport_source_paths(root)
    ebl_name = EBL_NAME_MAP[int(ebl_index)]
    return [
        source.manifest,
        source.runtime_dir / f"ebl_{ebl_name}.h5",
    ]


def validate_transport_source_root(
    root: str | Path,
    ebl_indices: list[int] | tuple[int, ...] | None = None,
) -> Path:
    """Check that a bundle contains the runtime HDF5 tables needed by the cascade."""
    source_root = Path(root).expanduser().resolve()
    if not source_root.is_dir():
        raise FileNotFoundError(f"GCascade bundle not found at {source_root}.")

    indices = [int(idx) for idx in (ebl_indices if ebl_indices is not None else sorted(EBL_NAME_MAP))]
    missing: list[str] = []
    for idx in indices:
        for path in required_transport_files(source_root, idx):
            if not path.exists():
                missing.append(str(path))
        runtime_path = transport_source_paths(source_root).runtime_dir / f"ebl_{EBL_NAME_MAP[idx]}.h5"
        if not runtime_path.exists():
            continue
        with h5py.File(runtime_path, "r") as handle:
            for dataset_name in REQUIRED_RUNTIME_DATASETS:
                if dataset_name not in handle:
                    missing.append(f"{runtime_path}:{dataset_name}")
    if missing:
        preview = "\n".join(f"- {path}" for path in missing[:12])
        extra = "" if len(missing) <= 12 else f"\n... and {len(missing) - 12} more"
        raise FileNotFoundError(
            "The selected GCascade bundle is missing transport tables required for gamma/electron cascades:\n"
            f"{preview}{extra}"
        )
    return source_root


def load_runtime_dataset(root: str | Path, ebl_index: int, dataset_name: str) -> np.ndarray:
    """Load one runtime dataset from the chosen EBL file as a NumPy array."""
    source_root = validate_transport_source_root(root, [int(ebl_index)])
    runtime_path = transport_source_paths(source_root).runtime_dir / f"ebl_{EBL_NAME_MAP[int(ebl_index)]}.h5"
    with h5py.File(runtime_path, "r") as handle:
        return np.asarray(handle[dataset_name], dtype=np.float64)

