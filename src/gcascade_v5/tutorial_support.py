from __future__ import annotations

"""Small helpers used by the tutorial notebook when locating the runtime bundle."""

from pathlib import Path

from . import bundle


def bundle_exists(path: str | Path) -> bool:
    """Check whether a path points to a valid GCascadeV5 runtime bundle."""
    return bundle.is_bundle_root(path)


def manifest_path(path: str | Path) -> Path:
    """Return the manifest file path for a given bundle directory."""
    return bundle.bundle_paths(path).manifest

