from __future__ import annotations

from pathlib import Path

from . import bundle


def bundle_exists(path: str | Path) -> bool:
    return bundle.is_bundle_root(path)


def manifest_path(path: str | Path) -> Path:
    return bundle.bundle_paths(path).manifest
