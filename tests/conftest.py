from __future__ import annotations

from pathlib import Path

import pytest

import gcascade_v5 as g


REPO_ROOT = Path(__file__).resolve().parents[1]
RUNTIME_LIBRARY = REPO_ROOT / "LibrariesV5"


def require_runtime_library() -> None:
    """Skip bundle-dependent tests when the runtime HDF5 library is unavailable."""
    if not RUNTIME_LIBRARY.exists():
        pytest.skip("LibrariesV5 is not available in this workspace.")


@pytest.fixture(scope="session")
def quick_bundle_path() -> Path:
    """Return the checked-in runtime bundle used by the modern cascade tests."""
    require_runtime_library()
    return RUNTIME_LIBRARY.resolve()


@pytest.fixture(autouse=True)
def reset_runtime_state():
    """Reset the singleton runtime state around every test for isolation."""
    g.reset_state()
    yield
    g.reset_state()

