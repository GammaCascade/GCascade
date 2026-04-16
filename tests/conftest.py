from __future__ import annotations

import os
from pathlib import Path

import pytest

import gcascade_v5 as g
import gcascade_v5.bundle as bundle
import gcascade_v5.legacy as legacy


REPO_ROOT = Path(__file__).resolve().parents[1]
LEGACY_LIBRARY = REPO_ROOT / "LibrariesV5-legacy"
LEGACY_GENERATED = REPO_ROOT / "generated_libraries"


def require_legacy_library() -> None:
    if not LEGACY_LIBRARY.exists():
        pytest.skip("Legacy GCascade tables are not available in this workspace.")


@pytest.fixture(scope="session")
def quick_bundle_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    require_legacy_library()
    target = tmp_path_factory.mktemp("gcascade_bundle_quick") / "bundle"
    return bundle.convert_legacy_library(
        LEGACY_LIBRARY,
        target,
        overwrite=False,
        ebl_indices=[1],
        include_builder=True,
        generated_source_path=LEGACY_GENERATED if LEGACY_GENERATED.exists() else None,
    )


@pytest.fixture(scope="session")
def slow_bundle_path(tmp_path_factory: pytest.TempPathFactory) -> Path:
    if os.getenv("GCASCADE_RUN_SLOW") != "1":
        pytest.skip("Set GCASCADE_RUN_SLOW=1 to build the slow regression bundle.")
    require_legacy_library()
    target = tmp_path_factory.mktemp("gcascade_bundle_slow") / "bundle"
    return bundle.convert_legacy_library(
        LEGACY_LIBRARY,
        target,
        overwrite=False,
        ebl_indices=[0, 1, 4],
        include_builder=True,
        generated_source_path=LEGACY_GENERATED if LEGACY_GENERATED.exists() else None,
    )


@pytest.fixture(autouse=True)
def reset_runtime_state():
    g.reset_state()
    legacy.reset_state()
    yield
    g.reset_state()
    legacy.reset_state()


def pytest_configure(config: pytest.Config) -> None:
    config.addinivalue_line("markers", "slow_regression: uses the legacy tables and slower cascade parity runs")


def skip_unless_slow_enabled() -> None:
    if os.getenv("GCASCADE_RUN_SLOW") != "1":
        pytest.skip("Set GCASCADE_RUN_SLOW=1 to run the slow legacy parity suite.")
