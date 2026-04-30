from pathlib import Path

import pytest

import gcascade_v5.physics_tables as physics_tables


def test_transport_source_runtime_tables_are_detected():
    root = Path("LibrariesV5").resolve()
    if not root.exists():
        pytest.skip("LibrariesV5 transport tables are not available in this workspace.")
    physics_tables.validate_transport_source_root(root, [1])
    required = physics_tables.required_transport_files(root, 1)
    assert required
    assert all(path.exists() for path in required)


def test_transport_source_ics_imfp_shape():
    root = Path("LibrariesV5").resolve()
    if not root.exists():
        pytest.skip("LibrariesV5 transport tables are not available in this workspace.")
    imfp = physics_tables.load_runtime_dataset(root, 1, "ics_imfp")
    assert imfp.shape == (1001, 300)
    assert imfp.min() >= 0.0
