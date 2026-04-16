from __future__ import annotations

from pathlib import Path
import shutil

import numpy as np

import gcascade_v5 as g
import gcascade_v5.bundle as bundle
import gcascade_v5.legacy as legacy
from gcascade_v5.redshift import redshift_cycle


def test_converted_bundle_contains_only_hdf5_and_json(quick_bundle_path: Path):
    suffixes = {path.suffix for path in quick_bundle_path.rglob("*") if path.is_file()}
    assert suffixes <= {".h5", ".json"}


def test_bundle_manifest_tracks_imported_generated_variant(quick_bundle_path: Path):
    info = bundle.bundle_info(quick_bundle_path)
    assert info["available_ebl_indices"] == [1]

    variants = bundle.generated_variants(quick_bundle_path, 1)
    assert variants
    assert variants[0]["variant_id"].startswith("imported_legacy_")
    active = info["active_generated_variants"]["1"]
    assert active["path"].endswith(".h5")


def test_active_cycle_path_matches_manifest(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)
    active = g.get_active_cycle_path(1)
    assert active == quick_bundle_path / "generated" / "imported_legacy_SL.h5"


def test_set_active_cycle_path_and_factory_reset(quick_bundle_path: Path, tmp_path: Path):
    local_bundle = tmp_path / "bundle"
    shutil.copytree(quick_bundle_path, local_bundle)
    g.set_library_path(local_bundle)

    generated_path = local_bundle / "generated" / "imported_legacy_SL.h5"
    selected = g.set_active_cycle_path(generated_path)
    assert selected == generated_path
    assert g.get_active_cycle_path() == generated_path

    g.reset_factory_settings()
    assert g.EBLindex == 1
    assert g.get_active_cycle_path() == local_bundle / "runtime" / "ebl_SL.h5"


def test_redshift_kernel_matches_legacy_on_every_runtime_window():
    spec = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    left_idx, weights, _ = bundle._build_redshift_tables()

    for window_idx in range(len(g.diffuseDistances)):
        z_hi = float(g.diffuseDistances[window_idx])
        z_lo = 0.0 if window_idx == 0 else float(g.diffuseDistances[window_idx - 1])
        legacy_out = legacy.RedshiftingCycle(spec, np.array([z_hi, z_lo], dtype=np.float64))
        new_out = redshift_cycle(spec, left_idx=left_idx[window_idx], weights=weights[window_idx], use_numba=False)
        assert np.allclose(new_out, legacy_out, rtol=1.0e-12, atol=0.0)


def test_bundle_point_attenuation_matches_legacy(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)
    legacy.set_library_path(Path("LibrariesV5-legacy").resolve())

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    for z_start in [0.1, 0.3]:
        new_out = g.AttenuatePoint(inj, z_start)
        legacy_out = legacy.AttenuatePoint(inj, z_start)
        assert np.allclose(new_out, legacy_out, rtol=1.0e-10, atol=0.0)


def test_bundle_point_cascade_matches_legacy_smoke(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)
    legacy.set_library_path(Path("LibrariesV5-legacy").resolve())

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    for z_start in [0.1, 0.3]:
        new_out = g.CascadePoint(inj, z_start)
        legacy_out = legacy.CascadePoint(inj, z_start)
        assert np.allclose(new_out, legacy_out, rtol=1.0e-9, atol=0.0)
