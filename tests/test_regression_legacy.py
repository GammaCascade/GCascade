from __future__ import annotations

from pathlib import Path
import os

import numpy as np
import pytest

import gcascade_v5 as g
import gcascade_v5.legacy as legacy


@pytest.mark.slow_regression
@pytest.mark.parametrize("ebl_index", [0, 1, 4])
@pytest.mark.parametrize("z_start", [0.1, 0.3, 1.0, 3.0])
@pytest.mark.parametrize("mode", ["point", "diffuse", "evolving"])
def test_cascade_regression_against_legacy(
    slow_bundle_path: Path,
    ebl_index: int,
    z_start: float,
    mode: str,
):
    if os.getenv("GCASCADE_RUN_SLOW") != "1":
        pytest.skip("Set GCASCADE_RUN_SLOW=1 to run the slow legacy parity suite.")

    g.set_library_path(slow_bundle_path)
    legacy.set_library_path(Path("LibrariesV5-legacy").resolve())
    if ebl_index != 1:
        g.changeEBLModel(ebl_index)
        legacy.changeEBLModel(ebl_index)

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    z_distrib = np.exp(-g.diffuseDistances / 2.0)
    inj2d = np.tile(inj, (len(g.diffuseDistances), 1))

    if mode == "point":
        new_out = g.CascadePoint(inj, z_start)
        legacy_out = legacy.CascadePoint(inj, z_start)
    elif mode == "diffuse":
        new_out = g.CascadeDiffuse(inj, z_start, z_distrib)
        legacy_out = legacy.CascadeDiffuse(inj, z_start, z_distrib)
    else:
        new_out = g.CascadeEvolving(inj2d, z_start, z_distrib)
        legacy_out = legacy.CascadeEvolving(inj2d, z_start, z_distrib)

    assert np.allclose(new_out, legacy_out, rtol=1.0e-7, atol=0.0)
    assert np.isclose(np.sum(new_out), np.sum(legacy_out), rtol=1.0e-5, atol=0.0)


@pytest.mark.slow_regression
def test_change_magnetic_field_creates_and_activates_hdf5_variant(quick_bundle_path: Path):
    if os.getenv("GCASCADE_RUN_SLOW") != "1":
        pytest.skip("Set GCASCADE_RUN_SLOW=1 to run the slow legacy parity suite.")

    g.set_library_path(quick_bundle_path)
    before = g.list_generated_variants(1)
    g.changeMagneticField(1.0e-18, 0.0, 1)
    after = g.list_generated_variants(1)

    assert len(after) == len(before) + 1
    active = g.get_bundle_info()["active_generated_variants"]["1"]
    assert active["path"].endswith(".h5")
    assert active["variant_id"] == after[-1]["variant_id"]
