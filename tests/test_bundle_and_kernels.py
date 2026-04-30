from __future__ import annotations

from pathlib import Path

import h5py
import numpy as np
import pytest

import gcascade_v5 as g
import gcascade_v5.bundle as bundle
from gcascade_v5.redshift import redshift_cycle, reference_redshift_cycle


def test_runtime_bundle_contains_only_runtime_hdf5_and_json(quick_bundle_path: Path):
    suffixes = {path.suffix for path in quick_bundle_path.rglob("*") if path.is_file()}
    assert suffixes <= {"", ".h5", ".json"}


def test_bundle_info_reports_runtime_files(quick_bundle_path: Path):
    info = bundle.bundle_info(quick_bundle_path)
    assert info["available_ebl_indices"] == [0, 1, 2, 3, 4, 5, 6]
    assert info["schema_version"] == 2
    assert info["files"]["runtime_common"] == "runtime/common.h5"
    assert info["files"]["runtime"]["1"] == "runtime/ebl_SL.h5"


def test_redshift_kernel_matches_reference_on_every_runtime_window():
    spec = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    left_idx, weights, _ = bundle._build_redshift_tables()

    for window_idx in range(len(g.diffuseDistances)):
        z_hi = float(g.diffuseDistances[window_idx])
        z_lo = 0.0 if window_idx == 0 else float(g.diffuseDistances[window_idx - 1])
        reference_out = reference_redshift_cycle(spec, np.array([z_hi, z_lo], dtype=np.float64))
        new_out = redshift_cycle(spec, left_idx=left_idx[window_idx], weights=weights[window_idx])
        assert np.allclose(new_out, reference_out, rtol=1.0e-12, atol=0.0)


def test_bundle_point_cascade_returns_gamma_and_electron_state(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    ele = g.cutoffPowerLaw(g.energies, gamma=2.4, cutoff=1e6, amp=1e35)
    result = g.CascadePoint(inj, 1.0e-6, electronSpectraPre=ele, return_state=True)
    assert isinstance(result, g.CascadeResult)
    assert result.gamma.shape == g.energies.shape
    assert result.electron.shape == g.energies.shape
    assert np.all(result.gamma >= 0.0)
    assert np.all(result.electron >= 0.0)
    assert result.diagnostics["initial_electron_energy"] > 0.0
    init = result.diagnostics["initial_gamma_energy"] + result.diagnostics["initial_electron_energy"]
    assert abs(result.diagnostics["energy_residual"]) <= 1.0e-9 * init
    assert result.diagnostics["relative_energy_residual"] == pytest.approx(
        result.diagnostics["energy_residual"] / init
    )


def test_point_cascade_tracks_redshift_energy_loss(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    result = g.CascadePoint(inj, 0.3, return_state=True)
    assert result.diagnostics["redshift_energy_lost"] > 0.0


def test_point_cascade_conserves_energy_with_synchrotron_and_electrons(quick_bundle_path: Path):
    g.set_library_path(quick_bundle_path)
    g.changeMagneticField(1.0e-9, 0.0, 1)
    try:
        gamma = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
        electron = g.cutoffPowerLaw(g.energies, gamma=2.4, cutoff=1e6, amp=3.0e39)
        result = g.CascadePoint(gamma, 0.3, electronSpectraPre=electron, return_state=True)
        init = result.diagnostics["initial_gamma_energy"] + result.diagnostics["initial_electron_energy"]
        assert result.diagnostics["synchrotron_energy_lost"] > 0.0
        assert abs(result.diagnostics["energy_residual"]) <= 1.0e-6 * init
        assert abs(result.diagnostics["relative_energy_residual"]) <= 1.0e-6
    finally:
        g.changeMagneticField(0.0, 0.0, 1)


def test_runtime_bundle_contains_transport_tables(quick_bundle_path: Path):
    runtime_path = bundle.resolve_runtime_ebl_path(quick_bundle_path, 1)
    with h5py.File(runtime_path, "r") as handle:
        for name in [
            "imfp",
            "extinction_coeffs",
            "attenuation_vectors",
            "pp_packed",
            "ics_imfp",
            "ics_extinction_coeffs",
            "ics_gamma_packed",
            "ics_electron_packed",
            "dEdt_ics",
        ]:
            assert name in handle
        assert handle["ics_imfp"].shape == (len(g.zReg), len(g.energies))
        assert handle["ics_gamma_packed"].shape == (len(g.zReg), bundle.PACKED_TRIANGULAR_SIZE)
        assert handle["ics_electron_packed"].shape == (len(g.zReg), bundle.PACKED_TRIANGULAR_SIZE)
