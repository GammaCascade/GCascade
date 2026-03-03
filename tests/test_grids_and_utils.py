import numpy as np
import gcascade_v5 as g


def test_grids_have_expected_shapes():
    assert g.energies.shape == (300,)
    assert g.diffuseDistances.shape == (1036,)
    assert g.diffuseSteps.shape == (1036,)
    assert g.zReg.shape == (1001,)


def test_grids_monotonic():
    assert np.all(np.diff(g.energies) > 0)
    assert np.all(np.diff(g.diffuseDistances) > 0)
    assert np.all(np.diff(g.zReg) > 0)


def test_cutoff_power_law_shape_and_positive():
    spec = g.cutoffPowerLaw(g.energies, gamma=2.0, cutoff=1e6, amp=1e40)
    assert spec.shape == g.energies.shape
    assert np.all(spec >= 0.0)


def test_hubble_scalar_and_array():
    assert g.hubble(0.0) > 0.0
    arr = g.hubble(np.array([0.0, 0.5, 1.0]))
    assert arr.shape == (3,)
    assert np.all(arr > 0.0)
