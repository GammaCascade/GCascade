import numpy as np
import pytest

import gcascade_v5 as g


def test_point_rejects_bad_injected_shape():
    with pytest.raises(ValueError):
        g.RedshiftPoint(np.zeros(10), 0.1)


def test_point_rejects_bad_electron_shape(quick_bundle_path):
    g.set_library_path(quick_bundle_path)
    with pytest.raises(ValueError):
        g.CascadePoint(np.zeros_like(g.energies), 0.1, electronSpectraPre=np.zeros(10))


def test_diffuse_rejects_bad_distribution_shape():
    inj = np.zeros_like(g.energies)
    with pytest.raises(ValueError):
        g.RedshiftDiffuse(inj, 0.1, np.zeros(10))


def test_evolving_rejects_bad_injected_shape():
    bad = np.zeros((100, 100))
    z = np.zeros_like(g.diffuseDistances)
    with pytest.raises(ValueError):
        g.RedshiftEvolving(bad, 0.5, z)


def test_evolving_rejects_bad_electron_shape(quick_bundle_path):
    g.set_library_path(quick_bundle_path)
    good = np.zeros((len(g.diffuseDistances), len(g.energies)))
    z = np.zeros_like(g.diffuseDistances)
    with pytest.raises(ValueError):
        g.CascadeEvolving(good, 0.5, z, electronSpectra=np.zeros((10, 10)))


def test_zstart_validation():
    inj = np.zeros_like(g.energies)
    with pytest.raises(ValueError):
        g.RedshiftPoint(inj, 11.0)
    with pytest.raises(ValueError):
        g.RedshiftPoint(inj, -1.0)
