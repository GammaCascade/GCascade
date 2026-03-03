import numpy as np
import gcascade_v5 as g


def test_redshift_point_is_deterministic():
    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    out1 = g.RedshiftPoint(inj, 0.3)
    out2 = g.RedshiftPoint(inj, 0.3)
    assert np.array_equal(out1, out2)
