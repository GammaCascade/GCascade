import numpy as np

import gcascade_v5 as g


if __name__ == "__main__":
    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    out = g.RedshiftPoint(inj, 0.3)

    print("Input shape:", inj.shape)
    print("Output shape:", out.shape)
    print("Min/Max output:", np.min(out), np.max(out))
