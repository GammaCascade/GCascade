from __future__ import annotations

"""Simple timing helper for the current GCascadeV5 runtime bundle."""

import argparse
import json
from pathlib import Path
import time

import numpy as np

import gcascade_v5 as g


def benchmark_once(func, *args, **kwargs):
    """Run one calculation once and return its wall time plus a simple spectral checksum."""
    start = time.perf_counter()
    output = func(*args, **kwargs)
    elapsed = time.perf_counter() - start
    gamma = output.gamma if hasattr(output, "gamma") else output
    return elapsed, float(np.sum(gamma))


def main() -> None:
    """Benchmark the point-source attenuation and cascade calls for one bundle."""
    parser = argparse.ArgumentParser(description="Benchmark the GCascadeV5 runtime bundle.")
    parser.add_argument("--bundle-root", type=Path, default=Path("LibrariesV5"))
    parser.add_argument("--z", type=float, default=0.3)
    parser.add_argument("--with-electrons", action="store_true")
    args = parser.parse_args()

    bundle_root = args.bundle_root.expanduser().resolve()
    g.reset_state()
    g.set_library_path(bundle_root)

    gamma_inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)
    electron_inj = None
    if args.with_electrons:
        electron_inj = g.cutoffPowerLaw(g.energies, gamma=2.4, cutoff=1e6, amp=3e39)

    result = {
        "z": args.z,
        "bundle_root": str(bundle_root),
        "with_electrons": bool(args.with_electrons),
        "benchmarks": {},
    }

    atten_time, atten_sum = benchmark_once(g.AttenuatePoint, gamma_inj, args.z)
    result["benchmarks"]["AttenuatePoint"] = {
        "seconds": atten_time,
        "gamma_sum": atten_sum,
    }

    cascade_time, cascade_sum = benchmark_once(
        g.CascadePoint,
        gamma_inj,
        args.z,
        electronSpectraPre=electron_inj,
        return_state=bool(args.with_electrons),
    )
    result["benchmarks"]["CascadePoint"] = {
        "seconds": cascade_time,
        "gamma_sum": cascade_sum,
    }

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

