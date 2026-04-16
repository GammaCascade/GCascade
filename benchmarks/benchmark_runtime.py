from __future__ import annotations

import argparse
import json
from pathlib import Path
import tempfile
import time

import numpy as np

import gcascade_v5 as g
import gcascade_v5.bundle as bundle
import gcascade_v5.legacy as legacy


def benchmark_once(func, *args):
    start = time.perf_counter()
    output = func(*args)
    elapsed = time.perf_counter() - start
    return elapsed, float(np.sum(output))


def main() -> None:
    parser = argparse.ArgumentParser(description="Benchmark GCascadeV5 legacy vs bundle runtime.")
    parser.add_argument("--legacy-library", type=Path, default=Path("LibrariesV5-legacy"))
    parser.add_argument("--generated-library", type=Path, default=Path("generated_libraries"))
    parser.add_argument("--bundle-root", type=Path, default=None)
    parser.add_argument("--z", type=float, default=0.3)
    args = parser.parse_args()

    legacy_library = args.legacy_library.expanduser().resolve()
    generated_library = args.generated_library.expanduser().resolve()
    bundle_root = args.bundle_root
    if bundle_root is None:
        bundle_root = Path(tempfile.mkdtemp(prefix="gcascade-benchmark-bundle-", dir=".")) / "bundle"
        bundle.convert_legacy_library(
            legacy_library,
            bundle_root,
            overwrite=False,
            ebl_indices=[1],
            include_builder=True,
            generated_source_path=generated_library if generated_library.exists() else None,
        )
    bundle_root = bundle_root.expanduser().resolve()

    inj = g.cutoffPowerLaw(g.energies, gamma=2.2, cutoff=1e7, amp=1e40)

    legacy.reset_state()
    legacy.set_library_path(legacy_library)
    g.reset_state()
    g.set_library_path(bundle_root)

    result = {"z": args.z, "bundle_root": str(bundle_root), "benchmarks": {}}
    for name in ["AttenuatePoint", "CascadePoint"]:
        legacy_time, legacy_sum = benchmark_once(getattr(legacy, name), inj, args.z)
        new_time, new_sum = benchmark_once(getattr(g, name), inj, args.z)
        result["benchmarks"][name] = {
            "legacy_seconds": legacy_time,
            "bundle_seconds": new_time,
            "speedup": legacy_time / new_time if new_time > 0 else None,
            "legacy_sum": legacy_sum,
            "bundle_sum": new_sum,
        }

    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
