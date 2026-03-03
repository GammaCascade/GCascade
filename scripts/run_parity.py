#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
import sys

import numpy as np

from gcascade_v5 import (
    RedshiftPoint,
    AttenuatePoint,
    CascadePoint,
    RedshiftDiffuse,
    AttenuateDiffuse,
    CascadeDiffuse,
    RedshiftEvolving,
    AttenuateEvolving,
    CascadeEvolving,
    changeEBLModel,
)


FUNCTIONS = {
    "RedshiftPoint": RedshiftPoint,
    "AttenuatePoint": AttenuatePoint,
    "CascadePoint": CascadePoint,
    "RedshiftDiffuse": RedshiftDiffuse,
    "AttenuateDiffuse": AttenuateDiffuse,
    "CascadeDiffuse": CascadeDiffuse,
    "RedshiftEvolving": RedshiftEvolving,
    "AttenuateEvolving": AttenuateEvolving,
    "CascadeEvolving": CascadeEvolving,
}


def load_csv(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", dtype=np.float64)


def compute_error(actual: np.ndarray, expected: np.ndarray, abs_floor: float) -> np.ndarray:
    denom = np.maximum(np.abs(expected), abs_floor)
    return np.abs(actual - expected) / denom


def run_case(case_dir: Path, rel_tol: float, abs_floor: float) -> tuple[bool, str]:
    meta_path = case_dir / "meta.json"
    if not meta_path.exists():
        return False, f"{case_dir.name}: missing meta.json"

    with meta_path.open("r", encoding="utf-8") as handle:
        meta = json.load(handle)

    fn_name = meta.get("function")
    if fn_name not in FUNCTIONS:
        return False, f"{case_dir.name}: unknown function '{fn_name}'"

    fn = FUNCTIONS[fn_name]

    ebl_index = meta.get("ebl_index")
    if ebl_index is not None:
        try:
            changeEBLModel(int(ebl_index))
        except ValueError:
            # already in this model
            pass

    z_start = float(meta["z_start"])

    if fn_name.endswith("Point"):
        inj = load_csv(case_dir / "inj.csv")
        args = [inj, z_start]
    elif fn_name.endswith("Diffuse"):
        inj = load_csv(case_dir / "inj.csv")
        z_distrib = load_csv(case_dir / "z_distrib.csv")
        args = [inj, z_start, z_distrib]
    elif fn_name.endswith("Evolving"):
        inj2d = load_csv(case_dir / "inj2d.csv")
        z_distrib = load_csv(case_dir / "z_distrib.csv")
        args = [inj2d, z_start, z_distrib]
    else:
        return False, f"{case_dir.name}: unhandled function category"

    expected = load_csv(case_dir / "expected.csv")
    actual = np.asarray(fn(*args), dtype=np.float64)

    if actual.shape != expected.shape:
        return False, f"{case_dir.name}: shape mismatch {actual.shape} != {expected.shape}"

    rel_err = compute_error(actual, expected, abs_floor=abs_floor)
    max_err = float(np.max(rel_err))

    (case_dir / "actual.csv").write_text(
        "\n".join(
            ",".join(f"{v:.18e}" for v in row)
            for row in np.atleast_2d(actual)
        )
        + "\n",
        encoding="utf-8",
    )

    (case_dir / "diff.csv").write_text(
        "\n".join(
            ",".join(f"{v:.18e}" for v in row)
            for row in np.atleast_2d(rel_err)
        )
        + "\n",
        encoding="utf-8",
    )

    ok = max_err <= rel_tol
    msg = (
        f"{case_dir.name}: {'PASS' if ok else 'FAIL'} | "
        f"function={fn_name} | max_rel_err={max_err:.3e} | tol={rel_tol:.3e}"
    )
    return ok, msg


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    fixtures_root = root / "benchmarks" / "v4_reference"

    rel_tol = float(sys.argv[1]) if len(sys.argv) > 1 else 1.0e-3
    abs_floor = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0e-45

    case_dirs = sorted([d for d in fixtures_root.iterdir() if d.is_dir()])
    if not case_dirs:
        print("No parity fixtures found. Add cases under benchmarks/v4_reference/<case_name>/")
        return 0

    any_fail = False
    for case_dir in case_dirs:
        ok, msg = run_case(case_dir, rel_tol=rel_tol, abs_floor=abs_floor)
        print(msg)
        if not ok:
            any_fail = True

    return 1 if any_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
