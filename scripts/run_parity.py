#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
from typing import Any

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
    changeMagneticField,
    set_generated_library_path,
    set_numba,
    set_progress,
    reset_state,
    is_numba_available,
)
import gcascade_v5.core as gc_core


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

REQUIRED_META_KEYS = {
    "function",
    "z_start",
    "inputs_kind",
    "pre_actions",
    "parity_targets",
    "v4_ref",
}


def load_csv(path: Path) -> np.ndarray:
    return np.loadtxt(path, delimiter=",", dtype=np.float64)


def stable_rel_error(actual: np.ndarray, expected: np.ndarray, abs_floor: float) -> np.ndarray:
    denom = np.maximum(np.abs(expected), abs_floor)
    return np.abs(actual - expected) / denom


def validate_meta(meta: dict[str, Any]) -> None:
    missing = sorted(REQUIRED_META_KEYS - set(meta.keys()))
    if missing:
        raise ValueError(f"Missing required meta keys: {missing}")

    fn_name = str(meta["function"])
    if fn_name not in FUNCTIONS:
        raise ValueError(f"Unsupported function in meta: {fn_name}")

    inputs_kind = str(meta["inputs_kind"])
    if inputs_kind not in {"point", "diffuse", "evolving"}:
        raise ValueError(f"Unsupported inputs_kind in meta: {inputs_kind}")

    if not isinstance(meta["pre_actions"], list):
        raise ValueError("meta.pre_actions must be a list")

    targets = meta["parity_targets"]
    if not isinstance(targets, list) or not targets:
        raise ValueError("meta.parity_targets must be a non-empty list")

    for target in targets:
        if target not in {"output", "cycle_table_sparse"}:
            raise ValueError(f"Unsupported parity target: {target}")


def apply_pre_actions(pre_actions: list[dict[str, Any]]) -> None:
    for idx, action in enumerate(pre_actions, start=1):
        if not isinstance(action, dict):
            raise ValueError(f"pre_actions[{idx}] must be an object")
        action_name = action.get("action")
        action_args = action.get("args", [])

        if action_name == "changeEBLModel":
            if len(action_args) != 1:
                raise ValueError("changeEBLModel pre_action requires 1 arg")
            changeEBLModel(int(action_args[0]))
            continue

        if action_name == "changeMagneticField":
            if len(action_args) != 3:
                raise ValueError("changeMagneticField pre_action requires 3 args")
            changeMagneticField(float(action_args[0]), float(action_args[1]), int(action_args[2]))
            continue

        raise ValueError(f"Unknown pre_action: {action_name}")


def load_case_inputs(case_dir: Path, inputs_kind: str) -> tuple[list[Any], float]:
    meta_path = case_dir / "meta.json"
    with meta_path.open("r", encoding="utf-8") as handle:
        meta = json.load(handle)

    z_start = float(meta["z_start"])

    if inputs_kind == "point":
        inj = load_csv(case_dir / "inj.csv")
        return [inj, z_start], z_start

    if inputs_kind == "diffuse":
        inj = load_csv(case_dir / "inj.csv")
        z_distrib = load_csv(case_dir / "z_distrib.csv")
        return [inj, z_start, z_distrib], z_start

    if inputs_kind == "evolving":
        inj2d = load_csv(case_dir / "inj2d.csv")
        z_distrib = load_csv(case_dir / "z_distrib.csv")
        return [inj2d, z_start, z_distrib], z_start

    raise ValueError(f"Unsupported inputs_kind: {inputs_kind}")


def compare_with_thresholds(
    actual: np.ndarray,
    expected: np.ndarray,
    *,
    rel_tol: float,
    abs_floor: float,
    near_zero_threshold: float,
) -> tuple[bool, float, float, np.ndarray]:
    if actual.shape != expected.shape:
        raise ValueError(f"shape mismatch {actual.shape} != {expected.shape}")

    abs_err = np.abs(actual - expected)
    mask = np.abs(expected) >= near_zero_threshold

    rel_err = np.zeros_like(actual, dtype=np.float64)
    if np.any(mask):
        rel_err[mask] = abs_err[mask] / np.abs(expected[mask])

    max_rel = float(np.max(rel_err[mask])) if np.any(mask) else 0.0
    max_abs_near_zero = float(np.max(abs_err[~mask])) if np.any(~mask) else 0.0

    pass_rel = bool(np.all(rel_err[mask] <= rel_tol)) if np.any(mask) else True
    pass_abs = bool(np.all(abs_err[~mask] <= abs_floor)) if np.any(~mask) else True

    # Stable relative error view for diff artifact files.
    diff_artifact = stable_rel_error(actual, expected, abs_floor=abs_floor)

    return pass_rel and pass_abs, max_rel, max_abs_near_zero, diff_artifact


def write_output_artifacts(case_dir: Path, actual: np.ndarray, diff: np.ndarray) -> None:
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
            for row in np.atleast_2d(diff)
        )
        + "\n",
        encoding="utf-8",
    )


def load_cycle_sparse(path: Path) -> np.ndarray:
    arr = load_csv(path)
    arr = np.atleast_2d(arr)
    if arr.shape[1] != 4:
        raise ValueError(f"expected cycle sparse CSV with 4 columns, got shape {arr.shape}")
    return arr


def sample_cycle_sparse(indices: dict[str, list[int]]) -> np.ndarray:
    state = gc_core._state()
    state.ensure_ebl_loaded()
    cycle = state.cycle_spec

    rows: list[list[float]] = []
    for z_idx in indices["z_idx"]:
        for e_in_idx in indices["e_in_idx"]:
            for e_out_idx in indices["e_out_idx"]:
                rows.append([
                    float(z_idx),
                    float(e_in_idx),
                    float(e_out_idx),
                    float(cycle[int(z_idx), int(e_in_idx), int(e_out_idx)]),
                ])

    return np.asarray(rows, dtype=np.float64)


def prepare_case_runtime(case_generated_dir: Path, numba_mode: str) -> None:
    reset_state()
    set_progress(False)
    set_generated_library_path(case_generated_dir)

    if numba_mode == "off":
        set_numba(False)
    elif numba_mode == "on":
        set_numba(True)


def run_case(
    case_dir: Path,
    *,
    output_rel_tol: float,
    output_abs_floor: float,
    cycle_rel_tol: float,
    cycle_abs_floor: float,
    near_zero_threshold: float,
    determinism_check: bool,
    runtime_generated_root: Path,
    numba_mode: str,
) -> dict[str, Any]:
    meta_path = case_dir / "meta.json"
    with meta_path.open("r", encoding="utf-8") as handle:
        meta = json.load(handle)

    validate_meta(meta)

    fn_name = str(meta["function"])
    fn = FUNCTIONS[fn_name]
    inputs_kind = str(meta["inputs_kind"])
    parity_targets = list(meta["parity_targets"])
    pre_actions = list(meta["pre_actions"])

    case_generated_dir = runtime_generated_root / case_dir.name
    case_generated_dir.mkdir(parents=True, exist_ok=True)
    prepare_case_runtime(case_generated_dir, numba_mode=numba_mode)

    if pre_actions:
        apply_pre_actions(pre_actions)
    elif "ebl_index" in meta and int(meta["ebl_index"]) != 1:
        changeEBLModel(int(meta["ebl_index"]))

    args, _ = load_case_inputs(case_dir, inputs_kind)

    actual_output = np.asarray(fn(*args), dtype=np.float64)
    expected_output = load_csv(case_dir / "expected.csv")

    output_ok, output_max_rel, output_max_abs_near, output_diff = compare_with_thresholds(
        actual_output,
        expected_output,
        rel_tol=output_rel_tol,
        abs_floor=output_abs_floor,
        near_zero_threshold=near_zero_threshold,
    )
    write_output_artifacts(case_dir, actual_output, output_diff)

    determinism_ok = True
    if determinism_check:
        actual_repeat = np.asarray(fn(*args), dtype=np.float64)
        determinism_ok = bool(np.array_equal(actual_output, actual_repeat))

    cycle_ok = True
    cycle_max_rel = 0.0
    cycle_max_abs_near = 0.0

    if "cycle_table_sparse" in parity_targets:
        if "cycle_sparse_indices" not in meta:
            raise ValueError(f"cycle_table_sparse target requires cycle_sparse_indices in {meta_path}")

        expected_cycle = load_cycle_sparse(case_dir / "expected_cycle_sparse.csv")
        actual_cycle = sample_cycle_sparse(meta["cycle_sparse_indices"])

        if actual_cycle.shape != expected_cycle.shape:
            raise ValueError(f"cycle sparse shape mismatch {actual_cycle.shape} != {expected_cycle.shape}")

        # Compare only value column; index columns should match exactly by construction.
        cycle_ok, cycle_max_rel, cycle_max_abs_near, cycle_diff = compare_with_thresholds(
            actual_cycle[:, 3],
            expected_cycle[:, 3],
            rel_tol=cycle_rel_tol,
            abs_floor=cycle_abs_floor,
            near_zero_threshold=near_zero_threshold,
        )

        np.savetxt(case_dir / "actual_cycle_sparse.csv", actual_cycle, delimiter=",")
        np.savetxt(
            case_dir / "diff_cycle_sparse.csv",
            np.column_stack([expected_cycle[:, :3], cycle_diff]),
            delimiter=",",
        )

    passed = bool(output_ok and cycle_ok and determinism_ok)

    return {
        "case": case_dir.name,
        "function": fn_name,
        "inputs_kind": inputs_kind,
        "passed": passed,
        "output_ok": bool(output_ok),
        "cycle_ok": bool(cycle_ok),
        "determinism_ok": bool(determinism_ok),
        "output_max_rel": float(output_max_rel),
        "output_max_abs_near_zero": float(output_max_abs_near),
        "cycle_max_rel": float(cycle_max_rel),
        "cycle_max_abs_near_zero": float(cycle_max_abs_near),
        "parity_targets": parity_targets,
    }


def write_summary_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def write_summary_markdown(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []
    lines.append("# GCascadeV5 Parity Summary")
    lines.append("")
    lines.append(f"- Total cases: {payload['total_cases']}")
    lines.append(f"- Passed: {payload['passed_cases']}")
    lines.append(f"- Failed: {payload['failed_cases']}")
    lines.append(f"- Output rel tol: {payload['tolerances']['output_rel_tol']:.3e}")
    lines.append(f"- Cycle rel tol: {payload['tolerances']['cycle_rel_tol']:.3e}")
    lines.append(f"- Near-zero threshold: {payload['tolerances']['near_zero_threshold']:.3e}")
    lines.append(f"- Numba mode: {payload['numba_mode']}")
    lines.append(f"- Numba available: {payload['numba_available']}")
    lines.append("")

    worst = payload.get("worst_cases", [])
    if worst:
        lines.append("## Worst Cases")
        lines.append("")
        for item in worst:
            lines.append(
                f"- {item['case']} | pass={item['passed']} | "
                f"output_max_rel={item['output_max_rel']:.3e} | "
                f"cycle_max_rel={item['cycle_max_rel']:.3e}"
            )
        lines.append("")

    if payload.get("failed_case_names"):
        lines.append("## Failed Cases")
        lines.append("")
        for name in payload["failed_case_names"]:
            lines.append(f"- {name}")
        lines.append("")

    path.write_text("\n".join(lines), encoding="utf-8")


def maybe_load_matrix_tolerances(matrix_path: Path) -> dict[str, float]:
    if not matrix_path.exists():
        return {}
    with matrix_path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return dict(payload.get("tolerances", {}))


def main() -> int:
    root = Path(__file__).resolve().parents[1]
    default_fixtures_root = root / "benchmarks" / "v4_reference"
    default_matrix = root / "scripts" / "parity_matrix.json"
    default_summary_json = root / "benchmarks" / "parity_reports" / "latest_summary.json"
    default_summary_md = root / "benchmarks" / "parity_reports" / "latest_summary.md"

    parser = argparse.ArgumentParser(description="Run GCascadeV5 parity checks against V4 fixtures")
    parser.add_argument("--fixtures-root", type=Path, default=default_fixtures_root)
    parser.add_argument("--matrix", type=Path, default=default_matrix)
    parser.add_argument("--summary-json", type=Path, default=default_summary_json)
    parser.add_argument("--summary-md", type=Path, default=default_summary_md)
    parser.add_argument("--output-rel-tol", type=float, default=None)
    parser.add_argument("--output-abs-floor", type=float, default=None)
    parser.add_argument("--cycle-rel-tol", type=float, default=None)
    parser.add_argument("--cycle-abs-floor", type=float, default=None)
    parser.add_argument("--near-zero-threshold", type=float, default=None)
    parser.add_argument("--numba", choices=["off", "on", "auto"], default="off")
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--case-filter", type=str, default=None)
    parser.add_argument("--determinism", action="store_true")
    args = parser.parse_args()

    matrix_tols = maybe_load_matrix_tolerances(args.matrix)
    output_rel_tol = float(args.output_rel_tol if args.output_rel_tol is not None else matrix_tols.get("output_rel_tol", 1.0e-3))
    output_abs_floor = float(args.output_abs_floor if args.output_abs_floor is not None else matrix_tols.get("output_abs_floor", 1.0e-45))
    cycle_rel_tol = float(args.cycle_rel_tol if args.cycle_rel_tol is not None else matrix_tols.get("cycle_rel_tol", 5.0e-3))
    cycle_abs_floor = float(args.cycle_abs_floor if args.cycle_abs_floor is not None else matrix_tols.get("cycle_abs_floor", 1.0e-40))
    near_zero_threshold = float(args.near_zero_threshold if args.near_zero_threshold is not None else matrix_tols.get("near_zero_threshold", 1.0e-35))

    case_dirs = sorted([d for d in args.fixtures_root.iterdir() if d.is_dir() and (d / "meta.json").exists()])
    if args.case_filter:
        case_dirs = [d for d in case_dirs if args.case_filter in d.name]
    if args.limit is not None:
        case_dirs = case_dirs[: args.limit]

    if not case_dirs:
        print("No parity fixtures found. Add cases under benchmarks/v4_reference/<case_name>/")
        return 0

    if args.numba == "on" and not is_numba_available():
        print("Requested --numba on, but numba is unavailable in this environment.", file=sys.stderr)
        return 2

    runtime_generated_root = args.fixtures_root / "_runtime_generated"
    runtime_generated_root.mkdir(parents=True, exist_ok=True)

    results: list[dict[str, Any]] = []
    any_fail = False

    for case_dir in case_dirs:
        try:
            detail = run_case(
                case_dir,
                output_rel_tol=output_rel_tol,
                output_abs_floor=output_abs_floor,
                cycle_rel_tol=cycle_rel_tol,
                cycle_abs_floor=cycle_abs_floor,
                near_zero_threshold=near_zero_threshold,
                determinism_check=args.determinism,
                runtime_generated_root=runtime_generated_root,
                numba_mode=args.numba,
            )
        except Exception as exc:
            detail = {
                "case": case_dir.name,
                "function": "<error>",
                "inputs_kind": "<error>",
                "passed": False,
                "output_ok": False,
                "cycle_ok": False,
                "determinism_ok": False,
                "output_max_rel": float("inf"),
                "output_max_abs_near_zero": float("inf"),
                "cycle_max_rel": float("inf"),
                "cycle_max_abs_near_zero": float("inf"),
                "parity_targets": [],
                "error": str(exc),
            }

        results.append(detail)
        status = "PASS" if detail["passed"] else "FAIL"
        msg = (
            f"{detail['case']}: {status} | function={detail['function']} | "
            f"output_max_rel={detail['output_max_rel']:.3e} | "
            f"cycle_max_rel={detail['cycle_max_rel']:.3e}"
        )
        if "error" in detail:
            msg += f" | error={detail['error']}"
        print(msg)

        if not detail["passed"]:
            any_fail = True

    sorted_worst = sorted(
        results,
        key=lambda d: max(
            d.get("output_max_rel", 0.0) if np.isfinite(d.get("output_max_rel", 0.0)) else 1.0e99,
            d.get("cycle_max_rel", 0.0) if np.isfinite(d.get("cycle_max_rel", 0.0)) else 1.0e99,
        ),
        reverse=True,
    )

    payload = {
        "total_cases": len(results),
        "passed_cases": int(sum(1 for d in results if d["passed"])),
        "failed_cases": int(sum(1 for d in results if not d["passed"])),
        "failed_case_names": [d["case"] for d in results if not d["passed"]],
        "numba_mode": args.numba,
        "numba_available": bool(is_numba_available()),
        "determinism": bool(args.determinism),
        "tolerances": {
            "output_rel_tol": output_rel_tol,
            "output_abs_floor": output_abs_floor,
            "cycle_rel_tol": cycle_rel_tol,
            "cycle_abs_floor": cycle_abs_floor,
            "near_zero_threshold": near_zero_threshold,
        },
        "worst_cases": sorted_worst[:20],
        "cases": results,
    }

    write_summary_json(args.summary_json, payload)
    write_summary_markdown(args.summary_md, payload)

    return 1 if any_fail else 0


if __name__ == "__main__":
    raise SystemExit(main())
