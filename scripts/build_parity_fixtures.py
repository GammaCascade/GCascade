#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil
from typing import Any

import numpy as np
import gcascade_v5 as gc


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MATRIX_PATH = ROOT / "scripts" / "parity_matrix.json"
DEFAULT_FIXTURES_ROOT = ROOT / "benchmarks" / "v4_reference"

POINT_FUNCTIONS = ["RedshiftPoint", "AttenuatePoint", "CascadePoint"]
DIFFUSE_FUNCTIONS = ["RedshiftDiffuse", "AttenuateDiffuse", "CascadeDiffuse"]
EVOLVING_FUNCTIONS = ["RedshiftEvolving", "AttenuateEvolving", "CascadeEvolving"]


def load_matrix(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def build_base_spectrum(point_case: dict[str, Any]) -> np.ndarray:
    return gc.cutoffPowerLaw(
        gc.energies,
        gamma=float(point_case["gamma"]),
        cutoff=float(point_case["cutoff"]),
        amp=float(point_case["amp"]),
    )


def build_diffuse_distribution(diffuse_case: dict[str, Any]) -> np.ndarray:
    z = gc.diffuseDistances
    kind = diffuse_case["kind"]
    norm = float(diffuse_case["norm"])

    if kind == "gaussian":
        mu = float(diffuse_case["mu"])
        sigma = float(diffuse_case["sigma"])
        return norm * np.exp(-((z - mu) / sigma) ** 2)

    if kind == "power_exp":
        power = float(diffuse_case["power"])
        z_decay = float(diffuse_case["z_decay"])
        return norm * np.power(1.0 + z, power) * np.exp(-z / z_decay)

    raise ValueError(f"Unsupported diffuse case kind: {kind}")


def build_evolving_factor(evolving_case: dict[str, Any]) -> np.ndarray:
    z = gc.diffuseDistances
    kind = evolving_case["kind"]

    if kind == "power":
        return np.power(1.0 + z, float(evolving_case["power"]))

    if kind == "power_exp":
        power = float(evolving_case["power"])
        z_decay = float(evolving_case["z_decay"])
        return np.power(1.0 + z, power) * np.exp(-z / z_decay)

    raise ValueError(f"Unsupported evolving case kind: {kind}")


def ensure_clean_root(fixtures_root: Path, clean: bool) -> None:
    fixtures_root.mkdir(parents=True, exist_ok=True)
    gitkeep = fixtures_root / ".gitkeep"
    if clean:
        for entry in fixtures_root.iterdir():
            if entry.name == ".gitkeep":
                continue
            if entry.is_dir():
                shutil.rmtree(entry)
            else:
                entry.unlink()
    gitkeep.touch(exist_ok=True)


def write_csv(path: Path, arr: np.ndarray) -> None:
    np.savetxt(path, np.asarray(arr, dtype=np.float64), delimiter=",")


def write_case(
    fixtures_root: Path,
    case_id: str,
    meta: dict[str, Any],
    inj: np.ndarray | None = None,
    z_distrib: np.ndarray | None = None,
    inj2d: np.ndarray | None = None,
) -> None:
    case_dir = fixtures_root / case_id
    case_dir.mkdir(parents=True, exist_ok=True)

    with (case_dir / "meta.json").open("w", encoding="utf-8") as handle:
        json.dump(meta, handle, indent=2, sort_keys=True)

    if inj is not None:
        write_csv(case_dir / "inj.csv", inj)
    if z_distrib is not None:
        write_csv(case_dir / "z_distrib.csv", z_distrib)
    if inj2d is not None:
        write_csv(case_dir / "inj2d.csv", inj2d)


def action(name: str, args: list[float | int]) -> dict[str, Any]:
    return {"action": name, "args": args}


def build_pre_actions_for_ebl(sequence: list[int], step_index: int) -> list[dict[str, Any]]:
    # sequence starts from default 1 (current model at package init).
    # step_index references destination at sequence[step_index], with step_index >= 1.
    return [action("changeEBLModel", [int(sequence[i])]) for i in range(1, step_index + 1)]


def make_meta(
    *,
    function: str,
    z_start: float,
    inputs_kind: str,
    ebl_index: int,
    pre_actions: list[dict[str, Any]],
    parity_targets: list[str],
    v4_ref: str,
    case_group: str,
    tags: list[str],
    cycle_sparse_indices: dict[str, list[int]] | None = None,
) -> dict[str, Any]:
    meta: dict[str, Any] = {
        "function": function,
        "z_start": float(z_start),
        "inputs_kind": inputs_kind,
        "ebl_index": int(ebl_index),
        "pre_actions": pre_actions,
        "parity_targets": parity_targets,
        "v4_ref": v4_ref,
        "case_group": case_group,
        "tags": tags,
    }
    if cycle_sparse_indices is not None:
        meta["cycle_sparse_indices"] = cycle_sparse_indices
    return meta


def build_fixtures(fixtures_root: Path, matrix: dict[str, Any]) -> int:
    count = 0
    v4_ref = str(matrix.get("v4_ref", "master"))

    point_cases = matrix["point_cases"]
    diffuse_cases = matrix["diffuse_cases"]
    evolving_cases = matrix["evolving_cases"]
    ebl_indices = [int(v) for v in matrix["ebl_indices"]]

    base_by_point_id: dict[str, np.ndarray] = {}
    diffuse_by_id: dict[str, np.ndarray] = {}
    evolving_factor_by_id: dict[str, np.ndarray] = {}

    for p in point_cases:
        base_by_point_id[p["id"]] = build_base_spectrum(p)
    for d in diffuse_cases:
        diffuse_by_id[d["id"]] = build_diffuse_distribution(d)
    for e in evolving_cases:
        evolving_factor_by_id[e["id"]] = build_evolving_factor(e)

    # Core propagation parity matrix (full stress over all EBL models).
    for ebl in ebl_indices:
        for i in range(3):
            p = point_cases[i]
            d = diffuse_cases[i]
            ev = evolving_cases[i]

            p_id = p["id"]
            d_id = d["id"]
            e_id = ev["id"]

            inj = base_by_point_id[p_id]
            z_distrib = diffuse_by_id[d_id]
            inj2d = evolving_factor_by_id[e_id][:, None] * inj[None, :]

            for fn in POINT_FUNCTIONS:
                case_id = f"prop_{fn}_ebl{ebl}_{p_id}"
                meta = make_meta(
                    function=fn,
                    z_start=float(p["z_start"]),
                    inputs_kind="point",
                    ebl_index=ebl,
                    pre_actions=[],
                    parity_targets=["output"],
                    v4_ref=v4_ref,
                    case_group="propagation",
                    tags=[p_id, f"ebl{ebl}"],
                )
                write_case(fixtures_root, case_id, meta, inj=inj)
                count += 1

            for fn in DIFFUSE_FUNCTIONS:
                case_id = f"prop_{fn}_ebl{ebl}_{p_id}_{d_id}"
                meta = make_meta(
                    function=fn,
                    z_start=float(p["z_start"]),
                    inputs_kind="diffuse",
                    ebl_index=ebl,
                    pre_actions=[],
                    parity_targets=["output"],
                    v4_ref=v4_ref,
                    case_group="propagation",
                    tags=[p_id, d_id, f"ebl{ebl}"],
                )
                write_case(fixtures_root, case_id, meta, inj=inj, z_distrib=z_distrib)
                count += 1

            for fn in EVOLVING_FUNCTIONS:
                case_id = f"prop_{fn}_ebl{ebl}_{p_id}_{d_id}_{e_id}"
                meta = make_meta(
                    function=fn,
                    z_start=float(p["z_start"]),
                    inputs_kind="evolving",
                    ebl_index=ebl,
                    pre_actions=[],
                    parity_targets=["output"],
                    v4_ref=v4_ref,
                    case_group="propagation",
                    tags=[p_id, d_id, e_id, f"ebl{ebl}"],
                )
                write_case(fixtures_root, case_id, meta, inj2d=inj2d, z_distrib=z_distrib)
                count += 1

    # changeEBLModel control sequences with probes.
    p2 = next(p for p in point_cases if p["id"] == "P2")
    d2 = next(d for d in diffuse_cases if d["id"] == "D2")
    inj_p2 = base_by_point_id["P2"]
    z_d2 = diffuse_by_id["D2"]

    for sequence_block in matrix["change_ebl_sequences"]:
        sequence_id = sequence_block["id"]
        sequence = [int(v) for v in sequence_block["sequence"]]
        for step_idx in range(1, len(sequence)):
            target_ebl = sequence[step_idx]
            pre_actions = build_pre_actions_for_ebl(sequence, step_idx)

            for fn in ["AttenuatePoint", "CascadePoint"]:
                case_id = f"ctrl_ebl_{sequence_id}_step{step_idx}_{fn}"
                meta = make_meta(
                    function=fn,
                    z_start=float(p2["z_start"]),
                    inputs_kind="point",
                    ebl_index=target_ebl,
                    pre_actions=pre_actions,
                    parity_targets=["output"],
                    v4_ref=v4_ref,
                    case_group="change_ebl",
                    tags=[sequence_id, f"target{target_ebl}"],
                )
                write_case(fixtures_root, case_id, meta, inj=inj_p2)
                count += 1

            case_id = f"ctrl_ebl_{sequence_id}_step{step_idx}_CascadeDiffuse"
            meta = make_meta(
                function="CascadeDiffuse",
                z_start=float(p2["z_start"]),
                inputs_kind="diffuse",
                ebl_index=target_ebl,
                pre_actions=pre_actions,
                parity_targets=["output"],
                v4_ref=v4_ref,
                case_group="change_ebl",
                tags=[sequence_id, "D2", "P2", f"target{target_ebl}"],
            )
            write_case(fixtures_root, case_id, meta, inj=inj_p2, z_distrib=z_d2)
            count += 1

    # changeMagneticField controls (output probes + sparse cycle parity).
    bfield_cfg = matrix["change_magnetic_field"]
    sparse_idx = bfield_cfg["cycle_sparse_indices"]

    for ebl in [int(v) for v in bfield_cfg["ebl"]]:
        for bfield in [float(v) for v in bfield_cfg["bfield"]]:
            for gamma in [float(v) for v in bfield_cfg["gamma"]]:
                pre_actions: list[dict[str, Any]] = []
                if ebl != 1:
                    pre_actions.append(action("changeEBLModel", [ebl]))
                pre_actions.append(action("changeMagneticField", [bfield, gamma, ebl]))

                combo_label = f"ebl{ebl}_B{bfield:.0e}_g{gamma:.1f}".replace("+", "")

                for i in range(3):
                    p = point_cases[i]
                    d = diffuse_cases[i]
                    ev = evolving_cases[i]

                    p_id = p["id"]
                    d_id = d["id"]
                    e_id = ev["id"]

                    inj = base_by_point_id[p_id]
                    z_distrib = diffuse_by_id[d_id]
                    inj2d = evolving_factor_by_id[e_id][:, None] * inj[None, :]

                    cycle_target = (p_id == "P2")

                    case_id = f"ctrl_bfield_{combo_label}_CascadePoint_{p_id}"
                    meta = make_meta(
                        function="CascadePoint",
                        z_start=float(p["z_start"]),
                        inputs_kind="point",
                        ebl_index=ebl,
                        pre_actions=pre_actions,
                        parity_targets=["output", "cycle_table_sparse"] if cycle_target else ["output"],
                        v4_ref=v4_ref,
                        case_group="change_bfield",
                        tags=[combo_label, p_id],
                        cycle_sparse_indices=sparse_idx if cycle_target else None,
                    )
                    write_case(fixtures_root, case_id, meta, inj=inj)
                    count += 1

                    case_id = f"ctrl_bfield_{combo_label}_CascadeDiffuse_{p_id}_{d_id}"
                    meta = make_meta(
                        function="CascadeDiffuse",
                        z_start=float(p["z_start"]),
                        inputs_kind="diffuse",
                        ebl_index=ebl,
                        pre_actions=pre_actions,
                        parity_targets=["output"],
                        v4_ref=v4_ref,
                        case_group="change_bfield",
                        tags=[combo_label, p_id, d_id],
                    )
                    write_case(fixtures_root, case_id, meta, inj=inj, z_distrib=z_distrib)
                    count += 1

                    case_id = f"ctrl_bfield_{combo_label}_CascadeEvolving_{p_id}_{d_id}_{e_id}"
                    meta = make_meta(
                        function="CascadeEvolving",
                        z_start=float(p["z_start"]),
                        inputs_kind="evolving",
                        ebl_index=ebl,
                        pre_actions=pre_actions,
                        parity_targets=["output"],
                        v4_ref=v4_ref,
                        case_group="change_bfield",
                        tags=[combo_label, p_id, d_id, e_id],
                    )
                    write_case(fixtures_root, case_id, meta, inj2d=inj2d, z_distrib=z_distrib)
                    count += 1

    return count


def main() -> int:
    parser = argparse.ArgumentParser(description="Build deterministic parity fixture inputs/meta from matrix spec.")
    parser.add_argument("--matrix", type=Path, default=DEFAULT_MATRIX_PATH, help="Path to parity matrix JSON")
    parser.add_argument("--fixtures-root", type=Path, default=DEFAULT_FIXTURES_ROOT, help="Fixture root directory")
    parser.add_argument("--clean", action="store_true", help="Delete existing fixture case directories before generation")
    args = parser.parse_args()

    matrix = load_matrix(args.matrix)
    ensure_clean_root(args.fixtures_root, clean=args.clean)
    count = build_fixtures(args.fixtures_root, matrix)

    manifest_path = args.fixtures_root / "_manifest_summary.json"
    summary = {
        "matrix": str(args.matrix),
        "fixtures_root": str(args.fixtures_root),
        "total_cases": count,
    }
    with manifest_path.open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)

    print(f"Generated {count} fixture case directories under {args.fixtures_root}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
