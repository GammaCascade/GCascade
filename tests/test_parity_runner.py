from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
RUN_PARITY_PATH = ROOT / "scripts" / "run_parity.py"


def load_run_parity_module():
    spec = importlib.util.spec_from_file_location("run_parity", RUN_PARITY_PATH)
    assert spec is not None
    assert spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_validate_meta_rejects_missing_required_fields():
    rp = load_run_parity_module()

    with pytest.raises(ValueError):
        rp.validate_meta({"function": "RedshiftPoint"})


def test_apply_pre_actions_rejects_unknown_action():
    rp = load_run_parity_module()

    with pytest.raises(ValueError):
        rp.apply_pre_actions([
            {"action": "unknownAction", "args": []},
        ])


def test_apply_pre_actions_calls_handlers_in_order(monkeypatch):
    rp = load_run_parity_module()
    calls: list[tuple[str, tuple[float | int, ...]]] = []

    monkeypatch.setattr(rp, "changeEBLModel", lambda ebl: calls.append(("changeEBLModel", (int(ebl),))))
    monkeypatch.setattr(
        rp,
        "changeMagneticField",
        lambda b, g, ebl: calls.append(("changeMagneticField", (float(b), float(g), int(ebl)))),
    )

    rp.apply_pre_actions(
        [
            {"action": "changeEBLModel", "args": [4]},
            {"action": "changeMagneticField", "args": [1e-16, 2.0, 4]},
        ]
    )

    assert calls == [
        ("changeEBLModel", (4,)),
        ("changeMagneticField", (1e-16, 2.0, 4)),
    ]
