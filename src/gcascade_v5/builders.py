from __future__ import annotations

from datetime import datetime, timezone
from pathlib import Path

import h5py
import numpy as np

from . import bundle, legacy
from .config import progress_marks, status


def _format_float_token(value: float) -> str:
    numeric = float(value)
    if numeric == 0.0:
        return "0e0"
    mantissa, exponent = f"{numeric:.2e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    return f"{mantissa}e{int(exponent)}"


def _variant_id(ebl_name: str, b_field: float, gamma: float) -> str:
    return f"ebl_{ebl_name}_B_{_format_float_token(b_field)}_gamma_{_format_float_token(gamma)}"


def generate_magnetic_field_variant(
    bundle_root: str | Path,
    generated_root: str | Path,
    *,
    ebl_index: int,
    b_field: float,
    gamma: float,
) -> Path:
    root = Path(bundle_root).expanduser().resolve()
    output_root = Path(generated_root).expanduser().resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    ebl_name = legacy.EBL_NAME_MAP[int(ebl_index)]
    builder_path = bundle.resolve_builder_ebl_path(root, int(ebl_index))

    new_b_tesla = float(b_field) / 10000.0
    z_factor = np.power(1.0 + legacy.zReg, float(gamma))
    synch_prefactor = ((new_b_tesla * z_factor) ** 2) / (2.0 * legacy.mu0)
    gamma_factor = np.sqrt((legacy.energies * 1.0e9) ** 2 - legacy.elecmass**2) / legacy.elecmass

    with h5py.File(builder_path, "r") as handle:
        d_edt_ics = np.asarray(handle["dEdt_ics"], dtype=np.float64)
        pp_packed = handle["pp_packed"]
        ots_packed = handle["ots_packed"]

        d_edt_sync = (
            (1.0 / legacy.echarge)
            * (4.0 / 3.0)
            * legacy.sigmaTe
            * synch_prefactor[:, None]
            * legacy.c
            * 1000.0
            * np.power(gamma_factor[None, :], 2.0)
        )

        f_ics = np.ones((len(legacy.zReg), len(legacy.energies)), dtype=np.float64)
        denom = d_edt_sync[:, 49:] + d_edt_ics[:, 49:]
        f_ics[:, 49:] = np.divide(
            d_edt_ics[:, 49:],
            denom,
            out=np.zeros_like(d_edt_ics[:, 49:]),
            where=denom != 0.0,
        )

        packed_cycle = np.empty((len(legacy.zReg), bundle.PACKED_TRIANGULAR_SIZE), dtype=np.float64)
        d_e = legacy.dEnergiesGamma[None, :]
        marks = progress_marks(len(legacy.zReg))
        for z_idx in range(len(legacy.zReg)):
            pp_slice = bundle._unpack_lower_slice(np.asarray(pp_packed[z_idx], dtype=np.float64))
            ots_slice = bundle._unpack_lower_slice(np.asarray(ots_packed[z_idx], dtype=np.float64))
            ots_weighted = ots_slice * f_ics[z_idx][:, None]

            a = pp_slice[:, :-1]
            b = pp_slice[:, 1:]
            o1 = ots_weighted[:-1, :]
            o2 = ots_weighted[1:, :]
            cycle_raw = 1.0e9 * ((a * d_e) @ o1 + (b * d_e) @ o2)
            cycle_runtime = cycle_raw * 1.0e9
            packed_cycle[z_idx] = bundle._pack_lower_cube(cycle_runtime[None, :, :])[0]
            if (z_idx + 1) in marks:
                status(f"changeMagneticField progress: {z_idx + 1}/{len(legacy.zReg)}")

    variant_id = _variant_id(ebl_name, float(b_field), float(gamma))
    target = output_root / f"{variant_id}.h5"
    with h5py.File(target, "w") as handle:
        handle.create_dataset("cycle_packed", data=packed_cycle, compression="lzf", shuffle=True)
        handle.attrs["ebl_index"] = int(ebl_index)
        handle.attrs["ebl_name"] = ebl_name
        handle.attrs["b_field_gauss"] = float(b_field)
        handle.attrs["gamma"] = float(gamma)
        handle.attrs["variant_id"] = variant_id
        handle.attrs["created_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    bundle.register_generated_variant(
        root,
        ebl_index=int(ebl_index),
        file_path=target,
        parameters={"b_field_gauss": float(b_field), "gamma": float(gamma)},
        variant_id=variant_id,
    )
    return target
