from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
import os
import time
from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.io import loadmat, savemat

# Physical constants (matched to GCascadeV4.wl)
Mpc = 3.08568e24  # cm / Mpc
c = 299792.458  # km / s
t = 2.725  # CMB temperature at z=0 [K]
mu0 = 1.256637062e-6  # N / A^2
elecmass = 0.51099895e6  # eV
echarge = 1.602176634e-19  # C
sigmaTe = 6.65245873e-29  # m^2

# Cosmological parameters from Planck [arXiv:1807.06209]
H0 = 67.4
OmegaLambda = 0.685
OmegaM = 0.315

DEFAULT_LIBRARY_PATH = Path(
    os.getenv("GCASCADE_LIB_PATH", "/Users/antonio/Desktop/Research/GCascade/LibrariesV4")
)
DEFAULT_GENERATED_LIBRARY_PATH = Path(
    os.getenv("GCASCADE_GENERATED_LIB_PATH", "generated_libraries")
)


def logspace(a: float, b: float, n: int) -> np.ndarray:
    """Generate n logarithmically spaced values from 10^a to 10^b."""
    return np.power(10.0, np.linspace(a, b, n, dtype=np.float64))


def hubble(z: float | np.ndarray) -> float | np.ndarray:
    """Hubble rate in (km/s)/Mpc."""
    z_arr = np.asarray(z, dtype=np.float64)
    out = H0 * np.sqrt(OmegaLambda + OmegaM * np.power(1.0 + z_arr, 3.0))
    if np.isscalar(z):
        return float(out)
    return out


# Setup diffuseSteps and diffuseDistances exactly as in V4
_first_step: list[float] = []
for i in range(5):
    count = 10 if i == 0 else 9
    _first_step.extend([1e-6 * (10**i)] * count)

diffuseSteps = np.concatenate(
    [np.array(_first_step, dtype=np.float64), np.full(990, 0.01, dtype=np.float64)]
)
diffuseDistances = np.cumsum(diffuseSteps)

# z-regions of constant CMB/EBL
a = np.arange(0.0, 10.0000001, 0.01, dtype=np.float64)
zReg = np.round(a, 2)

energies = logspace(-1, 12, 300)
dEnergiesGamma = np.diff(energies) / 2.0


EBL_NAME_MAP: dict[int, str] = {
    0: "CMB",
    1: "SL",
    2: "SLhigh",
    3: "SLlow",
    4: "Finke",
    5: "Franc",
    6: "Dom",
}

EBL_DESCRIPTION_MAP: dict[int, str] = {
    0: "CMB only",
    1: "Saldana-Lopez et al. (2021)",
    2: "Saldana-Lopez et al. (2021, 1 sigma upper uncertainty)",
    3: "Saldana-Lopez et al. (2021, 1 sigma lower uncertainty)",
    4: "Finke et al. (2022)",
    5: "Franceschini & Rodighiero (2018)",
    6: "Dominguez et al. (2011)",
}


@dataclass
class GCascadeState:
    library_path: Path
    generated_library_path: Path
    ebl_index: int = 1

    imfp_cmb: np.ndarray | None = None
    imfp_ebl: np.ndarray | None = None
    imfp: np.ndarray | None = None
    extinction_coeffs: np.ndarray | None = None
    cycle_spec: np.ndarray | None = None
    step_size_rows: list[np.ndarray] | None = None
    zreg_index_rows: list[np.ndarray] | None = None

    def reset(self) -> None:
        self.imfp_cmb = None
        self.imfp_ebl = None
        self.imfp = None
        self.extinction_coeffs = None
        self.cycle_spec = None
        self.step_size_rows = None
        self.zreg_index_rows = None

    def ensure_static_loaded(self) -> None:
        if self.step_size_rows is None:
            self.step_size_rows = self._load_step_size_rows()
        if self.zreg_index_rows is None:
            self.zreg_index_rows = self._load_zreg_index_rows()
        if self.imfp_cmb is None:
            self.imfp_cmb = self._load_csv(self.library_path / "pp-IMFPs" / "IMFPcmb.csv")

    def ensure_ebl_loaded(self) -> None:
        self.ensure_static_loaded()
        if self.extinction_coeffs is None or self.cycle_spec is None or self.imfp is None:
            self.load_ebl(self.ebl_index)

    def load_ebl(self, ebl_index: int) -> None:
        if ebl_index not in EBL_NAME_MAP:
            raise ValueError(f"Invalid EBL index: {ebl_index}. Valid values are {sorted(EBL_NAME_MAP)}")

        self.ensure_static_loaded()

        ebl_name = EBL_NAME_MAP[ebl_index]

        if ebl_index == 0:
            self.imfp_ebl = np.zeros_like(self.imfp_cmb, dtype=np.float64)
        else:
            imfp_ebl_path = self.library_path / "pp-IMFPs" / f"IMFPebl{ebl_name}.csv"
            self.imfp_ebl = self._load_csv(imfp_ebl_path)

        self.imfp = self.imfp_cmb + self.imfp_ebl
        self.extinction_coeffs = np.exp(-Mpc * self.imfp)

        cycle_path = self._resolve_cycle_spec_path(ebl_name)
        cycle_raw = self._load_mat_array(cycle_path)
        cycle_raw = self._ensure_cycle_shape(cycle_raw)
        self.cycle_spec = cycle_raw * 1.0e9
        self.ebl_index = ebl_index

    def _resolve_cycle_spec_path(self, ebl_name: str) -> Path:
        generated = self.generated_library_path / "cycle-spec" / f"cyclespec{ebl_name}.mat"
        if generated.exists():
            return generated

        fallback = self.library_path / "cycle-spec" / f"cyclespec{ebl_name}.mat"
        if not fallback.exists():
            raise FileNotFoundError(f"Could not find cycle spec file for EBL '{ebl_name}'")
        return fallback

    def _load_step_size_rows(self) -> list[np.ndarray]:
        arr = self._load_mat_array(self.library_path / "stepSizeArray.mat")
        if arr.shape == (10000, len(diffuseDistances)):
            arr = arr.T
        if arr.shape != (len(diffuseDistances), 10000):
            raise ValueError(f"Unexpected stepSizeArray shape: {arr.shape}")

        rows = [row[row > 0.0].astype(np.float64, copy=False) for row in arr]
        return rows

    def _load_zreg_index_rows(self) -> list[np.ndarray]:
        arr = self._load_mat_array(self.library_path / "zRegIndexArray.mat")
        if arr.shape == (10000, len(diffuseDistances)):
            arr = arr.T
        if arr.shape != (len(diffuseDistances), 10000):
            raise ValueError(f"Unexpected zRegIndexArray shape: {arr.shape}")

        rows = [(row[row > 0].astype(np.int64, copy=False) - 1) for row in arr]
        return rows

    @staticmethod
    def _load_csv(path: Path) -> np.ndarray:
        if not path.exists():
            raise FileNotFoundError(f"Missing CSV file: {path}")
        try:
            return np.loadtxt(path, delimiter=",", dtype=np.float64)
        except ValueError:
            # Some legacy V4 table files use whitespace delimiters despite .csv suffix.
            return np.loadtxt(path, dtype=np.float64)

    @staticmethod
    def _load_mat_array(path: Path) -> np.ndarray:
        if not path.exists():
            raise FileNotFoundError(f"Missing MAT file: {path}")

        payload = loadmat(path)
        keys = [k for k in payload.keys() if not k.startswith("__")]
        if not keys:
            raise ValueError(f"No non-metadata arrays found in MAT file: {path}")

        return np.asarray(payload[keys[0]], dtype=np.float64)

    @staticmethod
    def _ensure_cycle_shape(arr: np.ndarray) -> np.ndarray:
        target = (len(zReg), len(energies), len(energies))
        if arr.shape == target:
            return arr

        legacy = (len(energies), len(energies), len(zReg))
        if arr.shape == legacy:
            return np.transpose(arr, (2, 0, 1))

        raise ValueError(f"Unexpected cycle table shape: {arr.shape}")

    @staticmethod
    def ensure_e_loss_shape(arr: np.ndarray) -> np.ndarray:
        target = (len(zReg), len(energies))
        if arr.shape == target:
            return arr
        if arr.shape == (len(energies), len(zReg)):
            return arr.T
        raise ValueError(f"Unexpected energy-loss table shape: {arr.shape}")


_STATE: GCascadeState | None = None
EBLindex = 1


def _state() -> GCascadeState:
    global _STATE, EBLindex
    if _STATE is None:
        _STATE = GCascadeState(
            library_path=DEFAULT_LIBRARY_PATH,
            generated_library_path=DEFAULT_GENERATED_LIBRARY_PATH,
            ebl_index=EBLindex,
        )
    return _STATE


def reset_state() -> None:
    """Reset all cached runtime tables and state."""
    global _STATE, EBLindex
    _STATE = None
    EBLindex = 1


def set_library_path(path: str | os.PathLike[str]) -> None:
    """Set the read-only path to precomputed V4 libraries."""
    global _STATE
    p = Path(path).expanduser().resolve()
    if _STATE is None:
        _STATE = GCascadeState(p, DEFAULT_GENERATED_LIBRARY_PATH)
    else:
        _STATE.library_path = p
        _STATE.reset()


def set_generated_library_path(path: str | os.PathLike[str]) -> None:
    """Set path where V5-generated tables are written/read."""
    global _STATE
    p = Path(path).expanduser().resolve()
    if _STATE is None:
        _STATE = GCascadeState(DEFAULT_LIBRARY_PATH, p)
    else:
        _STATE.generated_library_path = p


def diffuseDistancesIndex(x: float) -> int:
    """Return the 0-based index in diffuseDistances nearest to x."""
    return int(np.argmin(np.abs(diffuseDistances - float(x))))


def _validate_z_start(z_start: float) -> float:
    z = float(z_start)
    if z < 0.0:
        raise ValueError("zStart must be non-negative")
    if z > 10.0:
        raise ValueError("GCascadeV5 does not allow sources with redshift above z=10")
    return z


def _validate_spectrum_1d(inj: np.ndarray | list[float]) -> np.ndarray:
    arr = np.asarray(inj, dtype=np.float64)
    if arr.shape != energies.shape:
        raise ValueError(
            f"Injected spectrum must have shape {energies.shape}, got {arr.shape}"
        )
    return arr


def _validate_spectrum_2d(inj: np.ndarray | list[list[float]]) -> np.ndarray:
    arr = np.asarray(inj, dtype=np.float64)
    target = (len(diffuseDistances), len(energies))
    if arr.shape != target:
        raise ValueError(f"Injected evolving spectrum must have shape {target}, got {arr.shape}")
    return arr


def _validate_z_distribution(z_distrib: np.ndarray | list[float]) -> np.ndarray:
    arr = np.asarray(z_distrib, dtype=np.float64)
    target = diffuseDistances.shape
    if arr.shape != target:
        raise ValueError(f"Comoving density distribution must have shape {target}, got {arr.shape}")
    return arr


def _inclusive_desc_range(start: float, stop: float, step: float = 1.0e-6) -> np.ndarray:
    if start < stop:
        raise ValueError(f"Expected start >= stop, got start={start}, stop={stop}")
    n = int(round((start - stop) / step))
    return start - step * np.arange(n + 1, dtype=np.float64)


def _build_z_windows(z_max_count: int, step_size: float = 1.0e-6) -> list[np.ndarray]:
    windows: list[np.ndarray] = []
    windows.append(_inclusive_desc_range(diffuseDistances[0], 0.0, step_size))
    for i in range(z_max_count - 1):
        windows.append(_inclusive_desc_range(diffuseDistances[i + 1], diffuseDistances[i], step_size))
    return windows


def _prepare_log_spectrum(spec: np.ndarray) -> np.ndarray:
    with np.errstate(divide="ignore", invalid="ignore"):
        log_spec = np.log10(spec)
    log_spec = np.where(np.isfinite(log_spec), log_spec, -200.0)
    return log_spec


def RedshiftingCycle(injSpectra: np.ndarray, zArrayLocal: np.ndarray) -> np.ndarray:
    """Internal redshifting cycle used by point/diffuse/evolving APIs."""
    stretched_energies = energies * ((1.0 + zArrayLocal[0]) / (1.0 + zArrayLocal[-1]))

    logfunc = _prepare_log_spectrum(np.asarray(injSpectra, dtype=np.float64))
    interp = interp1d(
        energies,
        logfunc,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True,
    )

    stretched = np.empty_like(logfunc)
    stretched[:-1] = interp(stretched_energies[:-1])
    stretched[-1] = interp(energies[-1])  # V4 behavior: do not redshift last energy bin

    return np.where(stretched >= -199.0, np.power(10.0, stretched), 0.0)


def AttenuationCycle(
    injSpectra: np.ndarray,
    zArrayLocal: np.ndarray,
    stepSizeArrayLocal: np.ndarray,
    zRegIndexArrayLocal: np.ndarray,
) -> np.ndarray:
    """Internal attenuation cycle used by point/diffuse/evolving APIs."""
    s = _state()
    s.ensure_ebl_loaded()

    final_result = np.asarray(injSpectra, dtype=np.float64).copy()
    for z_region_index, this_step_size in zip(zRegIndexArrayLocal, stepSizeArrayLocal, strict=True):
        final_result = np.power(s.extinction_coeffs[z_region_index], this_step_size) * final_result

    return RedshiftingCycle(final_result, zArrayLocal)


def CascadeCycle(
    injSpectra: np.ndarray,
    zArrayLocal: np.ndarray,
    stepSizeArrayLocal: np.ndarray,
    zRegIndexArrayLocal: np.ndarray,
) -> np.ndarray:
    """Internal cascade cycle used by point/diffuse/evolving APIs."""
    s = _state()
    s.ensure_ebl_loaded()

    final_result = np.asarray(injSpectra, dtype=np.float64).copy()

    for z_region_index, this_step_size in zip(zRegIndexArrayLocal, stepSizeArrayLocal, strict=True):
        attenuated_spec = np.power(s.extinction_coeffs[z_region_index], this_step_size) * final_result
        delta = final_result - attenuated_spec

        scaled = delta[:, None] * s.cycle_spec[z_region_index]
        trapezoid_kernel = scaled[1:, :] + scaled[:-1, :]
        result = np.sum(trapezoid_kernel * dEnergiesGamma[:, None], axis=0) + attenuated_spec

        final_result = result

    return RedshiftingCycle(final_result, zArrayLocal)


def _luminosity_distance_mpc(z_start: float) -> float:
    integral, _ = quad(
        lambda z_val: c / float(hubble(z_val)),
        0.0,
        z_start,
        epsabs=0.0,
        epsrel=1.0e-4,
        limit=500,
    )
    return (1.0 + z_start) * integral


def cutoffPowerLaw(enerGamma: np.ndarray | float, gamma: float, cutoff: float, amp: float) -> np.ndarray:
    """Evaluate cutoff power law: amp * E^-gamma * exp(-E/cutoff)."""
    e = np.asarray(enerGamma, dtype=np.float64)
    return amp * np.power(e, -gamma) * np.exp(-e / cutoff)


def specPlot(spec: np.ndarray) -> tuple[object, object]:
    """Plot E^2 f(E) on log-log axes; returns (figure, axes)."""
    try:
        import matplotlib.pyplot as plt  # local import to avoid global side effects
    except Exception as exc:  # pragma: no cover - optional dependency at runtime
        raise RuntimeError("matplotlib is not available; install it to use specPlot") from exc

    s = _validate_spectrum_1d(spec)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(energies, energies**2 * s)
    ax.set_xlabel(r"$E_\gamma$ [GeV]")
    ax.set_ylabel(r"$E_\gamma^2 \phi_\gamma$")
    ax.grid(True, which="both", alpha=0.25)
    return fig, ax


def RedshiftPoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)

    z_max_index = diffuseDistancesIndex(z_start) + 1
    z_array = _build_z_windows(z_max_index)

    final_result = inj_spectra.copy()
    for window in reversed(z_array):
        final_result = RedshiftingCycle(final_result, window)

    d_l = _luminosity_distance_mpc(z_start)
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def AttenuatePoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)

    s = _state()
    s.ensure_ebl_loaded()

    z_max_index = diffuseDistancesIndex(z_start) + 1
    z_array = _build_z_windows(z_max_index)
    params = list(
        zip(
            z_array,
            s.step_size_rows[:z_max_index],
            s.zreg_index_rows[:z_max_index],
            strict=True,
        )
    )

    final_result = inj_spectra.copy()
    for z_window, step_row, zreg_row in reversed(params):
        final_result = AttenuationCycle(final_result, z_window, step_row, zreg_row)

    d_l = _luminosity_distance_mpc(z_start)
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def CascadePoint(injSpectraPre: np.ndarray | list[float], zStart: float) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectraPre)
    z_start = _validate_z_start(zStart)

    s = _state()
    s.ensure_ebl_loaded()

    z_max_index = diffuseDistancesIndex(z_start) + 1
    z_array = _build_z_windows(z_max_index)
    params = list(
        zip(
            z_array,
            s.step_size_rows[:z_max_index],
            s.zreg_index_rows[:z_max_index],
            strict=True,
        )
    )

    final_result = inj_spectra.copy()
    for z_window, step_row, zreg_row in reversed(params):
        final_result = CascadeCycle(final_result, z_window, step_row, zreg_row)

    d_l = _luminosity_distance_mpc(z_start)
    return ((1.0 + z_start) ** 2 * final_result) / (4.0 * np.pi * (d_l * Mpc) ** 2)


def _volume_norms(z_start: float, z_distrib: np.ndarray) -> tuple[np.ndarray, list[np.ndarray]]:
    z_max_index = diffuseDistancesIndex(z_start) + 1
    z_interp = interp1d(
        diffuseDistances,
        z_distrib,
        kind="linear",
        bounds_error=False,
        fill_value="extrapolate",
        assume_sorted=True,
    )

    z_grid = diffuseDistances[:z_max_index]
    norms = (c * Mpc * z_interp(z_grid) * diffuseSteps[:z_max_index]) / hubble(z_grid)
    z_array = _build_z_windows(len(norms))
    return norms, z_array


def RedshiftDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    params = list(zip(volume_norms, z_array, strict=True))

    final_result = np.zeros_like(inj_spectra)
    for volume_norm, z_window in reversed(params):
        final_result = RedshiftingCycle(final_result + volume_norm * inj_spectra, z_window)

    return final_result / (4.0 * np.pi)


def AttenuateDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    s = _state()
    s.ensure_ebl_loaded()

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    z_max = len(volume_norms)
    params = list(
        zip(
            volume_norms,
            z_array,
            s.step_size_rows[:z_max],
            s.zreg_index_rows[:z_max],
            strict=True,
        )
    )

    final_result = np.zeros_like(inj_spectra)
    for volume_norm, z_window, step_row, zreg_row in reversed(params):
        final_result = AttenuationCycle(final_result + volume_norm * inj_spectra, z_window, step_row, zreg_row)

    return final_result / (4.0 * np.pi)


def CascadeDiffuse(
    injSpectra: np.ndarray | list[float],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_1d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    s = _state()
    s.ensure_ebl_loaded()

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    z_max = len(volume_norms)
    params = list(
        zip(
            volume_norms,
            z_array,
            s.step_size_rows[:z_max],
            s.zreg_index_rows[:z_max],
            strict=True,
        )
    )

    final_result = np.zeros_like(inj_spectra)
    for volume_norm, z_window, step_row, zreg_row in reversed(params):
        final_result = CascadeCycle(final_result + volume_norm * inj_spectra, z_window, step_row, zreg_row)

    return final_result / (4.0 * np.pi)


def RedshiftEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    params = list(
        zip(volume_norms, z_array, inj_spectra[: len(volume_norms)], strict=True)
    )

    final_result = np.zeros(len(energies), dtype=np.float64)
    for volume_norm, z_window, inj_row in reversed(params):
        final_result = RedshiftingCycle(final_result + volume_norm * inj_row, z_window)

    return final_result / (4.0 * np.pi)


def AttenuateEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    s = _state()
    s.ensure_ebl_loaded()

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    z_max = len(volume_norms)
    params = list(
        zip(
            volume_norms,
            z_array,
            s.step_size_rows[:z_max],
            s.zreg_index_rows[:z_max],
            inj_spectra[:z_max],
            strict=True,
        )
    )

    final_result = np.zeros(len(energies), dtype=np.float64)
    for volume_norm, z_window, step_row, zreg_row, inj_row in reversed(params):
        final_result = AttenuationCycle(final_result + volume_norm * inj_row, z_window, step_row, zreg_row)

    return final_result / (4.0 * np.pi)


def CascadeEvolving(
    injSpectra: np.ndarray | list[list[float]],
    zStart: float,
    zDistrib: np.ndarray | list[float],
) -> np.ndarray:
    inj_spectra = _validate_spectrum_2d(injSpectra)
    z_start = _validate_z_start(zStart)
    z_distrib = _validate_z_distribution(zDistrib)

    s = _state()
    s.ensure_ebl_loaded()

    volume_norms, z_array = _volume_norms(z_start, z_distrib)
    z_max = len(volume_norms)
    params = list(
        zip(
            volume_norms,
            z_array,
            s.step_size_rows[:z_max],
            s.zreg_index_rows[:z_max],
            inj_spectra[:z_max],
            strict=True,
        )
    )

    final_result = np.zeros(len(energies), dtype=np.float64)
    for volume_norm, z_window, step_row, zreg_row, inj_row in reversed(params):
        final_result = CascadeCycle(final_result + volume_norm * inj_row, z_window, step_row, zreg_row)

    return final_result / (4.0 * np.pi)


def changeEBLModel(EBL: int) -> None:
    global EBLindex

    s = _state()
    new_ebl = int(EBL)
    if new_ebl not in EBL_NAME_MAP:
        raise ValueError(f"Invalid EBL index: {new_ebl}")
    if new_ebl == s.ebl_index:
        raise ValueError(f"{EBL_DESCRIPTION_MAP[new_ebl]} is already the current EBL model")

    s.load_ebl(new_ebl)
    EBLindex = s.ebl_index


def _load_ebl_mat_3d(state: GCascadeState, subdir: str, stem: str, ebl_name: str) -> np.ndarray:
    raw = state._load_mat_array(state.library_path / subdir / f"{stem}{ebl_name}.mat")
    return state._ensure_cycle_shape(raw)


def _append_bfield_changelog(changelog_path: Path, ebl: int, b_field: float, gamma: float) -> None:
    changelog_path.parent.mkdir(parents=True, exist_ok=True)
    entry = (
        f"On {datetime.now().isoformat(timespec='seconds')}, the intergalactic magnetic field "
        f"for EBL model {ebl} was changed. "
        f"The new field strength is B(z) = {b_field}*(1+z)^{gamma} Gauss."
    )

    if changelog_path.exists():
        with changelog_path.open("r", encoding="utf-8") as handle:
            lines = [line.rstrip("\n") for line in handle.readlines() if line.strip()]
    else:
        lines = []

    lines.append(entry)
    with changelog_path.open("w", encoding="utf-8") as handle:
        handle.write("\n".join(lines) + "\n")


def changeMagneticField(BField: float, gamma: float, EBL: int) -> None:
    """
    Recompute PP+ICS cycle table for a selected EBL model and activate it.

    Notes:
    - Reads all base physics tables from read-only V4 Libraries.
    - Writes generated cycle tables under generated_library_path/cycle-spec.
    """
    global EBLindex

    s = _state()
    s.ensure_ebl_loaded()

    ebl = int(EBL)
    if ebl not in EBL_NAME_MAP:
        raise ValueError(f"Invalid EBL index: {ebl}")

    ebl_name = EBL_NAME_MAP[ebl]

    new_b_tesla = float(BField) / 10000.0

    dEdtICS = s._load_mat_array(s.library_path / "ics-E-loss-rates" / f"dEdt{ebl_name}.mat")
    dEdtICS = GCascadeState.ensure_e_loss_shape(dEdtICS)

    z_factor = np.power(1.0 + zReg, float(gamma))
    synch_prefactor = ((new_b_tesla * z_factor) ** 2) / (2.0 * mu0)
    gamma_factor = np.sqrt((energies * 1.0e9) ** 2 - elecmass**2) / elecmass

    dEdtsync = (
        (1.0 / echarge)
        * (4.0 / 3.0)
        * sigmaTe
        * synch_prefactor[:, None]
        * c
        * 1000.0
        * np.power(gamma_factor[None, :], 2.0)
    )

    fICS = np.ones((len(zReg), len(energies)), dtype=np.float64)
    denom = dEdtsync[:, 49:] + dEdtICS[:, 49:]
    fICS[:, 49:] = np.divide(
        dEdtICS[:, 49:],
        denom,
        out=np.zeros_like(dEdtICS[:, 49:]),
        where=denom != 0.0,
    )

    ppspec = _load_ebl_mat_3d(s, "pp-spec", "normalizedPPspec", ebl_name)
    otsspec = _load_ebl_mat_3d(s, "on-the-spot-ics-spec", "onthespotICSspec", ebl_name)
    otsspec_weighted = otsspec * fICS[:, :, None]

    cycle_export = np.empty((len(zReg), len(energies), len(energies)), dtype=np.float64)
    d_e = dEnergiesGamma[None, :]

    for i in range(len(zReg)):
        a = ppspec[i, :, :-1]
        b = ppspec[i, :, 1:]
        o1 = otsspec_weighted[i, :-1, :]
        o2 = otsspec_weighted[i, 1:, :]
        cycle_export[i, :, :] = 1.0e9 * ((a * d_e) @ o1 + (b * d_e) @ o2)

    cycle_dir = s.generated_library_path / "cycle-spec"
    cycle_dir.mkdir(parents=True, exist_ok=True)
    target = cycle_dir / f"cyclespec{ebl_name}.mat"

    if target.exists():
        backup = cycle_dir / f"bak_{target.name}"
        if backup.exists():
            backup = cycle_dir / f"bak_{int(time.time())}_{target.name}"
        target.rename(backup)

    savemat(target, {"Expression1": cycle_export})

    # Match V4 behavior by activating the updated cycle table immediately.
    if ebl == s.ebl_index:
        s.cycle_spec = cycle_export * 1.0e9
    else:
        s.load_ebl(ebl)
        EBLindex = s.ebl_index

    _append_bfield_changelog(
        s.generated_library_path / "BfieldChangelog.txt",
        ebl=ebl,
        b_field=float(BField),
        gamma=float(gamma),
    )


# Snake_case aliases
redshift_point = RedshiftPoint
attenuate_point = AttenuatePoint
cascade_point = CascadePoint

redshift_diffuse = RedshiftDiffuse
attenuate_diffuse = AttenuateDiffuse
cascade_diffuse = CascadeDiffuse

redshift_evolving = RedshiftEvolving
attenuate_evolving = AttenuateEvolving
cascade_evolving = CascadeEvolving

change_ebl_model = changeEBLModel
change_magnetic_field = changeMagneticField


__all__ = [
    "Mpc",
    "c",
    "t",
    "mu0",
    "elecmass",
    "echarge",
    "sigmaTe",
    "H0",
    "OmegaLambda",
    "OmegaM",
    "energies",
    "diffuseSteps",
    "diffuseDistances",
    "zReg",
    "dEnergiesGamma",
    "EBLindex",
    "logspace",
    "hubble",
    "cutoffPowerLaw",
    "specPlot",
    "diffuseDistancesIndex",
    "set_library_path",
    "set_generated_library_path",
    "reset_state",
    "RedshiftingCycle",
    "AttenuationCycle",
    "CascadeCycle",
    "RedshiftPoint",
    "AttenuatePoint",
    "CascadePoint",
    "RedshiftDiffuse",
    "AttenuateDiffuse",
    "CascadeDiffuse",
    "RedshiftEvolving",
    "AttenuateEvolving",
    "CascadeEvolving",
    "changeEBLModel",
    "changeMagneticField",
    "redshift_point",
    "attenuate_point",
    "cascade_point",
    "redshift_diffuse",
    "attenuate_diffuse",
    "cascade_diffuse",
    "redshift_evolving",
    "attenuate_evolving",
    "cascade_evolving",
    "change_ebl_model",
    "change_magnetic_field",
]
