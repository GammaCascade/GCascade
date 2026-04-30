from __future__ import annotations

"""Shared physical constants, cosmology helpers, and energy/redshift grids."""

import numpy as np
from scipy.integrate import quad


# Fundamental constants used by the cascade transport.
Mpc = 3.08568e24  # cm / Mpc
c = 299792.458  # km / s
t = 2.725  # Present-day CMB temperature [K]
mu0 = 1.256637062e-6  # N / A^2
elecmass = 0.51099895e6  # eV
echarge = 1.602176634e-19  # C
sigmaTe = 6.65245873e-29  # m^2

# Background cosmology used for redshift stepping and distance conversion.
H0 = 67.4  # km / s / Mpc
OmegaLambda = 0.685
OmegaM = 0.315


def hubble(z: float | np.ndarray) -> float | np.ndarray:
    """Return the Hubble expansion rate H(z) in km s^-1 Mpc^-1."""
    z_arr = np.asarray(z, dtype=np.float64)
    out = H0 * np.sqrt(OmegaLambda + OmegaM * np.power(1.0 + z_arr, 3.0))
    if np.isscalar(z):
        return float(out)
    return out


def luminosity_distance_mpc(z_start: float) -> float:
    """Integrate the standard luminosity distance out to the source redshift."""
    integral, _ = quad(
        lambda z_val: c / float(hubble(z_val)),
        0.0,
        float(z_start),
        epsabs=0.0,
        epsrel=1.0e-4,
        limit=500,
    )
    return (1.0 + float(z_start)) * integral


def cutoffPowerLaw(
    enerGamma: np.ndarray | float,
    gamma: float,
    cutoff: float,
    amp: float,
) -> np.ndarray:
    """Evaluate the injected spectrum amp * E^-gamma * exp(-E / cutoff)."""
    e = np.asarray(enerGamma, dtype=np.float64)
    return amp * np.power(e, -gamma) * np.exp(-e / cutoff)


def specPlot(spec: np.ndarray) -> tuple[object, object]:
    """Plot E^2 phi(E) on log-log axes for quick spectral inspection."""
    try:
        import matplotlib.pyplot as plt  # local import keeps plotting optional
    except Exception as exc:  # pragma: no cover - optional dependency at runtime
        raise RuntimeError("matplotlib is not available; install it to use specPlot") from exc

    arr = np.asarray(spec, dtype=np.float64)
    if arr.shape != energies.shape:
        raise ValueError(f"Spectrum must have shape {energies.shape}, got {arr.shape}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(energies, energies**2 * arr)
    ax.set_xlabel(r"$E_\gamma$ [GeV]")
    ax.set_ylabel(r"$E_\gamma^2 \phi_\gamma$")
    ax.grid(True, which="both", alpha=0.25)
    return fig, ax


def _round_significant(values: np.ndarray, digits: int) -> np.ndarray:
    """Round a floating grid to a fixed number of significant digits for stable display."""
    arr = np.asarray(values, dtype=np.float64)
    return np.array([float(f"{value:.{digits}g}") for value in arr], dtype=np.float64)


# Redshift windows used by the transport scheme.
_first_step: list[float] = []
for i in range(5):
    count = 10 if i == 0 else 9
    _first_step.extend([1e-6 * (10**i)] * count)

diffuseSteps = np.concatenate(
    [np.array(_first_step, dtype=np.float64), np.full(990, 0.01, dtype=np.float64)]
)
diffuseDistances = np.round(np.cumsum(diffuseSteps), 12)
diffuseDistances[-1] = 10.0
zReg = np.round(np.arange(0.0, 10.0000001, 0.01, dtype=np.float64), 2)

# Shared particle-energy grid.
energies = _round_significant(np.logspace(-1, 12, 300, dtype=np.float64), 15)
energies[0] = 0.1
energies[-1] = 1.0e12
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
