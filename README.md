# GCascadeV5

`GCascadeV5` is the Python implementation of the modern GCascade transport
scheme. The current runtime follows both gamma rays and electrons throughout
the cascade, including pair production, inverse-Compton scattering,
cosmological redshifting, and synchrotron losses.

## Install

From a terminal, move into the cloned repository:

```bash
cd /path/to/GCascadeV5
```

Create and activate a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

Install GCascadeV5 together with the notebook and test dependencies:

```bash
pip install -e '.[dev]'
```

## Download The Runtime Libraries

The `LibrariesV5` runtime bundle is distributed through Zenodo:

- record page: [https://zenodo.org/records/19926427](https://zenodo.org/records/19926427)
- direct ZIP download: [https://zenodo.org/records/19926427/files/LibrariesV5.zip](https://zenodo.org/records/19926427/files/LibrariesV5.zip)

The simplest setup is to download and extract the archive directly into the
repository root so that this folder exists:

```text
/path/to/GCascadeV5/LibrariesV5
```

With that layout, GCascadeV5 should find the library automatically when run
from the repository directory.

Open the tutorial notebook with:

```bash
jupyter notebook tutorial.ipynb
```

## Library Layout

The standard `LibrariesV5` bundle contains:

- `runtime/common.h5`
  - shared redshift-window tables
- `runtime/ebl_*.h5`
  - one HDF5 file per EBL model with:
    - photon attenuation tables
    - pair-production kernels
    - inverse-Compton kernels
    - inverse-Compton loss-rate tables
- `bundle_manifest.json`
  - metadata describing the bundle contents

Set the bundle location with:

```bash
export GCASCADE_LIB_PATH=/path/to/LibrariesV5
```

or at runtime:

```python
import gcascade_v5 as gc

gc.set_library_path("/path/to/LibrariesV5")
```

If you extracted the Zenodo archive into the repository root as
`GCascadeV5/LibrariesV5`, you normally do not need to set this manually.

## Quick Start

```python
import gcascade_v5 as gc

gamma_inj = gc.cutoffPowerLaw(gc.energies, gamma=2.2, cutoff=1e7, amp=1e40)
phi_gamma = gc.CascadePoint(gamma_inj, 0.3)
```

Electron injection is optional and uses the same energy grid:

```python
electron_inj = gc.cutoffPowerLaw(gc.energies, gamma=2.4, cutoff=1e6, amp=1e35)
result = gc.CascadePoint(gamma_inj, 0.3, electronSpectraPre=electron_inj, return_state=True)

phi_gamma = result.gamma
electron_state = result.electron
print(result.diagnostics)
```

For nonzero source redshift, the total energy budget closes as

\[
E_\mathrm{inj}
=
E_{\gamma,\mathrm{final}}
+
E_{e,\mathrm{final}}
+
E_\mathrm{below\ grid}
+
E_\mathrm{synchrotron}
+
E_\mathrm{redshift}.
\]

The diagnostic field `result.diagnostics["redshift_energy_lost"]` records the
energy removed by cosmological expansion.

## Diffuse And Evolving Sources

For diffuse cascades, use:

```python
phi = gc.CascadeDiffuse(gamma_inj, zStart, zDistrib)
```

Optional electron injection uses the same one-dimensional energy grid:

```python
result = gc.CascadeDiffuse(
    gamma_inj,
    zStart,
    zDistrib,
    electronSpectra=electron_inj,
    return_state=True,
)
```

For redshift-evolving injection histories, provide arrays with shape
`(len(gc.diffuseDistances), len(gc.energies))`:

```python
result = gc.CascadeEvolving(
    gamma_history,
    zStart,
    zDistrib,
    electronSpectra=electron_history,
    return_state=True,
)
```

## Magnetic Field

`changeMagneticField(Bfield, gamma, EBLindex)` sets the synchrotron-cooling
law used in the electron transport:

\[
B(z)=B_0(1+z)^\gamma.
\]

Example:

```python
import gcascade_v5 as gc

gc.changeMagneticField(1e-18, 0.0, 1)
```

Restore the default zero-field setup with:

```python
gc.reset_factory_settings()
```

## Inspecting The Bundle

```python
import gcascade_v5 as gc

print(gc.get_bundle_info())
```

## Tests

Run the test suite with:

```bash
PYTHONPATH=src python3 -m pytest -q
```

## Tutorial

All end-user examples are collected in:

- `tutorial.ipynb`
