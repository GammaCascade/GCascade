# GCascadeV5

`GCascadeV5` is a standalone Python port of `GCascadeV4`, developed with strict
function-by-function numerical parity goals before any physics updates.

## Scope And Repository Roles

- V5 development lives in this repository:
  `/path/to/GCascadeV5`
- V4 stays read-only and is used only for:
  - precomputed table reads from
    `/path/to/GCascade/LibrariesV4`
  - generation of reference outputs for parity checks

## Install

```bash
cd /path/to/GCascadeV5
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```

If editable install fails on some systems, try:

```bash
pip install -e . --no-build-isolation
```

## Data Paths

GCascadeV5 reads precomputed V4 tables from a configurable `library_path`.
Recommended: set this explicitly on each machine.

Set table input path before import:

```bash
export GCASCADE_LIB_PATH=/path/to/LibrariesV4
```

Or set it at runtime:

```python
import gcascade_v5 as gc
gc.set_library_path("/path/to/LibrariesV4")
```

`changeMagneticField` reads/writes generated cycle tables from/to
`generated_library_path` (default: `./generated_libraries`). Override with:

```bash
export GCASCADE_GENERATED_LIB_PATH=/custom/output/path
```

or:

```python
gc.set_generated_library_path("/custom/output/path")
```

Inspect active paths:

```python
print(gc.get_library_path())
print(gc.get_generated_library_path())
```

Progress/status printing is enabled by default (useful for long cascade runs).
Disable with:

```bash
export GCASCADE_PROGRESS=0
```

When enabled, all point/diffuse/evolving redshift, attenuation, and cascade APIs
show a dynamic `0% -> 100%` progress bar.

Numba acceleration is enabled by default when installed. Disable with:

```bash
export GCASCADE_NUMBA=0
```

Or at runtime:

```python
gc.set_numba(False)
```

## Quick Start (Minimal Example)

```python
import gcascade_v5 as gc

inj = gc.cutoffPowerLaw(gc.energies, gamma=2.2, cutoff=1e7, amp=1e40)
phi = gc.CascadePoint(inj, 0.3)
```

## Core Arrays And Their Meaning

- `energies`
  - gamma-ray energy grid in GeV
  - length 300
  - logarithmically spaced from `1e-1` to `1e12` GeV
- `diffuseDistances`
  - redshift grid used for diffuse/evolving source integration
  - length 1036
- `zReg`
  - coarser redshift bins (`0` to `10` in steps of `0.01`) used to index
    interaction tables

## Unit Conventions (Important)

- Injected spectrum `inj` (point + diffuse non-evolving):
  - units: `GeV^-1 s^-1`
  - shape: `(300,)`
  - evaluated at `energies`
- Redshift distribution `zDistrib` (diffuse + evolving):
  - units: `cm^-3`
  - shape: `(1036,)`
  - evaluated at `diffuseDistances`
- Evolving injected spectrum `inj2d`:
  - units: `GeV^-1 s^-1`
  - shape: `(1036, 300)`
  - interpreted as `inj2d[z_index, energy_index]`
- Typical outputs:
  - point functions: `GeV^-1 s^-1 cm^-2`
  - diffuse/evolving functions: `GeV^-1 s^-1 cm^-2 sr^-1`

## Public API (V4-Compatible Names)

### Point-source propagation

- `RedshiftPoint(inj, zStart)`
- `AttenuatePoint(inj, zStart)`
- `CascadePoint(inj, zStart)`

### Diffuse non-evolving population

- `RedshiftDiffuse(inj, zStart, zDistrib)`
- `AttenuateDiffuse(inj, zStart, zDistrib)`
- `CascadeDiffuse(inj, zStart, zDistrib)`

### Diffuse evolving population

- `RedshiftEvolving(inj2d, zStart, zDistrib)`
- `AttenuateEvolving(inj2d, zStart, zDistrib)`
- `CascadeEvolving(inj2d, zStart, zDistrib)`

### Advanced model controls

- `changeEBLModel(EBL_index)`
- `changeMagneticField(BField_gauss, gamma, EBL_index)`

Snake_case aliases are available for all public functions.

## EBL Model Index Map

- `0`: CMB only
- `1`: Saldana-Lopez et al. (2021) [default]
- `2`: Saldana-Lopez high
- `3`: Saldana-Lopez low
- `4`: Finke et al. (2022)
- `5`: Franceschini & Rodighiero (2018)
- `6`: Dominguez et al. (2011)

## Tutorial Notebook

For a user-friendly walkthrough with explanations before each command, use:

- `tutorial.ipynb`

It includes setup, units, array formatting, point/diffuse/evolving examples,
EBL/magnetic-field controls, plotting, and parity workflow instructions.

## Exporting Results

Example export pattern:

```python
import numpy as np
import gcascade_v5 as gc

inj = gc.cutoffPowerLaw(gc.energies, gamma=2.2, cutoff=1e7, amp=1e40)
phi = gc.CascadePoint(inj, 0.3)
np.savetxt(
    "point_cascade_flux.csv",
    np.column_stack([gc.energies, phi]),
    delimiter=",",
    header="E_GeV,Phi_GeV^-1_s^-1_cm^-2",
    comments="",
)
```

## Parity Workflow Against V4

1. Export V4 reference fixtures from Mathematica into
   `benchmarks/v4_reference/` (template: `scripts/export_v4_reference_template.wl`).
2. Run parity checks:

```bash
python3 scripts/run_parity.py
```

Each fixture directory must contain:

- `meta.json` (function name + parameters, including `z_start`)
- `expected.csv`
- input files:
  - point: `inj.csv`
  - diffuse: `inj.csv`, `z_distrib.csv`
  - evolving: `inj2d.csv`, `z_distrib.csv`

Default acceptance threshold is `1e-3` relative error (with a small absolute
floor near zero).

## Test Commands

```bash
python3 -m pytest -q
python3 scripts/run_parity.py
```

## Milestone Order

1. `RedshiftPoint`
2. `AttenuatePoint`
3. `CascadePoint`
4. `RedshiftDiffuse`
5. `AttenuateDiffuse`
6. `CascadeDiffuse`
7. `RedshiftEvolving`
8. `AttenuateEvolving`
9. `CascadeEvolving`
10. `changeEBLModel`
11. `changeMagneticField`

No physics updates should be introduced until parity milestones pass.
