# GCascadeV5

`GCascadeV5` is the Pythonized runtime for the GCascade V4 port.

The package now has two execution paths:

- `gcascade_v5`
  - the new HDF5-bundle runtime
- `gcascade_v5.legacy`
  - the preserved MAT/CSV-based implementation kept for regression checks and fallback

## Install

From a terminal, move into the cloned GCascadeV5 repository:

```bash
cd /path/to/GCascadeV5
```

Create a virtual environment inside the repository:

```bash
python3 -m venv .venv
```

Activate the virtual environment:

```bash
source .venv/bin/activate
```

Install GCascadeV5, its runtime dependencies, and the tutorial/test tools:

```bash
pip install -e '.[dev]'
```

This installs GCascadeV5 in editable mode, so local code changes are used immediately, together with the packages needed to run the code, the tests, and the tutorial notebook.

Required packages installed by the commands above:

- `numpy`
- `scipy`
- `h5py`
- `matplotlib`
- `notebook` for Jupyter Notebook and `tutorial.ipynb`

Open the tutorial with:

```bash
jupyter notebook tutorial.ipynb
```

Optional performance package:

- `numba`

Install the optional accelerator with:

```bash
pip install -e '.[performance]'
```

## Library Format

The canonical GCascadeV5 library format is now:

- HDF5 data files: `.h5`
- one manifest file: `bundle_manifest.json`

Legacy `.mat` and `.csv` files are no longer the canonical runtime format. They are supported only as import/export compatibility formats.

## Converting A Legacy Library

This repository now uses the following default layout:

- `LibrariesV5-legacy`
  - the preserved MAT/CSV library
- `LibrariesV5`
  - the new canonical HDF5 library used by the standard runtime

Convert a legacy `LibrariesV5-legacy` tree into the canonical HDF5 bundle before using the new runtime:

```python
from pathlib import Path
import gcascade_v5 as gc

gc.convert_legacy_library(
    Path("/path/to/legacy/LibrariesV5-legacy"),
    Path("/path/to/hdf5/LibrariesV5"),
    overwrite=False,
)
```

After conversion, point the runtime at the bundle root:

```python
gc.set_library_path("/path/to/hdf5/LibrariesV5")
```

Inspect the bundle:

```python
print(gc.get_bundle_info())
print(gc.list_generated_variants(1))
```

## Generated Magnetic-Field Variants

`changeMagneticField` now writes HDF5 generated variants and updates the bundle manifest to activate the new file.

```python
import gcascade_v5 as gc

gc.changeMagneticField(1e-18, 0.0, 1)
print(gc.list_generated_variants(1))
```

If you need a legacy MAT export for compatibility work, use:

```python
gc.export_active_cycle_to_legacy_mat(1, "/tmp/cyclespecSL.mat")
```

## Data Paths

Set the bundle root with:

```bash
export GCASCADE_LIB_PATH=/path/to/hdf5/LibrariesV5
```

Optional override for generated HDF5 variants:

```bash
export GCASCADE_GENERATED_LIB_PATH=/path/to/generated
```

Or set them at runtime:

```python
import gcascade_v5 as gc

gc.set_library_path("/path/to/hdf5/LibrariesV5")
gc.set_generated_library_path("/path/to/generated")
```

## Quick Start

```python
import gcascade_v5 as gc

inj = gc.cutoffPowerLaw(gc.energies, gamma=2.2, cutoff=1e7, amp=1e40)
phi = gc.CascadePoint(inj, 0.3)
```

## Legacy Fallback

The preserved pre-pythonization implementation remains available:

```python
import gcascade_v5.legacy as legacy

legacy.set_library_path("/path/to/legacy/LibrariesV5-legacy")
phi_legacy = legacy.CascadePoint(inj, 0.3)
```

This is the reference path used for regression validation.

## Tests

Default tests:

```bash
PYTHONPATH=src python3 -m pytest -q
```

Slow legacy-vs-new cascade regression suite:

```bash
GCASCADE_RUN_SLOW=1 PYTHONPATH=src python3 -m pytest -q
```

## Benchmark

Run the benchmark helper against a legacy library:

```bash
PYTHONPATH=src python3 benchmarks/benchmark_runtime.py --legacy-library LibrariesV5-legacy
```

## Tutorial

The full user workflow, including dependency installation, bundle conversion, generated variants, and legacy fallback checks, is documented in:

- `tutorial.ipynb`
