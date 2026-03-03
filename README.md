# GCascadeV5

`GCascadeV5` is a standalone Python port of `GCascadeV4` with strict function-by-function parity milestones.

## Scope
- V5 code lives in this repository: `/Users/antonio/Desktop/Research/GCascadeV5`
- V4 is read-only and used only for:
  - table reads from `/Users/antonio/Desktop/Research/GCascade/LibrariesV4`
  - reference outputs for parity checks

## Install
```bash
cd /Users/antonio/Desktop/Research/GCascadeV5
python3 -m venv .venv
source .venv/bin/activate
pip install -e '.[dev]'
```

## Data Path
By default, GCascadeV5 reads precomputed tables from:
`/Users/antonio/Desktop/Research/GCascade/LibrariesV4`

Override with:
```bash
export GCASCADE_LIB_PATH=/path/to/LibrariesV4
```

Generated V5 tables (for magnetic-field updates) are written under:
`generated_libraries/`

Override with:
```bash
export GCASCADE_GENERATED_LIB_PATH=/custom/output/path
```

## Public API (V4-compatible names)
- `RedshiftPoint`, `AttenuatePoint`, `CascadePoint`
- `RedshiftDiffuse`, `AttenuateDiffuse`, `CascadeDiffuse`
- `RedshiftEvolving`, `AttenuateEvolving`, `CascadeEvolving`
- `changeEBLModel`, `changeMagneticField`

Snake_case aliases are also provided.

## Parity Workflow
1. Export V4 references manually from Mathematica into `benchmarks/v4_reference/`.
2. Run parity checks:
```bash
python scripts/run_parity.py
```

Each fixture directory must contain:
- `meta.json` (function name + parameters)
- input files (`inj.csv`, optional `z_distrib.csv`, optional `inj2d.csv`)
- `expected.csv`

See `scripts/export_v4_reference_template.wl` for a Mathematica template.

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
