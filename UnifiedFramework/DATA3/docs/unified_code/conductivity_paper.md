# `conductivity_paper.py`

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


## Purpose

`conductivity_paper.py` provides conductivity-based transport/thermodynamic helper equations used by the unified loader’s conversion pipeline.

It contains three model functions:

- `_shedlovsky(...)`: base Shedlovsky relationship helper.
- `variant_shedlovsky(...)`: practical conductivity-to-concentration inversion variant used in this codebase.
- `msa_transport(...)`: mean spherical approximation (MSA)-based transport relation for multi-ion systems.

## How it is used

Primary caller:

- `apply_conductivity_to_concentration(...)` in `unified_codebase_library.py`

Typical flow:

1. Build/auto-fill model parameters from salt metadata.
2. Call `variant_shedlovsky` or `msa_transport`.
3. Map conductivity signals to concentration arrays for each vial.

## Function-level overview

### `_shedlovsky(conc, temp, epsilon, eta, lambda_0, a, z_1, z_2, lambda_0_cation, lambda_0_anion)`

Low-level equation helper implementing the base algebra used in the variant workflow.

### `variant_shedlovsky(...)`

Numerically stable wrapper designed for real experimental conversion workflows.
This is the preferred Shedlovsky interface for the unified code.

### `msa_transport(...)`

MSA transport calculation that accepts ion valency, diameter, diffusion, temperature, viscosity, and dielectric inputs.

## Inputs and units (convention notes)

The unified code currently assumes:

- internally consistent conductivity units with caller-provided parameterization,
- temperature in Kelvin where required,
- dielectric and viscosity values aligned with selected model assumptions.

When packaging, these should be normalized behind a single unit policy and validated schema.

## Numerical behavior notes

- Conversion can become sensitive at low/high concentration extremes.
- Parameter consistency (`lambda_0`, charges, ion-specific terms) strongly affects inversion quality.
- Caller-side safeguards in `unified_codebase_library.py` currently handle fallback and monotonic inversion logic.

## Refactor guidance for package version

- Keep this module pure-function and stateless.
- Add explicit unit checks and parameter dataclasses.
- Add reference tests against published benchmark points.
- Expose a stable API surface:
  - `variant_shedlovsky`
  - `msa_transport`

