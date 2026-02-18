# `experiment_load_v24.py`

## Purpose

`experiment_load_v24.py` is the user-facing runner script.  
It wires together:

- experiment loading (`.xlsx`/`.mat`),
- model-option configuration,
- parameter estimation via ParmEst,
- optional DoE/FIM analysis.

It is intentionally small and readable, with heavy lifting delegated to `experiment_dataload_OOP_v24.py`.

## Execution flow

1. Choose input file and selector (`path`, `selector`).
2. Supply optional metadata overrides (`specs`).
3. Toggle stages (`RUN_PARMEST`, `RUN_DOE`).
4. Load one experiment with `load_experiment_easy`.
5. Build `ModelOptions`.
6. Create labeled model via `DiafiltrationExperimentV24`.
7. Run `estimate_parameters_with_parmest_v24`.
8. Run `run_doe_with_pyomo_v24` (optional).

## Key configuration fields

- `mode`: experimental balance mode (`DATA`, `Lag`, `Overflow`).
- `run_mode`: simulation vs estimation.
- `b_form`: transport parameterization for solute flux (`single`, `pervial`, `convection`).
- `nfe`: DAE finite-element count (higher = finer + harder NLP).
- `solver`: ParmEst solver (`ipopt`).

## Typical output

- load status and validation issues,
- theta names (unknown parameter set),
- number of labeled measurements and inputs,
- estimated objective and theta,
- covariance (if converged),
- FIM and D-optimality metric.

## Operational guidance

- Start with `nfe=30` for stable estimation/covariance.
- Increase `nfe` only when needed for accuracy.
- Keep `solver="ipopt"` for ParmEst.
- If DoE struggles, tune NLP settings before increasing model complexity.

## Packaging/refactor notes

- This file should become a thin CLI entrypoint (`console_scripts`) in the package.
- Config values should move into a validated config object or YAML schema.
- Logging should replace `print` for reproducible runs.
