# refactored_ucb_library.py

This is the main library for the refactored diafiltration workflow.

## What it does

- Loads MATLAB `.mat` files into Python dictionaries
- Builds the first-principles Pyomo model
- Solves the model for parameter estimation
- Computes FIM and simple uncertainty metrics
- Recreates DATA1 and DATA2 paper plots
- Converts conductivity data to concentration when needed

## Main workflow

1. Load experimental data with `loadmat(...)`
2. If the file stores conductivity instead of concentration, call the conductivity helper
3. Fit the model with `solve_model(...)`
4. Optionally compute FIM with `calc_FIM(...)`
5. Recreate paper figures with the DATA1 or DATA2 helper functions

## Important functions

- `loadmat(...)`
- `solve_model(...)`
- `solve_model_B_fix(...)`
- `calc_FIM(...)`
- `run_campaign(...)`
- `run_paper_reproduction(...)`
- `_normalize_conductivity_measurements(...)`
- `_conductivity_to_concentration_series(...)`

## Conductivity conversion

If a dataset has `conductivity_cF = True`, the library:

- saves the original values in `cF_exp_conductivity`
- converts `cF_exp` to concentration
- uses `conductivity_paper.py` directly for the physics-based conversion

## Simple rule

This file is the place where we keep the workflow logic.
The user should not need to work inside the lower-level helper files.

