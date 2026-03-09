# `unified_codebase_library.py`

## Role in the system

This is the core unified module. It contains:

- data containers (`ExperimentalData`, `VialData`),
- file ingestion for `.xlsx` and `.mat`,
- validation and preprocessing helpers,
- conductivity-to-concentration transformation hooks,
- dynamic Pyomo model construction,
- ParmEst and DoE integration.

## High-level architecture

1. **Ingest** file -> normalized `ExperimentalData`.
2. **Validate** structure and required fields.
3. **(Optional)** convert conductivity signals to concentration.
4. **Build model** with v24 transport/thermo logic.
5. **Discretize** DAE to NLP.
6. **Label** outputs/errors/unknowns/inputs for ParmEst/DoE.
7. **Estimate parameters** with built-in ParmEst `SSE_weighted`.
8. **Compute covariance** (method fallback sequence).
9. **Run DoE/FIM** utilities.

## Data model

### `SourceType`

- `XLSX`: modern experiment workbooks.
- `MAT`: legacy MATLAB struct files.

### `VialData`

Per-vial segment container, including:

- time-aligned arrays (`time_s`, `mass_g`, retentate/permeate signals),
- optional transformed concentration arrays,
- assay scalars and sample metadata.

### `ExperimentalData`

Experiment-level container:

- provenance (`source`, `filename`, `sheet_name`),
- operating settings (`delP_bar`, `Temp_K`, `Am_cm2`, `rho_g_cm3`),
- feed/dialysate metadata (`M_F0_g`, `C_D_value`, etc.),
- components (`component_names`, `num_components`),
- initial guesses (`Lp0`, `B0`, `sigma0`, `theta0`),
- ordered vial list.

## Ingestion and parsing

### Entry points

- `load_experiment(...)`: file-type dispatch with validation.
- `load_experiment_easy(...)`: convenience wrapper with overrides, optional conversion, optional plotting.

### `.xlsx` path

Core routines:

- `read_excel_sheet`
- `parse_time_series_table`
- `parse_experiment_metadata_kv`
- `parse_icp_assay_table`
- `parse_calibration_block`
- `parse_vial_data_table`
- `build_vial_from_slice`

### `.mat` path

Core routines:

- `read_mat_file`
- `mat_struct_to_python`
- `load_from_mat`

### Initial guess transfer from legacy files

- `infer_initial_guess_db_from_mat_files`
- `apply_initial_guesses_from_db`

## Conductivity conversion utilities

### Main converter

- `apply_conductivity_to_concentration(...)`

This function can transform retentate/permeate conductivity series to concentration using either:

- variant Shedlovsky model (`conductivity_paper.variant_shedlovsky`),
- MSA transport model (`conductivity_paper.msa_transport`).

### Parameter autofill support

The module provides salt parsing and lookup dictionaries for:

- ionic charge,
- limiting equivalent conductivities at 25 C,
- salt-level derived values (`lambda_0`, `lambda_0_cation`, `lambda_0_anion`).

## Process model (v24)

### Core builder stack

- `model_construct_inter_v24` -> delegates to `model_construct_inter_v23` implementation body.
- `model_construct_for_parmest_v24` -> build + discretize + label (without custom objective).

### Equation groups implemented

State variables include:

- retentate concentration (`cF`),
- permeate-side concentration (`cH`),
- vial mass/accumulation (`mV`, `cVmV`, `cV`),
- water and solute fluxes (`Jw`, `Js`),
- interfacial and permeate concentrations (`cIn`, plus `cM`, `cP` in advanced branch).

ODE/DAE groups include:

- concentration dynamics,
- permeate hold-up dynamics,
- flux constitutive equations,
- vial linking constraints,
- boundary conditions.

### `.xlsx`-only advanced transport/thermo branch

Activated when:

- `exp.source == SourceType.XLSX`,
- `ModelOptions.use_advanced_xlsx_transport_thermo` is true.

Adds:

- advanced osmotic-pressure coupling (`delta_pi` expression),
- interfacial/permeate concentration states (`cM`, `cP`),
- thermodynamic diagnostic residuals (partition/electroneutral consistency expressions).

`.mat` experiments keep the legacy branch for backward compatibility.

## Estimation and covariance

### Current ParmEst design

- Uses `Estimator(..., obj_function="SSE_weighted")`.
- Uses Pyomo 6.9.5-style `theta_est` and `cov_est`.
- Exposes outputs through suffix labeling:
  - `experiment_outputs`,
  - `measurement_error`,
  - `unknown_parameters`,
  - `experiment_inputs`.

### Sigma stabilization strategy

To reduce boundary-pathology in covariance:

- estimation can use `sigma_logit` (unconstrained transformed parameter),
- physical sigma is mapped by logistic transform with interior floor/ceiling (`sigma_eps`),
- result postprocessing reattaches physical `sigma` alongside `sigma_logit`.

### Covariance fallback sequence

On `calc_cov=True`, tries:

1. `finite_difference`
2. `reduced_hessian`
3. `automatic_differentiation_kaug`

with temporary `ipopt.opt` settings to improve convergence robustness.

## DoE and FIM

Core helpers:

- `run_doe_with_pyomo_v24`
- `d_optimality`

Design-of-experiments path uses `pyomo.contrib.doe.DesignOfExperiments`, returning:

- FIM matrix,
- D-optimality (`log(det(FIM))`).

## Public APIs to preserve during refactor

For migration safety, treat these as package-facing entry points:

- `load_experiment`
- `load_experiment_easy`
- `apply_conductivity_to_concentration`
- `model_construct_for_parmest_v24`
- `estimate_parameters_with_parmest_v24`
- `run_doe_with_pyomo_v24`
- `DiafiltrationExperimentV24`
- `ModelOptions`, `ParameterGuess`, `ExperimentalData`, `VialData`

## Legacy code notes

Some v23-only helpers are intentionally commented and tagged `TODO(remove)` in current state.
These are retained as audit breadcrumbs during transition and should be removed after package API stabilization.

## Recommended package split (future)

- `unified_diafiltration.data` (containers + ingestion)
- `unified_diafiltration.measurement` (conductivity conversion)
- `unified_diafiltration.model` (Pyomo model construction)
- `unified_diafiltration.estimation` (ParmEst wrappers)
- `unified_diafiltration.design` (DoE/FIM wrappers)
- `unified_diafiltration.cli` (runner/entrypoints)

