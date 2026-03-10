# Workflow Algorithms (Active Unified Code)

## Scope

This document describes the active workflow algorithms for:
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/conductivity_paper.py`

Archived compatibility wrappers are intentionally excluded from the runtime algorithm.

## 1) Top-Level Runtime Algorithm (`unified_codebase_runfile.py`)

### Inputs
- CLI flags (`--profile`, `--file-path`, `--selector`, `--run-doe`, `--calc-cov`, `--plot`, `--convert-to-concentration`, `--nfe`, `--solver`).

### Algorithm
1. Parse CLI arguments.
2. Resolve `run_doe` with profile defaults:
   - `DATA1 -> False`
   - `DATA2 -> False`
   - `DATA3 -> True`
   - explicit CLI value overrides defaults.
3. Branch by profile:
   - `DATA1`: iterate over preset DATA1 `.mat` files (or user override), build `UnifiedPipelineConfigV24` per file.
   - `DATA2`: iterate over preset DATA2 `.mat` files (or user override), build `UnifiedPipelineConfigV24` per file.
   - `DATA3`: build one XLSX-oriented `UnifiedPipelineConfigV24` from preset/defaults or CLI overrides.
4. For each built config, call `run_unified_pipeline_v24(config)`.
5. Print run summary:
   - load/validation state
   - parameter labels and counts
   - optional ParmEst results
   - optional DoE/FIM results.

### Output
- Console output only; all substantive results come from `run_unified_pipeline_v24` return dictionary.

## 2) Unified Pipeline Algorithm (`run_unified_pipeline_v24` in `unified_codebase_library.py`)

### Inputs
- `UnifiedPipelineConfigV24` (file path, selector, load options, model options, stage toggles, solver controls, optional report hook).

### Algorithm
1. Load experiment through file-agnostic loader `load_experiment_easy(...)`:
   - detect source type (`.xlsx/.xls` vs `.mat`),
   - load source-specific structure,
   - validate loaded structure,
   - optionally apply spec overrides,
   - optionally convert conductivity signals to concentrations (XLSX path),
   - optionally quick-plot.
2. Resolve model options:
   - use explicit `config.model_options` if provided,
   - otherwise use source-aware defaults (`MAT` and `XLSX` defaults differ).
3. Build experiment wrapper object: `DiafiltrationExperimentV24(exp, options)`.
4. Build labeled Pyomo model via `get_labeled_model()`.
5. Build metadata record (`_build_unified_run_metadata_v24`).
6. Initialize base output dictionary with:
   - loaded experiment,
   - validation info,
   - resolved options,
   - model/object handles,
   - theta names,
   - label counts,
   - run metadata.
7. If `config.run_parmest` is true:
   - run `_run_unified_estimation_stage_v24` -> `estimate_parameters_with_parmest_v24(...)`.
8. If `config.run_doe` is true:
   - run `_run_unified_doe_stage_v24` -> `run_doe_with_pyomo_v24(...)`.
9. If `config.report_hook` exists, call it with final output dictionary.
10. Return output dictionary.

### Output
- A structured dictionary containing load, model, estimation, and DoE artifacts.

## 3) File-Agnostic Loading Algorithm (`load_experiment` + `load_experiment_easy`)

### 3.1 `load_experiment(path, selector, ...)`
1. Convert `path` to `Path` object.
2. Detect source type:
   - XLSX if `.xlsx/.xls`
   - MAT if `.mat`
   - otherwise fail.
3. Route to source-specific loader:
   - XLSX -> `load_from_xlsx(file_path, sheet_selector=selector, ...)`
   - MAT -> `load_from_mat(file_path, struct_key=selector)`
4. Run `validate_experiment(exp)`.
5. Return `(exp, validation)`.

### 3.2 `load_experiment_easy(...)`
1. If selector absent, apply convenience defaults:
   - XLSX -> first sheet (`0`)
   - MAT -> `"data_stru"`.
2. Call `load_experiment(...)`.
3. Apply optional field overrides (`specs`) directly to `ExperimentalData` object.
4. Optional conductivity-to-concentration conversion:
   - only for XLSX,
   - append informational warning if MAT conversion requested.
5. Optional quick plot.
6. Return `(exp, validation)`.

## 4) XLSX Loader Algorithm (`load_from_xlsx`)

1. Read workbook sheet via selector.
2. Parse metadata and tables:
   - time-series table,
   - vial data table,
   - calibration block,
   - optional assay/ICP table,
   - salt/component hints.
3. Canonicalize columns and data types.
4. Segment rows by vial swaps (`segment_indices_by_swap`).
5. Build one `VialData` object per segment (`build_vial_from_slice`).
6. Attach sheet-level metadata to `ExperimentalData`.
7. Fill defaults/provenance markers for missing optional metadata.
8. Apply initial-guess DB if provided and compatible.
9. Return fully assembled `ExperimentalData`.

## 5) MAT Loader Algorithm (`load_from_mat`)

1. Read `.mat` file (`read_mat_file`).
2. Select struct key (default key behavior if selector absent).
3. Recursively convert MATLAB structs to Python primitives (`mat_struct_to_python`).
4. Map legacy MAT fields into canonical `ExperimentalData` + `VialData` schema.
5. Normalize arrays/scalars/units and legacy naming conventions.
6. Return canonicalized `ExperimentalData`.

## 6) Conductivity-to-Concentration Algorithm (`apply_conductivity_to_concentration`)

### Inputs
- `ExperimentalData` object with conductivity signals in `uS/cm`.
- conversion model choice (`auto`, `variant_shedlovsky`, or `msa`).

### Algorithm
1. Determine temperature (`temp_K` override or `exp.Temp_K`).
2. Choose model:
   - `auto`: infer cation count from salt names,
   - <=2 cations -> variant Shedlovsky,
   - >=3 cations -> MSA.
3. Pre-autofill model parameters:
   - Shedlovsky parameters inferred from salt identity when possible,
   - MSA vector parameters inferred for known salts when possible.
4. For each vial and each signal type (`retentate`, `permeate`):
   - verify units are `uS/cm`,
   - call `_conductivity_to_concentration_series(...)`,
   - write converted concentration array + units + model tag back into vial fields.

### `_conductivity_to_concentration_series(...)` detailed algorithm
1. Import `conductivity_paper.py` helpers.
2. Convert conductivity units (`uS/cm -> mS/cm`).
3. Build pointwise inverse map by root-finding on forward conductivity model:
   - for Shedlovsky: invert `variant_shedlovsky([conc_M], ...)` over concentration bracket,
   - for MSA: invert `msa_transport(...)` over total concentration bracket.
4. For each signal point:
   - if NaN, keep NaN,
   - else invert forward function using `_invert_monotone_1d` and write concentration,
   - if inversion fails, write NaN (best-effort conversion).
5. Convert output units if needed (`M <-> mM`).
6. Return concentration series.

## 7) Conductivity Physics Algorithm (`conductivity_paper.py`)

### 7.1 `_shedlovsky`
1. Compute constants (`B`, `q`, `B1`, `B2`).
2. Compute ionic strength per concentration point.
3. Compute equivalent conductivity via Shedlovsky correction expression.
4. Return equivalent conductivity array.

### 7.2 `variant_shedlovsky`
1. Call `_shedlovsky` for equivalent conductivity.
2. Convert molar concentration to equivalent concentration.
3. Convert equivalent conductivity to specific conductivity.
4. Return specific conductivity in `mS/cm`.

### 7.3 `msa_transport`
1. Construct species number densities for 1, 2, or 3 salts.
2. For each concentration point:
   - compute mobilities and transport numbers,
   - solve alpha roots (bisection),
   - compute hydrodynamic and relaxation corrections,
   - compute species conductivity and sum bulk conductivity.
3. Convert to `mS/cm` and return.

## 8) Modeling + Estimation + DoE Algorithm (Core stage path)

### 8.1 Model assembly
1. Build dynamic model using source-aware v24 constructor wrappers.
2. Apply discretization (`nfe`, scheme).
3. Attach weighted objective and labeled suffixes for ParmEst/DoE.

### 8.2 Parameter estimation (`estimate_parameters_with_parmest_v24`)
1. Build ParmEst experiment list wrappers.
2. Run estimation with configured solver and options.
3. Optionally compute covariance.
4. If configured path fails, use internal fallback/restart helpers.
5. Return estimated parameters, objective, and optional covariance artifacts.

### 8.3 DOE (`run_doe_with_pyomo_v24`)
1. Build DoE object from experiment/model labels.
2. Execute finite-difference sensitivity/FIM workflow.
3. Compute D-optimality metric.
4. Return FIM and summary metrics.

## 9) Active Workflow Summary

1. `unified_codebase_runfile.py` builds profile-specific config.
2. `run_unified_pipeline_v24` orchestrates load -> model -> optional estimate -> optional DoE.
3. `conductivity_paper.py` is used indirectly through conversion wrappers for XLSX conductivity workflows.
4. Archived wrappers are not part of the active algorithmic execution path.
