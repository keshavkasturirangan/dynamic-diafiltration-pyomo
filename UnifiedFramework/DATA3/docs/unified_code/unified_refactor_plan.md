# Unified Refactor Plan

This note turns the desired unified architecture into a concrete plan against:

- [`unified_codebase_library.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- [`unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)

## Target architecture

The unified system should be organized around four shared layers:

1. input-file adapters
2. common process-model framework
3. common estimation pipeline with ParmEst
4. common FIM / sensitivity / DOE pipeline with Pyomo.DoE

The main scientific difference between `DATA1`, `DATA2`, and future `DATA3` should live in the process model. Paper recreation differences should be layered on top of that, not used as the primary pipeline split.

Object-oriented implementation note:

- the code should stay organized around explicit objects such as:
  - canonical experiment data containers
  - process-model profiles
  - fit-spec objects
  - pipeline configuration objects
- dataset behavior should not be spread across loose conditionals when a named object can own that behavior more cleanly

## Layer 1: Input-file adapters

Goal:

- one path for `.mat`
- one path for `.xlsx`
- both produce the same canonical `ExperimentalData` object

Current code location:

- source detection: [`detect_source_type(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- MAT loading: [`load_from_mat(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- XLSX loading: [`load_from_xlsx(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- unified entrypoint: [`load_experiment_easy(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

Assessment:

- this layer already exists in the right general shape
- the main remaining work is to keep source-specific quirks inside the adapters and not leak them deeper into model construction

Refactor direction:

- keep MAT/XLSX normalization logic at the load layer
- move any remaining source-specific scientific assumptions out of the process-model layer when they are only file-shape accommodations

## Layer 2: Common process-model framework

Goal:

- one common model builder
- process-model features activated or deactivated by profile
- examples:
  - `DATA1`: constant `B`
  - `DATA2`: `B` depends on concentration or transport formulation
  - `DATA3`: future process-model extensions

Current code location:

- mode enum: [`ExperimentMode`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- model options: [`ModelOptions`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- main model builder: [`model_construct_inter_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

Assessment:

- this layer exists, but the current split is still partly expressed in legacy terms such as `ExperimentMode.DATA`, `ExperimentMode.LAG`, and `ExperimentMode.OVERFLOW`
- the code already supports different process-model behavior through `mode`, `b_form`, and related options
- this should evolve toward a clearer process-model-profile abstraction
- the new `ProcessModelProfile` object is the first low-risk step in that direction

Refactor direction:

- keep one common model builder
- make the primary scientific split a process-model profile instead of a paper-profile or file-type branch
- treat `DATA1`, `DATA2`, and `DATA3` as named process-model profiles that configure:
  - governing equations
  - state set
  - parameterization
  - experiment structure

Recommended next step:

- introduce a process-model-profile vocabulary in docs and config without breaking the current enums yet
- carry the process-model profile through run metadata so reports can separate scientific model choice from file source and paper recreation overlays

## Layer 3: Common ParmEst estimation pipeline

Goal:

- one estimation engine using ParmEst
- profile-specific hooks define:
  - parameters to estimate
  - bounds
  - initial guesses
  - measurement mapping
  - any paper-recreation weighting rules

Current code location:

- ParmEst builder: [`model_construct_for_parmest_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- legacy grouped-fit builder: [`model_construct_for_legacy_fit_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- estimation runner: [`estimate_parameters_with_parmest_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

Assessment:

- the unified code already has one shared estimation entrypoint
- the main open design issue is separating:
  - process-model-specific estimation setup
  - paper-recreation-specific fit overrides

Refactor direction:

- keep one shared ParmEst stage
- make DATA1/DATA2 differences flow through process-model profile hooks
- keep paper-era parity rules as overlays on top of the estimation stage, not as their own independent pipeline family
- use a narrow model-setting name such as `ParameterTreatmentMode` for fixed-vs-estimated parameter handling, instead of a broader name that sounds like the full workflow

Recommended next step:

- define an explicit fit-spec layer for each process-model profile
- let paper profiles modify the fit spec rather than replace the whole estimation path

## Layer 4: Common FIM / sensitivity / DOE pipeline

Goal:

- one shared Pyomo.DoE stage
- profile-specific setup defines:
  - design variables
  - measured outputs
  - uncertain parameters
  - operating ranges

Current code location:

- DoE runner: [`run_doe_with_pyomo_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- pipeline hook: [`_run_unified_doe_stage_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

Assessment:

- this layer is already shared structurally
- it still needs profile-driven setup rather than ad hoc assumptions baked into the experiment/model state

Refactor direction:

- keep one common DoE/FIM stage
- make the process-model profile responsible for declaring the analysis problem

## Where paper profiles fit

Paper profiles are still useful, but they should be treated as recreation overlays.

Current code location:

- runner presets in [`unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)
- `paper_profile` fields in [`ModelOptions`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py) and [`UnifiedPipelineConfigV24`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

Role of paper profiles:

- recreate published figures and tables
- apply paper-specific bounds, initialization, weighting, and artifact expectations
- never become the main scientific split between pipelines

## Concrete changes to make next

1. Add a named object-oriented process-model-profile layer in config.
2. Map current legacy `ExperimentMode` settings onto that profile layer.
3. Move DATA1/DATA2-specific fit setup into profile-specific fit specs.
4. Keep `paper_profile` as a recreation overlay only.
5. Make run metadata report both:
   - input source profile
   - process-model profile
   - paper-recreation profile

## Implementation status

The refactor now has a dedicated orchestration folder:

- [`workflow/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow)
- primary module: [`unified_workflow.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py)

Current organization in that module:

- `DataLoader`: file-agnostic MAT/XLSX loading requests
- `DiafiltrationExperiment`: one dataset + one `ModelOptions`
- `UQEngine`: estimation, uncertainty summaries, model comparison, parameter-estimation DoE, and model-discrimination scoring

Compatibility behavior:

- [`run_unified_pipeline_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py) now routes through `UQEngine`
- existing tests and scripts can keep using the procedural API while the OO module becomes the maintained orchestration surface

## Refactor guardrails

Use the existing pytests as the acceptance check for architectural cleanup.

These are the minimum guardrails:

- [`test_utility_input_contract.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_utility_input_contract.py)
- [`test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)
- [`test_simulation_validation_figures.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py)

Current expected state during refactor:

- `test_utility_input_contract.py`: fully green
- `test_unified_codebase_pytest_validation.py`: green except for the known DATA1 concentration `xfail`
- `test_simulation_validation_figures.py`: green

Interpretation:

- if these tests stay in their current expected state, the refactor is reorganizing the code without breaking DATA1/DATA2 pipeline behavior
- if a new failure appears outside the known DATA1 concentration `xfail`, treat it as a refactor regression until proven otherwise

## What not to do

- do not fork separate end-to-end DATA1 and DATA2 pipelines
- do not let `.mat` versus `.xlsx` determine the scientific model
- do not let paper recreation logic become the primary architectural split

## Current north star

The unified code should read like:

- file-type adapter selects how experimental data are loaded
- process-model profile selects the scientific model
- shared ParmEst stage estimates parameters
- shared DoE stage performs FIM / sensitivity analysis
- optional paper profile recreates published outputs and validation behavior
