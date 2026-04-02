# Unified Workflow Module

This folder contains the refactored object-oriented orchestration layer for the
DATA1/DATA2/DATA3 diafiltration workflow.

Primary module:

- [`unified_workflow.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py)

Supporting artifact/output helpers:

- [`output_artifacts.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/output_artifacts.py)

## Main objects

- `DataLoader`
  - loads MAT and XLSX experiment files
  - normalizes them to the shared `ExperimentalData` representation
  - keeps source-specific parsing at the load layer

- `ModelOptions`
  - remains the lightweight configuration object for one candidate model
  - carries process-model profile, paper profile, estimation options, and DoE settings

- `DiafiltrationExperiment`
  - binds one dataset to one `ModelOptions`
  - builds labeled Pyomo models for simulation, ParmEst, and DoE
  - exposes measured/simulated data in plotting-friendly forms

- `UQEngine`
  - orchestrates estimation, uncertainty summaries, model comparison, plotting bundles, and experiment-design workflows
  - supports:
    - one dataset + many model options
    - many datasets + one model option
    - many datasets + many model options
  - makes shared vs dataset-specific regression treatment explicit through `ParameterScope` and `RegressionPlan`
  - distinguishes:
    - parameter-estimation DoE
    - model-discrimination DoE

## Workflow shape

The intended call pattern is:

`DataLoader -> UQEngine`

Typical usage:

1. create `DatasetRequest` objects
2. load them with `DataLoader`
3. configure one or more `ModelOptions`
4. run `UQEngine.fit_parameters(...)`, `compare_models(...)`, or `run_full_workflow(...)`
5. build DATA1/DATA2 validation artifacts through the workflow layer when needed

## Compatibility

The procedural compatibility API is still available:

- [`run_unified_pipeline_v24(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

That function now routes through `UQEngine`, so existing scripts and tests can
keep their current entrypoint while the workflow module serves as the maintained
orchestration surface.

## Validation hooks

Key regression coverage:

- [`test_unified_workflow_oop.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_workflow_oop.py)
- [`test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)
- [`test_data1_data2_nightly.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_data1_data2_nightly.py)
- [`test_simulation_validation_figures.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py)
