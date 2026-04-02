# Release Notes

## Unified Workflow Refactor

This release consolidates the diafiltration analysis stack into a single
maintainable workflow for DATA1, DATA2, and DATA3.

### Highlights

- Introduced a refactored object-oriented workflow centered on:
  - `DataLoader`
  - `ModelOptions`
  - `DiafiltrationExperiment`
  - `UQEngine`

- Unified data ingestion across MAT and XLSX experiment files.

- Routed canonical DATA1 and DATA2 reproduction through `UQEngine` instead of
  scattered script-local orchestration.

- Moved DATA1 and DATA2 validation artifact generation into the workflow layer.

- Added explicit support for:
  - one dataset + many model options
  - many datasets + one model option
  - many datasets + many model options

- Added explicit multi-dataset regression planning with clear separation between:
  - shared parameters
  - dataset-specific parameters

- Folded model comparison into `UQEngine` with AIC/AICc/BIC-style outputs.

- Made experiment-design workflows explicit and distinct:
  - parameter-estimation DoE
  - model-discrimination DoE

- Updated the nightly and regression validation framework to inspect
  workflow-owned DATA1/DATA2 artifact outputs directly.

### Key Files

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py)
- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/output_artifacts.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/output_artifacts.py)
- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/README.md`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/README.md)
- [`UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py)

### Validation Snapshot

Validated on the refactored branch with the following checks:

- `test_unified_workflow_oop.py`
  - workflow loading, model comparison, multi-dataset planning, DoE split, and artifact generation

- `test_unified_codebase_pytest_validation.py`
  - procedural compatibility and broad regression coverage

- `test_simulation_validation_figures.py -m nightly`
  - figure-level validation against extracted DATA1/DATA2 paper images

- `test_data1_data2_nightly.py`
  - canonical DATA1/DATA2 reproduction with workflow-owned artifact checks
  - nightly failures are now interpreted through the documented allowlist policy

### Known Limitations

- Some DATA1/DATA2 numeric parity targets remain on the nightly allowlist in
  [`known_nonpass.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv).
- The workflow and validation infrastructure are in place; remaining work is
  primarily scientific parity improvement rather than architectural cleanup.
