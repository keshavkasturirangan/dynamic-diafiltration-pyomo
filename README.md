# Diafiltration Comprehensive Workflow

This repository contains the refactored Pyomo-based diafiltration analysis
workflow for DATA1, DATA2, and DATA3.

Canonical development path:

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/)

Primary workflow module:

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/unified_workflow.py)

Workflow overview:

- `DataLoader`
  - loads MAT and XLSX experiment files
  - normalizes them into the shared `ExperimentalData` representation

- `ModelOptions`
  - stores one candidate model configuration
  - carries process-model profile, paper-profile overlay, estimation settings, and DoE settings

- `DiafiltrationExperiment`
  - binds one dataset to one `ModelOptions`
  - constructs the labeled Pyomo model used for simulation, ParmEst, and DoE

- `UQEngine`
  - orchestrates parameter estimation, uncertainty summaries, model comparison, plotting/validation artifacts, and next-experiment workflows
  - supports:
    - one dataset + many model options
    - many datasets + one model option
    - many datasets + many model options
  - distinguishes:
    - parameter-estimation DoE
    - model-discrimination DoE

## Where To Start

For the current workflow entrypoint summary, start here:

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/README.md`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/workflow/README.md)

For the compatibility procedural entrypoint:

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

For preset runs and paper-oriented configurations:

- [`UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)

For DATA1/DATA2 canonical reproduction:

- [`UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py)

## Validation

High-signal validation paths:

- [`UnifiedFramework/DATA3/tests/regression/test_unified_workflow_oop.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_workflow_oop.py)
- [`UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)
- [`UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py)
- [`UnifiedFramework/DATA3/tests/regression/test_data1_data2_nightly.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_data1_data2_nightly.py)

Nightly-validation notes:

- [`UnifiedFramework/DATA3/docs/validation/nightly/pytest_data_comparison_guide.md`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/pytest_data_comparison_guide.md)
- [`UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv)

## Environment

Suggested baseline environment:

- Python 3.11
- Pyomo
- NumPy
- pandas
- SciPy
- matplotlib
- CasADi

Example conda setup:

```bash
conda create -n dynamic-diafiltration -c anaconda -c conda-forge -c IDAES-PSE python=3.11 numpy matplotlib pandas scipy idaes-pse scikit-learn
conda activate dynamic-diafiltration
pip install casadi
```

## Legacy And Reference Paths

Legacy root-level scripts and notebooks are still present as reference or
compatibility paths, but they are no longer the primary development surface.

Examples:

- [`utility.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
- [`DATA1_model_demo.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_model_demo.ipynb)
- [`DATA2_model_demo.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_model_demo.ipynb)
- [`DATA2_visualization.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_visualization.ipynb)

## Current Status

The workflow refactor is in place and validated. Remaining work is primarily
scientific parity improvement for the still-allowlisted DATA1/DATA2 numeric
mismatches, not further workflow reorganization.
