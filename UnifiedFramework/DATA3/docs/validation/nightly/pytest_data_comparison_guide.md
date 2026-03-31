# Unified Pytest Data Comparison Guide

This guide describes the current validation split:

1. one-time generation of published-paper simulation CSV baselines
2. figure-level validation of those baselines against extracted paper figures
3. unified-code pytest comparisons that will reuse those baselines

## Current source of truth

One-time simulation-side CSV baselines now live in:

- [`simulation_validation_data_files/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files)

Key files there:

- [`README.md`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/README.md)
- [`simulation_validation_manifest.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/simulation_validation_manifest.csv)
- [`simulation_validation_metadata.json`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/simulation_validation_metadata.json)

These CSVs are generated from the published-paper legacy stack:

- [`utility.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
- [`DATA1_model_demo.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_model_demo.ipynb)
- [`DATA2_model_demo.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_model_demo.ipynb)
- [`DATA2_visualization.ipynb`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_visualization.ipynb)
- [`run_DATA2_model_variations.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/run_DATA2_model_variations.py)

Generation entrypoint:

- [`generate_simulation_validation_data_files.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/generate_simulation_validation_data_files.py)

## How figure validation works

Paper evidence root:

- [`pdf_extract/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract)

Figure-validation script:

- [`validate_simulation_validation_data_files.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/validate_simulation_validation_data_files.py)

The script:

1. reads the locked figure mapping from [`target_notebook_source_map.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/target_notebook_source_map.csv)
2. resolves the generated legacy figure PNGs tied to each published figure target
3. builds composite images for multi-panel figures
4. compares those composites against extracted paper images in `pdf_extract/.../images/`
5. writes a report with `best_score`, chosen layout, and `PASS`/`FAIL`

Outputs:

- [`simulation_validation_figure_report.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.csv)
- [`simulation_validation_figure_report.md`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.md)
- [`composites/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/composites)

Current result from the latest run:

- `23/23 PASS` for the mapped DATA1 main/SI and DATA2 main/SI figure targets

## Pytests

Unified code entrypoint/regression module:

- [`test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)

Figure-validation pytest module:

- [`test_simulation_validation_figures.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py)

What the figure-validation pytest does:

- checks that `simulation_validation_data_files` has the committed manifest
- reruns the paper-figure comparison report
- fails if any mapped figure target stops matching the extracted paper figure above threshold

Command:

```bash
MPLCONFIGDIR=/tmp/mplconfig pytest -q UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py -m "nightly"
```

## Unified-code defaults relevant to this flow

User-facing runfile:

- [`unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)

The default conductivity-to-concentration model is now `msa`.

If a user explicitly wants the older conductivity path, they must request:

- `--conductivity-model variant_shedlovsky`

## Practical interpretation

At this point:

- the paper-era simulation CSV baselines are stored in one reusable folder
- those baselines are validated against extracted paper figures
- pytest can rerun the figure validation and report whether the figure mappings still hold

The next layer, which is already partly in place, is to compare unified-code outputs against these same simulation validation CSV baselines instead of relying on ad hoc draft baseline files.
