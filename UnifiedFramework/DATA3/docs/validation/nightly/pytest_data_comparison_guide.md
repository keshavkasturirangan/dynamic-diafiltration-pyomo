# Unified Pytest Data Comparison Guide

This file documents the current pytest-based numeric comparison path for the unified codebase.

## What is being tested

Pytest module:
[`test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)

Current coverage in that module:

- `conductivity_paper.py` snapshot check for `variant_shedlovsky(...)`
- `unified_codebase_runfile.py` preset/config sanity check
- `unified_codebase_library.py` DATA1 filtration MAT path can run end-to-end
- DATA1 Figure 2 unified-vs-legacy numeric comparison using CSV baselines
- DATA1 Figure 2 unified-vs-paper digitized mass-trace comparison
- DATA2 Figure 3 unified-loader-vs-committed-measurement comparison

## Files involved

Unified code under test:

- [`conductivity_paper.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/conductivity_paper.py)
- [`unified_codebase_library.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)
- [`unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)

Legacy/model-side CSV baselines used for numeric comparison:

- [`legacy_paper_csvs/data1_main/fig2/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/unified/legacy_paper_csvs/data1_main/fig2)

Paper-side digitized CSVs created earlier:

- [`digitized_baselines/paper/data1_main/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/paper/data1_main)

One-time legacy exporter that generated the model-side baselines:

- [`generate_legacy_paper_plot_csvs.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/generate_legacy_paper_plot_csvs.py)

## What pytest writes

The legacy-baseline comparison fixture writes a numeric report here:

- [`unified_data1_fig2_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data1_fig2_comparison.csv)

The direct paper-digitization comparison fixture writes a second numeric report here:

- [`unified_data1_fig2_paper_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data1_fig2_paper_comparison.csv)

The current DATA2 comparison fixture writes a third numeric report here:

- [`unified_data2_measurement_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data2_measurement_comparison.csv)

Columns in that report:

- `panel_id`
- `series_id`
- `n_points`
- `mae`
- `max_abs_err`
- `rmse`

These are the main files to inspect if you want to see what the pytest actually compared.

## How the DATA1 comparison works

1. Pytest runs the unified DATA1 MAT pipeline for:
   - `data_stru-dataset501.1.mat`
   - `data_stru-dataset511.12.mat`
2. It extracts unified predicted trajectories from the solved Pyomo model:
   - mass in vial: `mV`
   - retentate concentration: `cF`
   - permeate concentration: `cH`
   - vial concentration: `cV`
3. For the legacy-baseline report:
   - it interpolates unified trajectories onto the time points stored in the legacy baseline CSVs
   - it computes `MAE`, `max_abs_err`, and `RMSE` for each predicted series
4. For the paper-digitized report:
   - it loads the refined paper-side digitized trace CSVs for DATA1 Figure 2 mass panels
   - it matches each paper trace to the best-fitting unified vial trajectory
   - it computes `MAE`, `max_abs_err`, and `RMSE` for that direct paper comparison
5. For the current DATA2 report:
   - it loads the DATA2 MAT experiment through the unified loader
   - it extracts the vial-4 measurement series used by the committed Figure 3 baseline
   - it compares those unified-loaded measurements to the committed CSV in `legacy_paper_csvs/data2_main/fig3`

## Current pytest status meaning

There are three kinds of outcomes in the current unified pytest module:

- `PASS`: the check succeeded
- `XFAIL`: the comparison ran, but the mismatch is a known current gap
- `FAIL`: something unexpected broke

Right now:

- DATA1 Figure 2 mass comparisons are expected to pass
- DATA1 Figure 2 paper-digitized mass comparisons are expected to pass
- DATA1 Figure 2 concentration comparisons are expected to `xfail`
  because unified concentration trajectories are not yet legacy-paper parity
- DATA2 Figure 3 measurement extraction comparison is expected to pass

That means pytest is already comparing against data, even when a result is marked as a known mismatch.

## Command to run

```bash
MPLCONFIGDIR=/tmp/mplconfig pytest -q UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py
```

## If you want to inspect the comparison manually

Open:

- [`unified_data1_fig2_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data1_fig2_comparison.csv)
- [`unified_data1_fig2_paper_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data1_fig2_paper_comparison.csv)
- [`unified_data2_measurement_comparison.csv`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation/unified_data2_measurement_comparison.csv)

Interpretation:

- low `mae` and low `max_abs_err` means the unified prediction is close to the legacy baseline
- in the paper comparison file, `matched_unified_vial` tells you which unified vial trace best aligned with each digitized paper trace
- large errors on the concentration rows mean the unified code is still diverging from the published-era baseline for those series
- the DATA2 file currently verifies one exact loader-backed measurement extraction, not yet a full unified-fit-vs-paper curve comparison

## Important limitation right now

This pytest data-comparison layer is strongest today for DATA1 Figure 2 because it has both legacy CSV baselines and paper-side digitized trace CSVs.

DATA2 is partially covered now through a fast measurement-extraction comparison, but full DATA2 unified-fit comparisons are still heavier because:

- unified DATA2 estimation takes longer
- some legacy DATA2 baseline generation paths are still expensive
- some DATA2 comparisons should probably be split into a slower pytest mark
