# Integration Review Index

Updated: 2026-04-01

## 1) Nightly pytest status (refreshed)

- `docs/validation/nightly/pytest_latest.xml`
- `docs/validation/nightly/pytest_latest_summary.csv`
- Latest run summary: 47 passed, 1 skipped, 290 deselected.

## 2) Digitized baseline expansion outputs

- Manifest (active panel mapping):
  - `docs/validation/nightly/digitized_baselines/manifest.csv`
- Comparison report (latest):
  - `docs/validation/digitized_baselines/comparison_latest.csv`
- New paper/unified panel baseline folders:
  - `docs/validation/nightly/digitized_baselines/paper/`
  - `docs/validation/nightly/digitized_baselines/unified/`

## 3) Visual comparison/share bundle

- Paper evidence extract:
  - `results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/`
- Composite images:
  - `results/simulation_validation/composites/`
- Figure report:
  - `results/simulation_validation/simulation_validation_figure_report.csv`
  - `results/simulation_validation/simulation_validation_figure_report.md`
- Share index:
  - `results/simulation_validation/SHARE_VISUAL_VALIDATION_INDEX.md`

## 4) Files intentionally left untouched in this cleanup

These are existing branch changes and were not modified/reverted by this integration cleanup:

- `tests/regression/test_simulation_validation_figures.py`
- `tests/regression/test_unified_codebase_pytest_validation.py`
- `refactor_prompt.md`
- `refactor_pseudocode.md`

## 5) Cleanup performed

- Removed accidental nested output directory:
  - `UnifiedFramework/DATA3/UnifiedFramework/`
- Removed obsolete, unreferenced baseline files:
  - `docs/validation/nightly/digitized_baselines/paper/data2_f7_panelA.csv`
  - `docs/validation/nightly/digitized_baselines/unified/data2_f7_panelA.csv`
