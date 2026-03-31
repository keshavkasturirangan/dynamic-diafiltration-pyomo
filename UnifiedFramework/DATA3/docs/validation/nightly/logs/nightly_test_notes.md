# Pytest Notes

This running log captures key outcomes from nightly validation runs.

## 2026-03-06 03:00 local | run_id: 20260306-nightly-local
- Context: Local post-setup validation after adding nightly gate/workflow.
- Validation gate: PASS (`UnifiedFramework/DATA3/scripts/validation/nightly/nightly_validation_gate.py`).
- Gate summary: no unexpected `FAIL`/`MISSING_VALUE`; 7 known non-pass targets allowed as warnings.
- Full nightly pytest: TIMEOUT (bounded local run; solver-heavy execution).
- Action: keep gate as nightly health signal; optimize full regression runtime separately.

## 2026-03-30 local | simulation validation baseline refresh
- Context: Replaced the earlier draft CSV baseline path with a one-time published-paper simulation baseline folder:
  - `UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/`
- Source stack used for the baseline export:
  - `utility.py`
  - `DATA1_model_demo.ipynb`
  - `DATA2_model_demo.ipynb`
  - `DATA2_visualization.ipynb`
  - `run_DATA2_model_variations.py`
- Figure validation report:
  - `UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.csv`
- Figure-validation outcome: `23/23 PASS` for mapped DATA1 main/SI and DATA2 main/SI figure targets using composite-image comparison against extracted paper figures in `20260306-023356-paper-pdf-extract/pdf_extract/`.
- New pytest:
  - `UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py`
- Workflow update:
  - `.github/workflows/nightly-validation.yml` now installs `opencv-python-headless` and runs the figure-validation pytest in addition to the nightly status-gate pytest.
