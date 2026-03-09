# DATA1/DATA2 Reproduction Target Matrix

This document tracks published-result reproduction targets for validating the unified framework before full package refactoring.

## Validation policy

- Only figures and tables that appear in `UnifiedFramework/DATA3/published_works/` are valid reproduction targets.
- Generated notebook/script plots are candidate reproductions only after they are mapped back to a published-paper figure/table.
- `PASS` / `FAIL` are reserved for numeric validation against published table values or digitized data extracted from published figures.
- `MISSING_VALUE` means a published target exists but the numeric comparison is not yet complete.
- “Locked” / “Page-matched” in the status column mean published-paper mapping evidence exists; they are not equivalent to numeric pass/fail.

- Parameter tolerance:
  - primary: relative error <= 5%
  - acceptable with explanation: <= 7%
- Objective/WLS/SSE tolerance: relative error <= 5%
- Curves/plots tolerance: normalized RMSE <= 5% of observed range
- Table values tolerance: absolute tolerance based on reported table precision (2-3 significant digits)
- Nightly regression mode: fail on tolerance miss

## Canonical datasets (initial stages)

- DATA1 Stage A canonical MAT:
  - `DATA1_matlab/data_library/data_stru-dataset511.12.mat`
- DATA2 Stage C canonical MAT:
  - `DATA1_matlab/data_library/data_stru-dataset270511.123.mat`
  - `DATA1_matlab/data_library/data_stru-dataset270611.123.mat`

## Source papers

- DATA1 main: `UnifiedFramework/DATA3/published_works/DATA1_main.pdf`
- DATA1 SI: `UnifiedFramework/DATA3/published_works/DATA1_SI.pdf`
- DATA2 main: `UnifiedFramework/DATA3/published_works/DATA2_main.pdf`
- DATA2 SI: `UnifiedFramework/DATA3/published_works/DATA2_SI.pdf`

## DATA1 main targets

| Target ID | Paper item | Artifact type | Dataset scope | Unified artifact output (planned) | Status |
|---|---|---|---|---|---|
| D1-M-F2 | Fig. 2 (all four plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_fig2_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-M-T1 | Table 1 | Table | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data1_table1_side_by_side.csv` | FAIL [numeric run: 20260221-refresh] (artifact: Preliminary FAIL (`20260218-185741`)) |
| D1-M-F3 | Fig. 3 | Figure | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_fig3_model_overlay.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-M-T2 | Table 2 | Table | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data1_table2_side_by_side.csv` | MISSING_VALUE (artifact: Preliminary FAIL (`20260218-185741`)) |
| D1-M-F4 | Fig. 4 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_fig4_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-M-F5 | Fig. 5 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_fig5_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-M-F6 | Fig. 6 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_fig6_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |

## DATA1 SI targets

| Target ID | Paper item | Artifact type | Dataset scope | Unified artifact output (planned) | Status |
|---|---|---|---|---|---|
| D1-S-FS2 | Fig. S2 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS2_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-S-FS3 | Fig. S3 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS3_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-S-FS4 | Fig. S4 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS4_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-S-FS5 | Fig. S5 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS5_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-S-FS6 | Fig. S6 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS6_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |
| D1-S-FS7 | Fig. S7 (all plots) | Figure set | DATA1 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data1_si_figS7_*.png` | MISSING_VALUE (artifact: Locked (`20260221-data1-notebook`)) |

## DATA2 main targets

| Target ID | Paper item | Artifact type | Dataset scope | Unified artifact output (planned) | Status |
|---|---|---|---|---|---|
| D2-M-VARS | Verify manipulated/measured variables | Verification report | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_variable_verification.csv` | NOT_APPLICABLE (reference-only: experimental inputs/measurements) |
| D2-M-T2 | Table 2 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/paper_vs_unified_side_by_side.csv` | FAIL [numeric run: 20260218-190331] (artifact: Preliminary FAIL (`20260218-190331`)) |
| D2-M-F6 | Fig. 6 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_fig6_*.png` | FAIL [numeric run: 20260221-refresh] (artifact: Page-matched (`20260306-023356-paper-pdf-extract`)) |
| D2-M-T3 | Table 3 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_table3_side_by_side.csv` | MISSING_VALUE (artifact: WARNING-baseline (`20260306-mainpaper-table-baseline`)) |
| D2-M-T4 | Table 4 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_table4_side_by_side.csv` | MISSING_VALUE (artifact: WARNING-baseline (`20260306-mainpaper-table-baseline`)) |
| D2-M-F7 | Fig. 7 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_fig7_*.png` | MISSING_VALUE (artifact: Page-matched (`20260306-023356-paper-pdf-extract`)) |
| D2-M-F8 | Fig. 8 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_fig8_*.png` | MISSING_VALUE (artifact: Locked (`20260306-023356-paper-pdf-extract`)) |
| D2-M-T5 | Table 5 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_table5_side_by_side.csv` | MISSING_VALUE (artifact: WARNING-baseline (`20260306-mainpaper-table-baseline`)) |
| D2-M-T6 | Table 6 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_table6_side_by_side.csv` | MISSING_VALUE (artifact: WARNING-baseline (`20260306-mainpaper-table-baseline`)) |
| D2-M-F9 | Fig. 9 | Figure | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_fig9.png` | MISSING_VALUE (artifact: Locked (`20260306-023356-paper-pdf-extract`)) |

## DATA2 SI targets

| Target ID | Paper item | Artifact type | Dataset scope | Unified artifact output (planned) | Status |
|---|---|---|---|---|---|
| D2-S-TS1 | Table S1 | Table | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/data2_si_tableS1_side_by_side.csv` | MISSING_VALUE (artifact: Page-matched (`20260306-023356-paper-pdf-extract`)) |
| D2-S-FS1 | Fig. S1 (both plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS1_*.png` | MISSING_VALUE (artifact: Page-matched (`20260306-023356-paper-pdf-extract`)) |
| D2-S-FS2 | Fig. S2 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS2_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS3 | Fig. S3 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS3_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS4 | Fig. S4 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS4_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS5 | Fig. S5 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS5_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS6 | Fig. S6 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS6_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS7 | Fig. S7 (both plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS7_*.png` | MISSING_VALUE (artifact: Locked (`20260306-cross-verification`)) |
| D2-S-FS8 | Fig. S8 (all plots) | Figure set | DATA2 | `UnifiedFramework/DATA3/results/reproduction/<run_id>/UnifiedFramework/DATA3/figures/data2_si_figS8_*.png` | MISSING_VALUE (artifact: Page-matched (`20260306-023356-paper-pdf-extract`)) |

## Staged execution roadmap

- Stage A: DATA1 canonical dataset (`511.12`) for initial figure/table pipeline.
- Stage B: DATA1 full target set.
- Stage C: DATA2 canonical datasets (`270511.123`, `270611.123`) for variable verification + core tables.
- Stage D: DATA2 full target set.
- Stage E: Convert accepted reproductions into nightly `pytest -m "slow and nightly"` regression tests.
