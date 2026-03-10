# DATA2 Panel Checklist (Post-Extraction Status)

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


Date: 2026-03-06
Status: PDF extraction and page matching complete; panel lock completed for Figures 8, 9, S2-S7 (S8 already mapped).
Evidence run: `20260306-023356-paper-pdf-extract`

## Paper Extraction Evidence
- Main paper extraction: `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_main/`
- SI extraction: `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_si/`
- Target page index: `UnifiedFramework/DATA3/docs/validation/target_pdf_page_index.csv`

## Available DATA2 Generated Artifacts (Repo Root)
- `pressure_change_lag.png`
- `pressure_change_overflow.png`
- `calib_curve.png`
- `mass_tc-dat270611.123.png`
- `Js_predict.png`
- `Js_predict0.png`
- `Js_predict1.png`
- `Jw_predict.png`
- `concentrating_residuals_boxplot.png`
- `diluting_residuals_boxplot.png`
- `startup_barplot.png`
- `partition_sensitivity.png`

## Target-by-Target Map

| Target ID | Paper item | Paper pages | Candidate artifact(s) | Confidence | Basis |
|---|---|---|---|---|---|
| D2-M-F6 | Fig. 6 (all plots) | 9,10,11,12,14 | `mass_tc-dat270611.123.png` (plus numeric nRMSE tie-in) | Medium | notebook sectioning + tolerance linkage + paper page match |
| D2-M-F7 | Fig. 7 (all plots) | 13 | `Js_predict.png`, `Js_predict0.png`, `Js_predict1.png`, `Jw_predict.png` | Medium | `DATA2_visualization.ipynb` markdown `Figure 7, S8` + paper page match |
| D2-M-F8 | Fig. 8 (all plots) | 13 | `partition_sensitivity.png` | High | panel-level visual match to paper Figure 8 (A/B/C) |
| D2-M-F9 | Fig. 9 | 15 | `startup_barplot.png` (A), `concentrating_residuals_boxplot.png` (B) | High | panel-level visual match to paper Figure 9 structure and labels |
| D2-S-FS1 | Fig. S1 (both plots) | 12 | `calib_curve.png` | High | `DATA2_visualization.ipynb` section `Figure S1` + paper page match |
| D2-S-FS2 | Fig. S2 (all plots) | 11,13,18 | `mass-dat270611.121.png`, `concentration-dat270611.121.png`, `mass-dat270711.121.png`, `concentration-dat270711.121.png` | High | cross-verification A.1/A.2 outputs (dataset-publication map) |
| D2-S-FS3 | Fig. S3 (all plots) | 14 | `mass-dat270511.221.png`, `concentration-dat270511.221.png`, `mass-dat270511.321.png`, `concentration-dat270511.321.png` | High | cross-verification B.1/B.2 outputs |
| D2-S-FS4 | Fig. S4 (all plots) | 15 | `mass-dat270511.421.png`, `concentration-dat270511.421.png`, `mass-dat270511.921.png`, `concentration-dat270511.921.png` | High | cross-verification C.1/C.2 outputs |
| D2-S-FS5 | Fig. S5 (all plots) | 16 | `mass-dat270511.521.png`, `concentration-dat270511.521.png`, `mass-dat270511.621.png`, `concentration-dat270511.621.png` | High | cross-verification D.1/D.2 outputs |
| D2-S-FS6 | Fig. S6 (all plots) | 17 | `mass-dat270511.721.png`, `concentration-dat270511.721.png`, `mass-dat270511.821.png`, `concentration-dat270511.821.png` | High | cross-verification E.1/E.2 outputs |
| D2-S-FS7 | Fig. S7 (both plots) | 18 | `Bpervial.png` (A), `Js_Jw_cin.png` (B) | High | direct visual match to Figure S7 panel content |
| D2-S-FS8 | Fig. S8 (all plots) | 19 | `Js_predict.png`, `Js_predict0.png`, `Js_predict1.png`, `Jw_predict.png` | Medium | `DATA2_visualization.ipynb` markdown `Figure 7, S8` + paper page match |

## Table Targets (Current State)

| Target ID | Paper item | Paper pages | Status |
|---|---|---|---|
| D2-M-T2 | Table 2 | 9,10,11,13 | Numeric comparison exists (`paper_vs_unified_side_by_side.csv`), currently FAIL |
| D2-M-T3 | Table 3 | 11 | Baseline side-by-side extracted (`UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table3_side_by_side.csv`), WARNING: unified values missing |
| D2-M-T4 | Table 4 | 12 | Baseline side-by-side extracted (`UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table4_side_by_side.csv`), WARNING: unified values missing |
| D2-M-T5 | Table 5 | 14 | Baseline side-by-side extracted (`UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table5_side_by_side.csv`), WARNING: unified values missing |
| D2-M-T6 | Table 6 | 14,15 | Baseline side-by-side extracted (`UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table6_side_by_side.csv`), WARNING: unified values missing |
| D2-S-TS1 | Table S1 | 11 | Missing dedicated side-by-side extraction artifact |
| D2-M-VARS | Verify manipulated/measured variables | 7 | Reference-only target (Table 1 metadata integrity check); not part of numeric reproduction |

## Remaining Lock Steps
1. Populate unified-side values for existing main-paper warning baselines (`D2-M-T3/T4/T5/T6`).
2. Generate missing side-by-side CSV output for `D2-S-TS1`.
