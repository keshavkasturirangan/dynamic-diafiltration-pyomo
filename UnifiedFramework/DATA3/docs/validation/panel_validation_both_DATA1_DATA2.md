# Panel Validation Checkpoint (DATA1 + DATA2)

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


Date: 2026-03-06
Scope: target-level figure/table validation using extracted paper artifacts plus regenerated notebook outputs.

Interpretation note:
- “Locked” in this document means the generated artifact has been mapped to a figure/table that exists in `UnifiedFramework/DATA3/published_works/`.
- It does not mean the target has numerically passed reproduction.
- Numeric validation remains separate and is reported through `target_validation_status_consolidated.csv`.

## Inputs Used
- Target matrix: `UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md`
- Consolidated status: `UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv`
- Source map: `UnifiedFramework/DATA3/docs/validation/target_notebook_source_map.csv`
- Paper page index: `UnifiedFramework/DATA3/docs/validation/target_pdf_page_index.csv`
- DATA1 locked mapping: `UnifiedFramework/DATA3/docs/validation/data1_panel_checklist.md`
- DATA2 provisional mapping: `UnifiedFramework/DATA3/docs/validation/data2_panel_checklist_provisional.md`

## Extraction Run (Completed)
Run ID: `20260306-023356-paper-pdf-extract`

Artifacts:
- `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data1_main/`
- `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data1_si/`
- `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_main/`
- `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_si/`

Extraction totals:
- DATA1 main: 12 pages, 9 embedded images
- DATA1 SI: 10 pages, 30 embedded images (+OCR page text)
- DATA2 main: 20 pages, 12 embedded images
- DATA2 SI: 20 pages, 17 embedded images

## Current Validation State
- Target-to-paper-page matching: complete for all 32 targets (`target_pdf_page_index.csv`).
- DATA1: locked figure mapping remains valid and unchanged.
- DATA2: paper page matches are complete and most figure panel mappings are now locked.

## DATA1 Result
Status: COMPLETE for published-paper mapping scope (page match + published-paper artifact mapping already in place), but not complete for full numeric reproduction.

## DATA2 Result
Status: PARTIAL (published-paper figure mapping is largely complete; numeric reproduction still has table gaps and failures).

Ready/strong for published-paper mapping:
- `D2-S-FS1`: explicit notebook section and paper page match.
- `D2-S-FS2` to `D2-S-FS6`: locked via cross-verification outputs for A.1-E.2.
- `D2-S-FS7`: locked to `Bpervial.png` + `Js_Jw_cin.png` by panel-level visual match.
- `D2-M-F7`, `D2-S-FS8`: identified notebook outputs and paper page match.
- `D2-M-F8`: locked to `partition_sensitivity.png` by panel-level visual match.
- `D2-M-F9`: locked to `startup_barplot.png` + `concentrating_residuals_boxplot.png` by panel-level visual match.
- `D2-M-T2`, `D2-M-F6`: numeric validation rows already wired (currently FAIL).
- `D2-M-VARS`: reference-only metadata integrity target (not numeric reproduction).

Still unresolved for numeric reproduction or missing published-paper comparison artifacts:
- `D2-M-T3`, `D2-M-T4`, `D2-M-T5`, `D2-M-T6`, `D2-S-TS1`

## Next Work Items
1. Populate unified-side values for `D2-M-T3/T4/T5/T6` (warning baselines already exist).
2. Generate dedicated side-by-side artifact for `D2-S-TS1`.
3. Refresh `target_validation_status_consolidated.csv` after numeric table value generation.
