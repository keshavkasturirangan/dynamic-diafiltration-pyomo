# Panel Validation Checkpoint (DATA1 + DATA2)

Date: 2026-03-06
Scope: target-level figure/table validation using extracted paper artifacts plus regenerated notebook outputs.

## Inputs Used
- Target matrix: `docs/validation/paper_target_matrix.md`
- Consolidated status: `docs/validation/target_validation_status_consolidated.csv`
- Source map: `docs/validation/target_notebook_source_map.csv`
- Paper page index: `docs/validation/target_pdf_page_index.csv`
- DATA1 locked mapping: `docs/validation/data1_panel_checklist.md`
- DATA2 provisional mapping: `docs/validation/data2_panel_checklist_provisional.md`

## Extraction Run (Completed)
Run ID: `20260306-023356-paper-pdf-extract`

Artifacts:
- `results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data1_main/`
- `results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data1_si/`
- `results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_main/`
- `results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/data2_si/`

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
Status: COMPLETE for current scope (page match + locked notebook mapping already in place).

## DATA2 Result
Status: PARTIAL (figure panel locking largely complete; remaining gaps are table artifacts and variable report).

Ready/strong:
- `D2-S-FS1`: explicit notebook section and paper page match.
- `D2-S-FS2` to `D2-S-FS6`: locked via cross-verification outputs for A.1-E.2.
- `D2-S-FS7`: locked to `Bpervial.png` + `Js_Jw_cin.png` by panel-level visual match.
- `D2-M-F7`, `D2-S-FS8`: identified notebook outputs and paper page match.
- `D2-M-F8`: locked to `partition_sensitivity.png` by panel-level visual match.
- `D2-M-F9`: locked to `startup_barplot.png` + `concentrating_residuals_boxplot.png` by panel-level visual match.
- `D2-M-T2`, `D2-M-F6`: numeric validation rows already wired (currently FAIL).
- `D2-M-VARS`: reference-only metadata integrity target (not numeric reproduction).

Still unresolved targets (non-figure artifacts):
- `D2-M-T3`, `D2-M-T4`, `D2-M-T5`, `D2-M-T6`, `D2-S-TS1`

## Next Work Items
1. Populate unified-side values for `D2-M-T3/T4/T5/T6` (warning baselines already exist).
2. Generate dedicated side-by-side artifact for `D2-S-TS1`.
3. Refresh `target_validation_status_consolidated.csv` after numeric table value generation.
