# Validation Docs

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


This folder tracks DATA1/DATA2 published-result reproduction for unified-code validation.

## Validation semantics

- Source of truth for validation targets: figures and tables that appear in `UnifiedFramework/DATA3/published_works/` only.
- `PASS` / `FAIL`: reserved for numeric validation against published table values or digitized data extracted from published paper figures.
- `MISSING_VALUE`: a published target exists, but the numeric comparison is incomplete because either the paper-side numeric baseline or unified-side value is still missing.
- `NOT_APPLICABLE`: the target is not currently part of numeric gating. In practice this often means the target is published-paper mapped/page-matched but has not yet been converted into a numeric comparison.
- “Locked” or “visual/artifact validation” should be read as: the generated artifact has been mapped back to a figure/table that exists in `UnifiedFramework/DATA3/published_works/`; it is not, by itself, a numeric pass.

## Files

- `paper_target_matrix.md`: required figure/table targets and stage plan.
- `paper_reference_values_template.csv`: side-by-side paper/unified reference template.
- `target_validation_status_consolidated.csv`: target-level status (`PASS/FAIL/MISSING_VALUE/NOT_APPLICABLE`).
- `target_notebook_source_map.csv`: target-to-artifact evidence mapping.
- `target_pdf_page_index.csv`: target-to-paper-page evidence index.
- `panel_validation_both_DATA1_DATA2.md`: current high-level checkpoint.
- `artifact_map_DATA1_DATA2.md`: short path map for current DATA1/DATA2 artifacts under the DATA3 subtree.
- `data2_panel_checklist_provisional.md`: DATA2 panel/table mapping status (includes locked + pending items).
- `nightly/README.md`: nightly validation folder index and purpose.
- `nightly/config/known_nonpass.csv`: allowlist of accepted nightly non-pass numeric targets.
- `nightly/nightly_ops.md`: runbook for nightly gate and prune helpers.
- `nightly/logs/nightly_test_notes.md`: rolling nightly run notes.
- `nightly/digitized_baselines/`: WebPlotDigitizer-based paper figure baseline folder (`manifest.csv`, `thresholds.csv`, `comparison_latest.csv`).

## Reproduction scaffold

Run canonical datasets:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --solver ipopt \
  --nfe 30
```

Outputs are written under:

- `UnifiedFramework/DATA3/results/reproduction/<run_id>/summaries/*.json`
- `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/unified_metrics.csv`
- `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/paper_vs_unified_side_by_side.csv`

Main-paper DATA2 table baselines currently live under:

- `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table3_side_by_side.csv`
- `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table4_side_by_side.csv`
- `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table5_side_by_side.csv`
- `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table6_side_by_side.csv`
- `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/table_baseline_extraction_notes.txt`

Notes:
- These table artifacts are `WARNING` baselines until unified-side values are populated.
- `D2-M-VARS` is a reference-only metadata integrity target (Table 1 inputs/measurements), not a numeric reproduction target.

## Storage Helpers

Use dry-run first:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive
```

Apply prune (archives old `UnifiedFramework/DATA3/figures/` first):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive \
  --apply
```

Restore one pruned run's figures:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/restore_reproduction_figures.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --run-id 20260221-004354 \
  --apply
```

If `UnifiedFramework/DATA3/figures/` already exists and you want to replace it, add `--overwrite`.

Refresh consolidated numeric target status (`PASS/FAIL/MISSING_VALUE/NOT_APPLICABLE`):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/update_target_validation_status_consolidated.py \
  --consolidated-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv \
  --tolerance-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/tolerance_eval_refresh_20260221.csv \
  --evidence-run-id 20260221-refresh
```

Nightly status gate (unexpected non-pass targets fail, known ones warn):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/nightly_validation_gate.py \
  --status-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv \
  --exceptions-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv
```

Digitized panel comparison (paper curve vs unified curve):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/compare_digitized_baselines.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --manifest /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/manifest.csv \
  --thresholds /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/thresholds.csv \
  --out-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/comparison_latest.csv
```

Current blocker for full DATA2 numeric closure:
- `run_DATA2_model_variations.py` fails in this runtime (Pyomo NL writer overflow under Python 3.12), so unified-side values for main-paper Tables 3-6 are not yet auto-populated.
