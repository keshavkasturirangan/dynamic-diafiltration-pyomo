# Validation Docs

This folder tracks DATA1/DATA2 published-result reproduction for unified-code validation.

## Files

- `paper_target_matrix.md`: required figure/table targets and stage plan.
- `paper_reference_values_template.csv`: side-by-side paper/unified reference template.
- `target_validation_status_consolidated.csv`: target-level status (`PASS/FAIL/MISSING_VALUE/NOT_APPLICABLE`).
- `target_notebook_source_map.csv`: target-to-artifact evidence mapping.
- `target_pdf_page_index.csv`: target-to-paper-page evidence index.
- `panel_validation_both_DATA1_DATA2.md`: current high-level checkpoint.
- `data2_panel_checklist_provisional.md`: DATA2 panel/table mapping status (includes locked + pending items).
- `nightly_known_nonpass.csv`: allowlist of accepted nightly non-pass numeric targets.
- `nightly_ops.md`: runbook for nightly gate and prune helpers.

## Reproduction scaffold

Run canonical datasets:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/reproduce_data1_data2.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --solver ipopt \
  --nfe 30
```

Outputs are written under:

- `results/reproduction/<run_id>/summaries/*.json`
- `results/reproduction/<run_id>/tables/unified_metrics.csv`
- `results/reproduction/<run_id>/tables/paper_vs_unified_side_by_side.csv`

Main-paper DATA2 table baselines currently live under:

- `results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table3_side_by_side.csv`
- `results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table4_side_by_side.csv`
- `results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table5_side_by_side.csv`
- `results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table6_side_by_side.csv`
- `results/reproduction/20260306-mainpaper-table-baseline/tables/table_baseline_extraction_notes.txt`

Notes:
- These table artifacts are `WARNING` baselines until unified-side values are populated.
- `D2-M-VARS` is a reference-only metadata integrity target (Table 1 inputs/measurements), not a numeric reproduction target.

## Storage Helpers

Use dry-run first:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive
```

Apply prune (archives old `figures/` first):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive \
  --apply
```

Restore one pruned run's figures:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/restore_reproduction_figures.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --run-id 20260221-004354 \
  --apply
```

If `figures/` already exists and you want to replace it, add `--overwrite`.

Refresh consolidated numeric target status (`PASS/FAIL/MISSING_VALUE/NOT_APPLICABLE`):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/validation/update_target_validation_status_consolidated.py \
  --consolidated-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/target_validation_status_consolidated.csv \
  --tolerance-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/tolerance_eval_refresh_20260221.csv \
  --evidence-run-id 20260221-refresh
```

Nightly status gate (unexpected non-pass targets fail, known ones warn):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/validation/nightly_validation_gate.py \
  --status-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/target_validation_status_consolidated.csv \
  --exceptions-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/nightly_known_nonpass.csv
```

Current blocker for full DATA2 numeric closure:
- `run_DATA2_model_variations.py` fails in this runtime (Pyomo NL writer overflow under Python 3.12), so unified-side values for main-paper Tables 3-6 are not yet auto-populated.
