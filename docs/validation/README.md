# Validation Docs

This folder tracks DATA1/DATA2 published-result reproduction for unified-code validation.

## Files

- `paper_target_matrix.md`: required figure/table targets and stage plan.
- `paper_reference_values_template.csv`: side-by-side paper/unified reference template.

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
