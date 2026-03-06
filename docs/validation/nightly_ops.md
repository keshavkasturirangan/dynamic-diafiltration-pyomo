# Nightly Ops

This guide describes the nightly validation gate and storage-aware pruning helpers.

## What runs nightly

Workflow: `.github/workflows/nightly-validation.yml`

Nightly steps:
1. Run status gate on `docs/validation/target_validation_status_consolidated.csv`.
2. Treat only unexpected `FAIL`/`MISSING_VALUE` targets as hard failures.
3. Keep documented known non-pass targets as warnings via `docs/validation/nightly_known_nonpass.csv`.
4. Run storage prune in dry-run mode for visibility into reclaimable space.

## Local commands

Run the nightly status gate locally:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/validation/nightly_validation_gate.py \
  --status-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/target_validation_status_consolidated.csv \
  --exceptions-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/nightly_known_nonpass.csv
```

Prune dry-run (no deletions):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive
```

Prune apply (archives figures then removes pruned run figures):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive \
  --apply
```

## How to update known non-pass targets

Edit `docs/validation/nightly_known_nonpass.csv`:
- Keep one row per currently accepted non-pass target.
- Set `expected_status` to `FAIL` or `MISSING_VALUE`.
- Set `active=true` to allow warning behavior.
- Use `review_after` to force periodic cleanup.

When a target improves to `PASS` or `NOT_APPLICABLE`, remove or deactivate that row.

## Current behavior boundaries

- `NOT_APPLICABLE` is never treated as a failure in nightly gating.
- `PASS_WITH_EXPLANATION` is currently treated as pass.
- The gate fails only on unexpected `FAIL`/`MISSING_VALUE` or status mismatches against the exception list.
