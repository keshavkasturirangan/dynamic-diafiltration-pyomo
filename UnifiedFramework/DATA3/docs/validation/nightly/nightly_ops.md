# Nightly Ops

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


This guide describes the nightly validation gate and storage-aware pruning helpers.

## What runs nightly

Workflow: `.github/workflows/nightly-validation.yml`

Nightly steps:
1. Run status gate on `UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv`.
2. Treat only unexpected `FAIL`/`MISSING_VALUE` targets as hard failures.
3. Keep documented known non-pass targets as warnings via `UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv`.
4. Run digitized panel comparison report (`UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/comparison_latest.csv`).
5. Run storage prune in dry-run mode for visibility into reclaimable space.

## Local commands

Run the nightly status gate locally:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/nightly_validation_gate.py \
  --status-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv \
  --exceptions-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv
```

Run digitized paper-vs-unified figure comparison:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/compare_digitized_baselines.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --manifest /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/manifest.csv \
  --thresholds /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/thresholds.csv \
  --out-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/comparison_latest.csv
```

Append an entry to the running nightly notes document:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/append_nightly_note.py \
  --run-id 20260306-nightly-local \
  --context "Local post-setup validation" \
  --gate-status PASS \
  --gate-summary "No unexpected FAIL/MISSING_VALUE; 7 known warnings allowed." \
  --pytest-status TIMEOUT \
  --action "Keep gate as nightly health signal; optimize full pytest runtime."
```

Prune dry-run (no deletions):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive
```

Prune apply (archives figures then removes pruned run figures):

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive \
  --apply
```

## How to update known non-pass targets

Edit `UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv`:
- Keep one row per currently accepted non-pass target.
- Set `expected_status` to `FAIL` or `MISSING_VALUE`.
- Set `active=true` to allow warning behavior.
- Use `review_after` to force periodic cleanup.

When a target improves to `PASS` or `NOT_APPLICABLE`, remove or deactivate that row.

## Current behavior boundaries

- `NOT_APPLICABLE` is never treated as a failure in nightly gating.
- `PASS_WITH_EXPLANATION` is currently treated as pass.
- The gate fails only on unexpected `FAIL`/`MISSING_VALUE` or status mismatches against the exception list.

## Running notes document

Use this file as the rolling nightly checkpoint log:

- `UnifiedFramework/DATA3/docs/validation/nightly/logs/nightly_test_notes.md`

## Digitized baseline activation

To enable a panel target:
1. Export digitized paper points (`paper/*.csv`) from WebPlotDigitizer.
2. Export matching unified outputs (`unified/*.csv`) on same axes/units.
3. Set the row to `active=true` in `UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/manifest.csv`.
4. Tune tolerances in `UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/thresholds.csv`.

When you want strict nightly enforcement, add flags to the compare command:
- `--fail-on-fail`
- `--fail-on-missing`
