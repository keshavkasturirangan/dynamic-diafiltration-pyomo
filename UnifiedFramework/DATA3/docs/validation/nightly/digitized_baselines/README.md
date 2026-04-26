# Digitized Figure Baselines (WebPlotDigitizer)

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


This folder stores numeric baselines extracted from published figure panels (e.g., via https://automeris.io/).

## Purpose

Use digitized panel data as paper-reference truth when table values are unavailable.
Then compare unified-pipeline outputs against those digitized baselines with reproducible metrics.

## Files

- `manifest.csv`: one row per figure/panel comparison target.
- `thresholds.csv`: target-specific pass/fail tolerances.
- `paper/`: digitized CSV exports from published figures.
- `unified/`: unified-pipeline CSV exports aligned to same x/y definitions.

## Minimal workflow

1. Digitize paper figure/panel and save CSV in `paper/`.
2. Export matching unified curve as CSV in `unified/`.
3. Add row in `manifest.csv` with file paths and column mapping.
4. Set tolerances in `thresholds.csv`.
5. Run:

```bash
python /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/validation/nightly/compare_digitized_baselines.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --manifest /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/manifest.csv \
  --thresholds /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/thresholds.csv \
  --out-csv /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/comparison_latest.csv
```

## Notes

- Keep axes/units consistent (`x` and `y` definitions must match between paper and unified).
- Use `x_min`/`x_max` in `manifest.csv` to compare only the trusted panel range.
- `active=false` rows are ignored.
