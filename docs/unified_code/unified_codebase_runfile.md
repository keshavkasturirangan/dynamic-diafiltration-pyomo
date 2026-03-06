# `unified_codebase_runfile.py`

## Purpose

`unified_codebase_runfile.py` is the user-facing runner script.  
It now supports profile-based switching for:

- `DATA1` reproduction,
- `DATA2` reproduction,
- `DATA3` experiments (custom file/sheet workflow).

Heavy lifting remains delegated to `unified_codebase_library.py`.

## Execution flow

1. Select `--profile DATA1|DATA2|DATA3`.
2. Resolve default file/selector/model-options from profile.
3. Apply optional overrides (`--file-path`, `--selector`, `--nfe`, flags).
4. Build `UnifiedPipelineConfigV24`.
5. Execute `run_unified_pipeline_v24`.
6. Print run summary (load status, theta labels, estimation/DoE output).

## Profiles

- `DATA1`
  - file: `DATA1_matlab/data_library/data_stru-dataset511.12.mat`
  - defaults: `mode=DATA`, `b_form=single`, estimation run mode
- `DATA2`
  - files:
    - `DATA1_matlab/data_library/data_stru-dataset270511.123.mat`
    - `DATA1_matlab/data_library/data_stru-dataset270611.123.mat`
  - defaults: `mode=Lag`, `b_form=convection`, `fix_sigma_in_estimation=True`
- `DATA3`
  - file: `UnifiedFramework/ExperimentalDataFiles/NF270_MC2.xlsx` (override recommended)
  - default selector/specs enabled for quick XLSX experiments
  - default `run_doe=True` (can be overridden with `--no-run-doe`)

## Key CLI options

- `--profile DATA1|DATA2|DATA3`
- `--file-path <path>`
- `--selector <xlsx_sheet_name>`
- `--nfe <int>`
- `--run-doe`
- `--no-run-doe`
- `--calc-cov`
- `--convert-to-concentration`
- `--plot`
- `--solver ipopt`

## Typical output

- load status and validation issues,
- theta names (unknown parameter set),
- number of labeled measurements and inputs,
- estimated objective and theta,
- covariance (if converged),
- FIM and D-optimality metric.

## Example commands

```bash
# DATA1 reproduction baseline
python UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/unified_codebase_runfile.py \
  --profile DATA1 --nfe 30

# DATA2 reproduction baseline with DoE
python UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/unified_codebase_runfile.py \
  --profile DATA2 --run-doe --nfe 30

# DATA3 experiment from custom XLSX
python UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/unified_codebase_runfile.py \
  --profile DATA3 \
  --file-path /abs/path/to/experiment.xlsx \
  --selector "SheetName" \
  --nfe 30
```

## Operational guidance

- Start with `nfe=30` for stable estimation/covariance.
- Increase `nfe` only when needed for accuracy.
- Keep `solver="ipopt"` for ParmEst.
- If DoE struggles, tune NLP settings before increasing model complexity.

## Packaging/refactor notes

- This file should become a thin CLI entrypoint (`console_scripts`) in the package.
- Config values should move into a validated config object or YAML schema.
- Logging should replace `print` for reproducible runs.
