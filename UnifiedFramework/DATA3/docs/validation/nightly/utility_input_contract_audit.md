# `utility.py` Input Contract Audit

This note documents what [`utility.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py) actually expects from `data_stru`, based on code inspection and representative MAT files.

The goal is not to add fallbacks casually. The goal is to understand:

- what fields are truly required to run the DATA1 and DATA2 pipelines
- which requirements are mode-dependent
- where current legacy files do and do not satisfy those requirements

## Scope

This audit covers:

- [`solve_model(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
- [`solve_model_B_fix(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
- [`model_construct_inter(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
- plotting/interpolation helpers that consume the same `data_stru`

## Core `data_stru` Shape

`utility.py` expects a dict with:

- `data_stru['data_config']`
- `data_stru['data_raw']`

### `data_config` fields referenced in `utility.py`

Always referenced in the main solve path:

- `Am`
- `B0`
- `C_F0`
- `Lp0`
- `M_F0`
- `Temp`
- `delP`
- `n`
- `namec`
- `ni`
- `rho`
- `sigma0`

Referenced in gating / objective / DATA2-style branches:

- `n_v0`
- `n_extra`
- `M_O`
- `C_D`
- `n_h`
- `n_A`

Also present in many MAT files but not central to this solve path:

- `theta0`
- `nc`
- `nr`

### `data_raw` fields referenced in `utility.py`

- `time`
- `mass`
- `cV_avg`
- `cF_exp`

Other raw keys may exist, like `number` or `nr`, but `utility.py` does not depend on them for model construction.

## Mode-Specific Requirements

### 1. DATA1 mode

Legacy branch mapping:
- [`model_construct_inter(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
  with `mode == 'DATA'`

What the code actually uses in DATA1 mode:

- all always-required core fields listed above
- `n_v0` in mass/retentate objective gating and in plotting/interpolation
- `n_extra` in the objective for permeate/vial counting
- `C_D` in the DATA1 retentate ODE:
  - `dcF/dt = Am * rho / M_F0 * (cD * Jw - Js)`
- `M_O` is read early in model construction, but not meaningfully used in the DATA branch
- `n_h` and `n_A` are not used in DATA1 mode

What this means:

- `C_D` is structurally required by the DATA1 equations as written
- `n_v0` and `n_extra` are structurally required by the DATA1 objective/gating as written
- `M_O` is not logically required for DATA1 mode, even though the function currently reads it unconditionally near the top

### 2. DATA2 mode

Legacy branch mapping:
- [`model_construct_inter(...)`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py)
  with `mode != 'DATA'`, especially `Lag` / `Overflow`

What the code uses additionally in DATA2 mode:

- `M_O`
- `C_D`
- `n_v0`
- `n_extra`
- `n_h`
- `n_A`

Why:

- `M_O` is used in overflow/lag mass balance and initialization
- `C_D` is used in retentate dynamics
- `n_v0`, `n_extra`, `n_h`, `n_A` control startup/holdup/overflow vial logic

What this means:

- for DATA2 mode, these fields are not optional
- they are part of the actual model definition, not just plotting metadata

## Representative File Audit

### DATA1 `501.1`

File:
- [`data_stru-dataset501.1.mat`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset501.1.mat)

Observed `data_config` keys:

- `Am`
- `B0`
- `C_F0`
- `Lp0`
- `M_F0`
- `Temp`
- `delP`
- `n`
- `namec`
- `nc`
- `ni`
- `nr`
- `rho`
- `sigma0`
- `theta0`

Missing relative to the full `utility.py` contract:

- `C_D`
- `M_O`
- `n_v0`
- `n_extra`
- `n_h`
- `n_A`

Observed first-vial raw field types:

- `cF_exp`: `float`
- `cV_avg`: `float`
- `mass`: `list`
- `time`: `list`

Implication:

- this file does not satisfy the full contract that `utility.py` assumes
- in particular, `cF_exp` is scalar here, while several `utility.py` paths assume iterable retentate measurements

### DATA1 `511.12`

File:
- [`data_stru-dataset511.12.mat`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_matlab/data_library/data_stru-dataset511.12.mat)

Observed `data_config` keys:

- all of the `501.1` keys above
- plus `C_D`
- plus `M_O`

Still missing:

- `n_v0`
- `n_extra`
- `n_h`
- `n_A`

Observed first-vial raw field types:

- `cF_exp`: `list`
- `cV_avg`: `float`
- `mass`: `list`
- `time`: `list`

Implication:

- this file is closer to what `utility.py` expects, but still does not satisfy the full gating contract

### DATA2 `270511.121`

File:
- [`data_stru-dataset270511.121.mat`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/data_library/data_stru-dataset270511.121.mat)

Observed `data_config` keys:

- all core keys
- `C_D`
- `M_O`
- `n_extra`
- `n_v0`

Observed first-vial raw field types:

- `cF_exp`: `list`
- `cV_avg`: `list`
- `mass`: `list`
- `time`: `list`

Implication:

- this shape is much closer to the full contract assumed by `utility.py`
- this explains why the DATA2 staged workflow is more naturally compatible with the legacy code assumptions

## Main Findings

### Finding 1: `utility.py` encodes more than one input contract

It is not one universal schema.

There is a smaller DATA1 contract and a larger DATA2 contract.

### Finding 2: `utility.py` currently reads some DATA2 fields too early

Even before mode-specific logic is applied, `model_construct_inter(...)` reads:

- `M_O`
- `C_D`
- `n_v0`

That makes the code stricter than the DATA branch really needs.

### Finding 3: DATA1 MAT files are heterogeneous

`501.1` and `511.12` do not even have exactly the same config shape:

- `501.1` lacks both `C_D` and `M_O`
- `511.12` includes them

So “the DATA1 contract” is not fully uniform across files.

### Finding 4: scalar `cF_exp` is a real legacy shape

For `501.1`, `cF_exp` is scalar in the first vial.

That means any path that assumes `for item in cF_exp` or `len(cF_exp)` without normalization will fail on legitimate DATA1 inputs.

## Practical Guidance

### Safe conclusion

Before changing behavior, the code should distinguish:

- fields required for DATA1 mode
- fields required only for DATA2 mode
- fields that may be scalar or vector depending on dataset family

### Recommended coding rule

Do not treat missing fields as arbitrary “fallback” opportunities.

Instead:

- only require fields that the active mode actually needs
- normalize scalar-vs-vector measurement shapes explicitly where legacy files use both forms
- preserve strict requirements for DATA2 mode where those fields affect real model behavior

### Implication for current work

For the current DATA1 concentration audit:

- using `utility.py` as the guide is still correct
- but we should interpret it as a legacy implementation that spans DATA1 and DATA2 pipelines, not a single clean schema
- any hardening should be mode-aware and justified by the active equations/objective, not by convenience
