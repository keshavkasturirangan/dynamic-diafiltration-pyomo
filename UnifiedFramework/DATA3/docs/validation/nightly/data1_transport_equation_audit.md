# DATA1 Transport Equation Audit

Date: 2026-03-30

Goal:
- Compare the DATA1 transport-model behavior in `utility.py` and the unified code.
- Allow one intentional difference only: the unified code may use MSA instead of the legacy conductivity-to-concentration correlation.

Scope:
- `DATA1` MAT / `ExperimentMode.DATA`
- single-`B` estimation path
- focus on Fig. 2 concentration mismatch

## Main Result

The remaining DATA1 concentration mismatch does not appear to come from the core DATA-mode transport equations. After auditing the equations and state-linking rules, and after the later fixed-theta state-trace comparison, the unified DATA-mode model appears to be aligned with the legacy `utility.py` formulation closely enough for fixed-theta simulation parity.

The one confirmed structural bug that was found and fixed was:
- scalar MAT `cF_exp` assays were being broadcast and counted as full time-series residuals in the unified objective
- legacy `utility.py` treats those scalar assays as endpoint-style measurements

Fixing that changed dataset `501.1` from:
- unified `sigma ~= 0.001`

to:
- unified `sigma ~= 0.999`

which is much closer to the legacy MATLAB fit:
- legacy `sigma ~= 1.0`

## Equation Audit Summary

For `ExperimentMode.DATA`, the following were checked against `utility.py`:

- `ode_cF`
  - legacy: `dcF/dt = Am * rho / M_F0 * (cD * Jw - Js)`
  - unified: same

- `ode_cH`
  - legacy DATA branch uses constant `mH`
  - unified DATA branch uses constant `MH_ML`
  - behavior is aligned

- `ode_mV`
  - legacy and unified both use `dmV/dt = Jw * Am * rho`

- `ode_cVmV`
  - legacy DATA branch uses `Jw * cH * Am * rho`
  - unified DATA branch uses the same

- `eqn_cIn`
  - legacy and unified both use the same exponential film relation in DATA mode

- `eqn_Jw`
  - legacy and unified both use
    - `Jw * 36000 = Lp * (delP - (cIn - cH) * ni * sigma * R * T)`

- `eqn_Js`
  - legacy and unified both use
    - `Js * 10000 = B * (cIn - cH)`
    for single-`B`

- `eqn_cV`
  - legacy and unified both use
    - `mV * cV = cVmV`

- inter-vial linking
  - `cF_linking`, `cH_linking`, `mV_linking`, `cVmV_linking`
  - DATA-mode behavior matches the legacy reset/continuity pattern used for Fig. 2

- initial conditions
  - `cH[1,0] = first_vial cV_avg * 0.8`
  - `cF[1,0] = first non-NaN cF_exp`
  - `mV[1,0] = 1e-6`
  - `cVmV[1,0] = 1e-12`
  - unified matches legacy here too

## Confirmed Fix

File:
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`

Change:
- add `retentate_signal_endpoint_assay`
- when MAT `cF_exp` is scalar/length-1, store it compatibly but count it once at the vial endpoint in the legacy grouped objective

Why this matters:
- legacy `utility.py` counts scalar `cF_exp` once per vial
- the previous unified path effectively counted it many times
- that distorted DATA1 concentration weighting

## What Was Tested

### 1. Parameter drift after objective fix

Artifacts:
- `UnifiedFramework/DATA3/results/pytest_validation/data1_concentration_gap_summary.csv`
- `UnifiedFramework/DATA3/results/pytest_validation/data1_concentration_gap_endpoints.csv`

Observation:
- `501.1` now recovers `sigma ~= 0.999`
- but `Lp` and `B` still differ from legacy MATLAB
- concentration mismatch remains

### 2. Discretization sensitivity

Unified DATA1 `501.1` fit was rerun with:
- `nfe = 30`
- `nfe = 60`
- `nfe = 120`
- `nfe = 300`

Result:
- all runs converged to essentially the same parameter set
- increasing `nfe` did not move the solution toward the legacy MATLAB fit

This means the concentration gap is not explained by coarse discretization in the current unified validation path.

## Current Best Interpretation

After fixing the assay-weighting bug, the remaining mismatch is likely not due to:
- the simulation validation CSVs
- the DATA-mode transport equations listed above
- the MAT initial guess values
- the DATA-mode initialization/linking rules
- using `nfe=30` instead of `300`

The remaining gap is therefore most likely caused by one of:
- a legacy-vs-unified estimation workflow difference
- a solver/optimization-path difference that is not captured by `nfe`
- a subtle difference in how the unified model is initialized before optimization versus legacy MATLAB

## Important Constraint

The intended modeling difference should be:
- unified code may use MSA instead of the legacy conductivity-to-concentration correlation

For DATA1 MAT concentration fits, that difference should not be active, because the MAT workflow is already using concentration-like values from the dataset rather than converting conductivity online.

## Recommended Next Step

Focus on the estimation path, not the fixed-theta simulation path.

Most useful follow-up checks:

- compare the legacy and unified fitting objectives term-by-term
- compare fitted `Lp`, `B`, and `sigma` under the same DATA1 paper recipe
- treat new fixed-theta simulation mismatches as suspicious, because the current
  state-trace evidence indicates that fixed-theta parity is already essentially
  achieved
