# Unified Paper-Profile Parity Checklist

This note turns the current DATA1/DATA2 paper-parity work into a concrete checklist for the unified codebase.

The goal is simple:

- `utility.py` and the legacy MATLAB / notebook workflow remain the scientific reference for published DATA1 and DATA2 results
- `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py` should expose the corresponding behavior explicitly
- parity decisions should be visible in code, tests, and validation artifacts rather than implied by scattered special cases

## Current explicit profiles

The unified runner now exposes these named paper or workflow intents:

- process-model profiles:
  - `DATA1`
  - `DATA2`
  - `DATA3`
- `DATA1_PAPER`
- `DATA2_PAPER`
- `DATA3_UNIFIED`

These are wired through:

- [`unified_codebase_runfile.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py)
- [`unified_codebase_library.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py)

## Already represented in unified code

### DATA1 paper path

These legacy DATA1 behaviors are already encoded in the unified codebase:

- explicit `paper_profile="DATA1_PAPER"` in the runner presets
- companion `_cF.csv` override for `C_F0`, matching the MATLAB DATA1 paper workflow
- DATA1 legacy-fit bounds for `Lp`, `B`, and `sigma`
- 300-point discretization in the DATA1 legacy-fit path
- grouped measurement-objective handling for legacy-style MAT fitting
- endpoint-assay tagging for scalar DATA1 retentate measurements
- DATA1 filtration fallback `C_D = 0.0` for legacy MAT files that omit dialysate concentration
- DATA1 vial boundary carryover and reset repair after simulator initialization
- run metadata now records `paper_profile`

### DATA2 paper path

These unified pieces are already in place for DATA2:

- explicit `paper_profile="DATA2_PAPER"` in the runner presets
- committed simulation validation CSV corpus for DATA2 main and SI figure/table targets
- figure-level validation against extracted paper figures
- fast DATA2 pytest coverage for committed measurement-side comparisons

## What is validated today

### Pytest coverage

Active pytests already confirm several parts of the paper-parity scaffolding:

- [`test_unified_codebase_pytest_validation.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_unified_codebase_pytest_validation.py)
- [`test_simulation_validation_figures.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py)
- [`test_utility_input_contract.py`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_utility_input_contract.py)

Current practical status:

- unified code now has an explicit process-model-profile concept in addition to paper-profile labels
- unified code recognizes paper profiles explicitly
- simulation validation baselines are committed and reusable
- extracted-paper figure matching is green
- DATA1 concentration parity still has one known `xfail`

### Validation artifacts

Useful reference outputs live under:

- [`simulation_validation_data_files/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files)
- [`results/simulation_validation/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation)
- [`results/pytest_validation/`](/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/pytest_validation)

## Remaining open gaps

### 1. DATA1 estimation-path parity is still incomplete

This is the most visible remaining gap.

What is true now:

- fixed-theta DATA1 simulation parity is close when unified and legacy are run at the same legacy parameters
- the unified DATA1 estimator still lands at different fitted `Lp`, `B`, and sometimes `sigma`
- that difference is what keeps the DATA1 Figure 2 concentration comparison in `xfail`

Implication:

- the remaining DATA1 gap is mainly in fitting workflow parity, not in the basic simulation equations alone

### 2. MATLAB paper-fit workflow is not yet first-class in unified code

The legacy DATA1 workflow uses a concrete path through:

- `run_data_analysis.m`
- `model_config.m`
- `update_model_config.m`
- `weight_config(...)`
- `fitting(...)`

The unified code mirrors part of that behavior, but not the whole operational workflow as a named, inspectable path.

Needed:

- make the DATA1 paper-fit recipe explicit enough that a user can say "run the DATA1 paper workflow" without relying on undocumented legacy assumptions

### 3. DATA1 paper-model recipe should be named directly

The DATA1 paper workflow is not just "MAT input."

It is specifically tied to the legacy paper recipe:

- `mod = '201cvmv'`
- `concpolar = true`

Needed:

- expose this recipe as a first-class DATA1 paper-model configuration in unified code and docs

### 4. DATA2 paper-fit parity is less mature than DATA1 paper-fit parity

DATA2 has good baseline and figure-validation coverage, but the unified code still needs a clearer paper-fit workflow mirror.

Needed:

- explicit mapping from DATA2 legacy scripts and notebooks into unified-paper configuration
- direct unified-fit-vs-legacy-fit validation for the key DATA2 paper figures and tables

### 5. Native unified paper-figure generation is still incomplete

The current validation flow proves that committed baseline CSVs and generated composites match the extracted paper figures.

What is still missing:

- a first-class unified figure-generation path that owns the final DATA1 and DATA2 paper-style plots directly, rather than validating against legacy-generated outputs after the fact

## Recommended next implementation order

1. Finish DATA1 estimation-path parity under `DATA1_PAPER`.
2. Keep the process-model profile as the primary scientific split and move more DATA1/DATA2 behavior into that object-oriented layer.
3. Encode the exact DATA1 paper model recipe explicitly in unified config.
4. Add direct DATA2 paper-fit workflow mapping under `DATA2_PAPER`.
5. Add native unified figure-output adapters for the paper plots.
6. Promote key paper-figure targets from known-gap status into strict numeric pytest comparisons.

## Advisor-facing summary

If you need to explain the current status quickly:

- the unified code now knows when it is trying to behave like the DATA1 or DATA2 paper workflow
- the paper-era simulation baseline corpus is committed and validated against extracted paper figures
- DATA1 still has an estimation-parity gap for concentration panels
- DATA2 still needs the same depth of explicit paper-fit workflow mirroring that DATA1 now has started
