# Unified Codebase Validation Development Plan

## Objective
Validate the unified framework by reproducing required DATA1 and DATA2 published results using:
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/experiment_dataload_OOP_v24.py`
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/experiment_load_v24.py`
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/conductivity_paper.py`

## Scope and Rules
- Reproduction target: all required figures/tables listed for DATA1 main/SI and DATA2 main/SI.
- Solver policy: use `ipopt` (no `ef_ipopt`).
- Statistical equivalence is acceptable.
- Tolerances:
  - Parameters: <= 5% relative error (<= 7% with explanation)
  - Objective/WLS/SSE: <= 5% relative difference
  - Curves/plots: NRMSE <= 5% of observed range
  - Table values: absolute tolerance based on reported precision/significant digits
- Plot outputs must match paper styling.
- Nightly tests: fail on tolerance miss.

## Dataset Priorities
- Stage A (first concrete benchmark): `data_stru-dataset511.12.mat`
- Stage C core DATA2 first files: `270511.123`, `270611.123`
- Canonical file mapping sourced from:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA1_model_demo.ipynb`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_model_demo.ipynb`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_visualization.ipynb`

## Workstreams

### WS1: Reproduction Matrix and References
- Maintain required target list and status in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/paper_target_matrix.md`
- Maintain side-by-side paper reference values in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/docs/validation/paper_reference_values_template.csv`

### WS2: Unified Execution + Artifact Generation
- Primary script:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/scripts/reproduce/reproduce_data1_data2.py`
- Responsibilities:
  - Load MAT datasets through unified loader
  - Fit parameters via ParmEst (`ipopt`)
  - Generate paper-style figures
  - Generate unified metrics and side-by-side CSV tables

### WS3: Plotting Integration in Unified Module
- Add/expand plotting utilities in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataLoader/UnifiedCode/experiment_dataload_OOP_v24.py`
- Ensure figure generation is reusable by scripts and future package API.

### WS4: Regression Test Automation
- Nightly regression test:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/tests/regression/test_data1_data2_nightly.py`
- Markers and test config in `pytest.ini`
- CI mode: slow/nightly for DATA2-heavy workloads

### WS5: Repo Cleanup Proposal (No Moves Yet)
- Deliver a concrete rename/move proposal first.
- Execute moves only after approval.

## Stage Plan

### Stage A (In Progress)
Goal:
- DATA1 main: Fig.2 (all), Table 1, Fig.3, Table 2

Deliverables:
- Generated figures
- Side-by-side CSV for tables/metrics
- Tolerance status report

Exit Criteria:
- All Stage A outputs generated from unified code path
- Numerical comparison table produced
- Any misses documented with explanation

### Stage B
Goal:
- Remaining DATA1 main + DATA1 SI figures/tables

Deliverables:
- Complete DATA1 artifact set
- Updated matrix status and residual gap notes

Exit Criteria:
- DATA1 targets fully covered and tolerance-assessed

### Stage C
Goal:
- DATA2 core first (`270511.123`, `270611.123`), then full DATA2 main/SI

Deliverables:
- Manipulated/measured variable verification report
- Required DATA2 tables/figures
- Side-by-side numeric comparisons

Exit Criteria:
- DATA2 core cases pass targets; full set generated and assessed

### Stage D
Goal:
- Stabilize nightly regression runs and reporting

Deliverables:
- Deterministic nightly job outputs
- Failure summaries tied to target IDs and tolerances

Exit Criteria:
- Nightly test fails only on true tolerance regressions

### Stage E
Goal:
- Package-readiness prep

Deliverables:
- Final repo reorg proposal and accepted structure
- Mapping from scripts to future package modules/CLI

Exit Criteria:
- Approved transition plan to Python package layout

## File/Artifact Conventions
- Run outputs: `results/reproduction/<run_id>/`
  - `figures/`
  - `tables/`
  - `summaries/`
- Comparison table must include:
  - target ID
  - paper value
  - unified value
  - absolute/relative error
  - pass/fail status
  - explanation field when using 5–7% relaxed parameter window

## Immediate Next Actions
1. Populate paper reference values for Stage A rows in `paper_reference_values_template.csv`.
2. Run Stage A script with `ipopt` and generate artifact bundle.
3. Validate Stage A against tolerance thresholds and update status matrix.
4. Resolve any Stage A mismatches (model options, plotting parity, data mapping).
5. Begin Stage B expansion once Stage A is accepted.

## Open Items to Track
- Covariance behavior when `k_aug` is unavailable in environment (documented limitation path).
- Any solver option tuning needed for stable covariance/fit across datasets.
- Final plotting parity details (fonts/line styles/layout per paper).
- Repo cleanup rename/move proposal document (pending).
