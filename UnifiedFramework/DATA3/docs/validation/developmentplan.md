# Unified Codebase Validation Development Plan

> **Canonical code path (source of truth):** `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/`  
> Runner entrypoint: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`  
> Core pipeline: `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`  
> Legacy root-level scripts (for example `utility.py`, `run_*.py`) are compatibility/reference paths, not primary development paths.


## Objective
Validate the unified framework by reproducing required DATA1 and DATA2 published results using:
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/conductivity_paper.py`

## Validation Update (2026-03-30)
- The committed simulation-side CSV baseline corpus is now:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/`
- These CSVs are the one-time published-paper baseline export assembled from:
  - `utility.py`
  - `DATA1_model_demo.ipynb`
  - `DATA2_model_demo.ipynb`
  - `DATA2_visualization.ipynb`
  - `run_DATA2_model_variations.py`
- The current figure-validation report comparing mapped/generated legacy figure composites against extracted paper figures is:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/results/simulation_validation/simulation_validation_figure_report.csv`
- Latest status from that report: `23/23 PASS` for mapped DATA1 main/SI and DATA2 main/SI figure targets.
- Dedicated figure-validation pytest:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py`
- User-facing conductivity default in `unified_codebase_runfile.py` is now `msa`; `variant_shedlovsky` must be requested explicitly with `--conductivity-model variant_shedlovsky`.

## Validation Baseline Corpus (2026-03-05)

Use these sources as the reference corpus for validation and parity decisions.

1. Current unified code + docs:
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/conductivity_paper.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/unified_code/conductivity_paper.md`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/unified_code/unified_codebase_library.md`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/unified_code/unified_codebase_runfile.md`
2. Last v24 compatibility wrappers and pseudocode lineage:
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/archives/experiment_dataload_OOP_v24.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/experiment_load_v24.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/Pseudocodes/DataLoader/versions/pseudocode_v2.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/archives/processmodel_pseudocode_v8_updated_conductivity.py`
3. Legacy reproduction stack (DATA1/DATA2 parity):
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/utility.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_model_demo.ipynb`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/DATA2_visualization.ipynb`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/run_cross_verification.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/run_DATA2_model_variations.py`
   - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/run_pre_B_dependence.py`
4. Published references:
   - `UnifiedFramework/DATA3/published_works/` (paper UnifiedFramework/DATA3/figures/tables used for acceptance checks)
5. Digitization fallback when panel/table source data are not available:
   - `https://automeris.io/` (WebPlotDigitizer) with source path + extraction notes recorded per target.

## Scope and Rules
- Reproduction target: all required UnifiedFramework/DATA3/figures/tables listed for DATA1 main/SI and DATA2 main/SI.
- Solver policy: use `ipopt` (no `ef_ipopt`).
- Statistical equivalence is acceptable.
- Tolerances:
  - Parameters: <= 5% relative error (<= 7% with explanation)
  - Objective/WLS/SSE: <= 5% relative difference
  - Curves/plots: NRMSE <= 5% of observed range
  - Table values: absolute tolerance based on reported precision/significant digits
- Plot outputs must match paper styling.
- Nightly tests: fail on tolerance miss.

## Execution Reset (2026-02-21)
- Problem: work expanded into structural refactors without clear stop conditions.
- Reset principle: reproduction-first, refactor-second.
- Hard rule: no additional architecture work unless it directly unblocks a required paper target.

### What “Done” means right now
1. DATA1 main + SI targets are locked to concrete artifacts and pass/fail tolerance status is computed.
2. DATA2 core (`270511.123`, `270611.123`) has source-backed baselines and pass/fail tolerance status.
3. Unified run path supports both MAT and XLSX from one external API (`run_unified_pipeline_v24`), with compatibility wrappers retained.
4. Nightly regression scaffold is runnable and tied to target IDs.
5. Lightweight nightly validation is exercised through committed `pytest` checks against the consolidated status CSV and extracted paper-linked figure artifacts.

### Explicitly deferred until after above is done
- Further module splits/renames.
- Packaging polish not required for immediate reproduction.
- Non-essential API redesign.

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
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md`
- Maintain side-by-side paper reference values in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv`

### WS2: Unified Execution + Artifact Generation
- Primary script:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`
- Responsibilities:
  - Load MAT datasets through unified loader
  - Fit parameters via ParmEst (`ipopt`)
  - Generate paper-style figures
  - Generate unified metrics and side-by-side CSV tables

### WS3: Plotting Integration in Unified Module
- Add/expand plotting utilities in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
- Ensure figure generation is reusable by scripts and future package API.

### WS4: Regression Test Automation
- Nightly regression test:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/tests/regression/test_data1_data2_nightly.py`
- Lightweight nightly validation pytests:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_nightly_validation_assets.py`
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/tests/regression/test_simulation_validation_figures.py`
- Markers and test config in `pytest.ini`
- CI mode: slow/nightly for DATA2-heavy workloads
- Scope of lightweight nightly pytests:
  - validate `target_validation_status_consolidated.csv` through the allowlisted nightly gate policy,
  - verify mapped figure targets still resolve to committed generated artifacts and extracted paper pages under `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/`,
  - rerun figure-composite comparisons against extracted published-paper figures using `simulation_validation_data_files` as the committed simulation-side baseline corpus,
  - keep the solver-heavy full reproduction test separate until runtime is acceptable for scheduled CI.

### WS5: Repo Cleanup Proposal (No Moves Yet)
- Deliver a concrete rename/move proposal first.
- Execute moves only after approval.

### WS6: Profile-Switching UX for Unified Runner
- Maintain a simple user entrypoint to switch between:
  - DATA1 reproduction
  - DATA2 reproduction
  - DATA3 experiment workflows
- Current implementation lives in:
  - `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
- User contract:
  - `--profile DATA1|DATA2|DATA3`
  - optional `--file-path` / `--selector` overrides
  - optional stage toggles (`--run-doe`, `--calc-cov`)

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
- Remaining DATA1 main + DATA1 SI UnifiedFramework/DATA3/figures/tables

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
- Stabilize nightly `pytest` validation runs and reporting

Deliverables:
- Deterministic nightly `pytest` outputs
- Failure summaries tied to target IDs and tolerances
- Paper-page and artifact existence checks tied to locked/page-matched figure targets

Exit Criteria:
- Lightweight nightly `pytest` suite fails only on true unexpected regressions or broken paper/artifact mappings
- Full solver-heavy reproduction `pytest` remains available but is not required in scheduled CI until runtime is reduced

### Stage E
Goal:
- Package-readiness prep

Deliverables:
- Final repo reorg proposal and accepted structure
- Mapping from scripts to future package modules/CLI

Exit Criteria:
- Approved transition plan to Python package layout

## File/Artifact Conventions
- Run outputs: `UnifiedFramework/DATA3/results/reproduction/<run_id>/`
  - `UnifiedFramework/DATA3/figures/`
  - `tables/`
  - `summaries/`
- Comparison table must include:
  - target ID
  - paper value
  - unified value
  - absolute/relative error
  - pass/fail status
  - explanation field when using 5–7% relaxed parameter window

## Continuity Breadcrumbs (Storage-Aware)

Purpose:
- Preserve only the minimum context needed to resume quickly when disk space is constrained.
- Avoid re-scanning large historical runs before coding or validation work.

Primary breadcrumb index:
- `/Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo/UnifiedFramework/DATA3/docs/validation/continuity_ledger.csv`

Nightly validation assets (as of 2026-03-09):
- Docs root: `UnifiedFramework/DATA3/docs/validation/nightly/`
- Nightly runbook: `UnifiedFramework/DATA3/docs/validation/nightly/nightly_ops.md`
- Known non-pass allowlist: `UnifiedFramework/DATA3/docs/validation/nightly/config/known_nonpass.csv`
- Running nightly notes: `UnifiedFramework/DATA3/docs/validation/nightly/logs/nightly_test_notes.md`
- Digitized baselines: `UnifiedFramework/DATA3/docs/validation/nightly/digitized_baselines/`
- Nightly scripts: `UnifiedFramework/DATA3/scripts/validation/nightly/`

Minimum keep-set (do not delete):
- `UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md`
- `UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv`
- `UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv`
- `UnifiedFramework/DATA3/docs/validation/tolerance_eval_refresh_20260221.csv` (or latest equivalent)
- `UnifiedFramework/DATA3/docs/validation/data1_panel_checklist.md`
- For canonical reproduction run(s):
  - `UnifiedFramework/DATA3/results/reproduction/<run_id>/run_metadata.json`
  - `UnifiedFramework/DATA3/results/reproduction/<run_id>/tables/*.csv`
  - `UnifiedFramework/DATA3/results/reproduction/<run_id>/summaries/*.json`
- For locked DATA1 notebook mapping run:
  - `UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/figures/data1_notebook_regen/manifest_data1_notebook_regen.json`

Space-management policy:
- Keep only the latest 3 full reproduction run folders with `UnifiedFramework/DATA3/figures/` preserved.
- For older runs, keep `run_metadata.json`, `tables/`, and `summaries/`; archive or remove `UnifiedFramework/DATA3/figures/`.
- Do not archive/delete any run explicitly referenced by:
  - `target_validation_status_consolidated.csv`
  - `continuity_ledger.csv`
  - `data1_panel_checklist.md`
- Use the helper:
  - `UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py` (dry-run by default)

Resume-fast command set:
```bash
# Canonical unified rerun (DATA1 + DATA2 core)
python UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --solver ipopt \
  --nfe 30

# Numeric-only refresh (skip stage figure regeneration for speed)
python UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --solver ipopt \
  --nfe 30 \
  --skip-stage-artifacts

# Deterministic DATA1 notebook artifact regeneration
python UnifiedFramework/DATA3/scripts/reproduce/regenerate_data1_notebook_artifacts.py --run-id 20260221-data1-notebook

# Storage-aware prune preview (no deletion)
python UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive

# Apply storage-aware prune
python UnifiedFramework/DATA3/scripts/reproduce/prune_reproduction_artifacts.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --keep-full 3 \
  --mode archive \
  --apply

# Restore a pruned run's figures from archive
python UnifiedFramework/DATA3/scripts/reproduce/restore_reproduction_figures.py \
  --repo-root /Users/kkasturi/GitHub/keshav-dev-dynamic-diafiltration-pyomo \
  --run-id 20260221-004354 \
  --apply

# Recompute consolidated numeric status with NOT_APPLICABLE semantics
python UnifiedFramework/DATA3/scripts/validation/update_target_validation_status_consolidated.py \
  --consolidated-csv UnifiedFramework/DATA3/docs/validation/target_validation_status_consolidated.csv \
  --tolerance-csv UnifiedFramework/DATA3/docs/validation/tolerance_eval_refresh_20260221.csv \
  --evidence-run-id 20260221-refresh
```

Handoff rule:
- After each materially relevant run, update `UnifiedFramework/DATA3/docs/validation/continuity_ledger.csv` with:
  - run ID,
  - purpose,
  - retained evidence paths,
  - next action,
  - blockers (if any).
- If a run is archived/pruned, keep a pointer to the archive artifact path in the ledger.

## Execution Tracker

| Step ID | Task | Status | Date | Evidence |
|---|---|---|---|---|
| EX-001 | Create validation matrix + scaffold script | Completed | 2026-02-18 | `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`, `UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md`, `UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv` |
| EX-002 | Migrate active workflow to ipopt-only user-facing solver policy | Completed | 2026-02-18 | Commit `4317ce6`; updated `unified_codebase_library.py`, `unified_codebase_runfile.py`, regression/docs |
| EX-003 | Execute Stage A benchmark on `DATA1_511.12` | Completed | 2026-02-18 | Run IDs `20260218-182946`, `20260218-184013`; Stage A UnifiedFramework/DATA3/figures/tables generated |
| EX-004 | Execute Stage C core benchmarks (`270511.123`, `270611.123`) | Completed | 2026-02-18 | Run ID `20260218-183025`; unified metrics + side-by-side CSV generated |
| EX-005 | Add ParmEst adapter for Pyomo 6.9.5 solver-token compatibility | Completed | 2026-02-18 | Commit `bfe7992`; `unified_codebase_library.py` |
| EX-006 | Populate paper baseline values for Stage A/Stage C | Completed | 2026-02-18 | Updated `UnifiedFramework/DATA3/docs/validation/paper_reference_values_template.csv` (DATA1 from legacy fit file; DATA2 from published-era notebook/script artifacts) |
| EX-007 | Compute tolerance pass/fail status from populated baselines | Completed (Rerun with source-backed DATA2) | 2026-02-18 | Run `20260218-190331`; comparison table at `UnifiedFramework/DATA3/results/reproduction/20260218-190331/tables/paper_vs_unified_side_by_side.csv` |
| EX-008 | Expand to Stage B full DATA1 targets | Completed (covered by EX-016) | 2026-02-20 | Stage B generation implemented and executed in run `20260220-223843` |
| EX-009 | Improve DATA1/DATA2 fit parity via dataset-specific model forms and parameter mapping | In Progress | 2026-02-18 | Updated DATA2 runs to `Lag + convection`; switched comparison to `beta_0/beta_1`; exploratory constrained run `20260218-191502` |
| EX-010 | Legacy execution-parity audit against upstream `utility.py` + notebook flow | Completed | 2026-02-19 | Mapped objective/suffix rules and patched unified gating (`n_v0`, `n_extra`, zero-skip, scalar/series `cV_avg`) in `unified_codebase_library.py`; smoke run `20260219-004414` |
| EX-011 | Implement interpolation-consistent measurement mapping for ParmEst suffix workflow | Completed | 2026-02-19 | Added interpolation-linked measurement proxy variables/constraints and switched suffix labeling to proxy outputs; runs `20260219-121031`, `20260219-121117` |
| EX-012 | Align DATA1 fit model form/parameterization with legacy utility.py | Completed (major parity gains) | 2026-02-19 | Added DATA-mode legacy boundary/linking behavior, legacy grouped objective fit path, MAT `ni` mapping, and DATA1 non-logit sigma setting; run `20260219-164626` |
| EX-013 | Use MAT-native initialization by default + ipopt-only ParmEst path + XLSX multistart hook | Completed | 2026-02-19 | Added `build_guess_from_experiment_v24`, removed `ef_ipopt` bridging, and integrated optional `pyomo.contrib.multistart` for XLSX single-experiment fallback; run `20260219-171416` |
| EX-014 | Replace external multistart dependency with internal restart loop and keep common MAT/XLSX estimation flow | Completed | 2026-02-20 | Swapped fallback from `pyomo.contrib.multistart` solve call to internal randomized restart loop; preserved MAT-native initialization and shared estimation pipeline; run `20260220-221507` |
| EX-015 | Enable utility-style restart sweep for DATA1 MAT legacy-fit path | Completed | 2026-02-20 | Added optional MAT-legacy restart flag and ran DATA1 sweep (`multistart_iterations=30`); run `20260220-223530` |
| EX-016 | Execute Stage B DATA1 artifact generation (main Fig.4/5/6 + SI S2-S7) | Completed | 2026-02-20 | Added Stage B generator + manifest in scaffold and produced artifacts in run `20260220-223843` |
| EX-017 | Deterministically regenerate DATA1 notebook outputs and produce panel checklist | Completed | 2026-02-21 | Added `UnifiedFramework/DATA3/scripts/reproduce/regenerate_data1_notebook_artifacts.py`; run `20260221-data1-notebook` generated 44 PNGs with manifest + panel checklist |
| EX-018 | Lock DATA1 main/SI figure mappings by PDF visual comparison and update matrix | Completed | 2026-02-21 | Locked SI/main mappings in `data1_panel_checklist.md`; updated `paper_target_matrix.md` DATA1 figure statuses to `Locked (20260221-data1-notebook)` |
| EX-019 | Rename unified modules + add shared MAT/XLSX pipeline entrypoint | Completed | 2026-02-21 | Renamed to `unified_codebase_library.py` and `unified_codebase_runfile.py`, added `UnifiedPipelineConfigV24` + `run_unified_pipeline_v24`, and preserved legacy wrappers |
| EX-020 | Deepen unified MAT/XLSX orchestration internals behind one API contract | Completed | 2026-02-21 | Added source-aware option resolver, stage helpers (`_run_unified_estimation_stage_v24`, `_run_unified_doe_stage_v24`), run metadata builder, and optional reporting hook in `run_unified_pipeline_v24` |
| EX-021 | Implement profile-based runner switching (DATA1/DATA2/DATA3) and document usage | Completed | 2026-03-05 | Updated `unified_codebase_runfile.py` with CLI profile switching and refreshed `UnifiedFramework/DATA3/docs/unified_code/unified_codebase_runfile.md` |
| EX-022 | Separate `NOT_APPLICABLE` from `MISSING_VALUE` for consolidated numeric validation status | Completed | 2026-03-05 | Added `UnifiedFramework/DATA3/scripts/validation/update_target_validation_status_consolidated.py` and refreshed `target_validation_status_consolidated.csv` from tolerance rows |
| EX-023 | Enable full paper extraction runtime and build target-to-page evidence index | Completed | 2026-03-06 | Installed PDF/OCR tooling; extraction run `20260306-023356-paper-pdf-extract`; `UnifiedFramework/DATA3/docs/validation/target_pdf_page_index.csv` (`32/32` targets matched) |
| EX-024 | Lock DATA2 main/SI figure panel mappings (prioritize main paper first) | Completed (figures) | 2026-03-06 | Locked `D2-M-F8`, `D2-M-F9`, `D2-S-FS2..FS7` in source map and consolidated status using generated artifacts + paper page evidence |
| EX-025 | Generate dedicated main-paper table artifacts for DATA2 Tables 3-6 with warning semantics | Completed (warning baseline) | 2026-03-06 | Added `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/data2_table3_side_by_side.csv` ... `data2_table6_side_by_side.csv`; set `numeric_validation_status=MISSING_VALUE` with WARNING notes |
| EX-026 | Reclassify `D2-M-VARS` as reference-only integrity target | Completed | 2026-03-06 | Updated matrix/source-map/consolidated/checklists to `NOT_APPLICABLE (reference-only: experimental inputs/measurements)` and removed PASS/FAIL gating expectation |

## Immediate Next Actions
1. Close numeric reproduction gaps for main-paper targets first:
   - DATA1: `D1-M-T1` (`FAIL`)
   - DATA2: `D2-M-T2` (`FAIL`), `D2-M-F6` (`FAIL`), `D2-M-T3..T6` (`MISSING_VALUE`)
2. Resolve `run_DATA2_model_variations.py` runtime blocker (Pyomo overflow under Python 3.12) to populate unified-side values for Tables 3-6.
3. Refresh tolerance evaluation and consolidated status after table-value generation.
4. Then address SI table `D2-S-TS1` (if still unresolved) and finalize one consolidated validation report keyed by target ID.
5. Add nightly guardrails that fail on `MISSING_VALUE` only for numeric-applicable targets (ignore `NOT_APPLICABLE` rows, including `D2-M-VARS`).

## Open Items to Track
- Covariance behavior when `k_aug` is unavailable in environment (documented limitation path).
- Any solver option tuning needed for stable covariance/fit across datasets.
- Final plotting parity details (fonts/line styles/layout per paper).
- Repo cleanup rename/move proposal document (pending).
- Investigate DATA1 Stage A objective and parameter mismatches (objective, `Lp`, `B`) vs legacy baseline.
- Align DATA2 comparison metrics with paper model form (`beta_0`, `beta_1`, optional `S0`) since scalar `theta.B` is not directly comparable.
- Add dataset-specific initial guesses/fixed-parameter controls in scaffold to avoid boundary convergence for DATA2 (`beta_0`, `beta_1`, `sigma`).
- Runtime blocker: `run_DATA2_model_variations.py` currently fails in this environment with Pyomo NL writer overflow (`OverflowError: Python integer -1 out of bounds for uint8`), preventing unified-side numeric fill for DATA2 Tables 3-6.

## Progress Log
- 2026-02-18: Stage A executed with `ipopt` on `DATA1_511.12` via `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`.
  - Run ID: `20260218-182946`
  - Generated: Fig.2 set, Table 1 side-by-side CSV, Fig.3 overlay, Table 2 side-by-side CSV.
  - Note: side-by-side tolerance status remains `MISSING_VALUE` until paper baseline values are populated in `paper_reference_values_template.csv`.
- 2026-02-18: Stage C core executed with `ipopt` on `DATA2_270511.123` and `DATA2_270611.123`.
  - Run ID: `20260218-183025`
  - Generated: unified metrics and side-by-side comparison CSV (baseline paper values still pending).
- 2026-02-18: EX-006 baseline population updated and EX-007 first pass executed.
  - Run ID: `20260218-185741`
  - Status summary: `PASS=7`, `FAIL=5`, `MISSING_VALUE=3`.
  - Failures observed in DATA1 Stage A (`objective`, `theta.Lp`, `theta.B`) and objective rows for DATA2 (with provisional baseline values).
- 2026-02-18: EX-007 rerun after replacing DATA2 provisional baselines with source-backed values.
  - Run ID: `20260218-190331`
  - Status summary: `PASS=1`, `FAIL=8`, `MISSING_VALUE=6`.
  - DATA2 `theta.B` and some objective rows remain `MISSING_VALUE` where no scalar/table-equivalent baseline exists in current artifacts.
- 2026-02-18: EX-009 exploratory improvement run with dataset-specific DATA2 form.
  - Run ID: `20260218-191502`
  - Changes tested: DATA2 `mode=Lag`, `b_form=convection`, tightened `beta` bounds, fixed `sigma` in estimation.
  - Outcome: objective dropped substantially, but parameter parity remains outside tolerance; indicates local-minimum/weighting mismatch rather than simple solver configuration issue.
- 2026-02-18: EX-009 redesign run with cross-validation metric switch.
  - Run ID: `20260218-192333`
  - Changes implemented in scaffold:
    - dataset-specific execution mode (`estimation` for `DATA2_270511.123`; fixed-parameter cross-validation for `DATA2_270611.123`)
    - DATA2 metric mapping changed to curve-quality targets (`nrmse.mass`, `nrmse.cF`, `nrmse.cV`)
    - side-by-side merge changed to left-join on references only (prevents non-target metric inflation).
  - Outcome: comparison semantics improved, but tolerance failures remain significant (`FAIL=6`, `PASS=1`, `MISSING_VALUE=3`).
- 2026-02-18: EX-009 stability fixes and re-run.
  - Run ID: `20260218-194541` (DATA1 only quick check)
  - Fixes:
    - removed unsupported IPOPT option (`acceptable_iter`) in curve-simulation retry path
    - removed duplicate `objective_legacy` insertion in metric table assembly
  - Outcome: pipeline stability improved; numeric mismatch for DATA1 remains and requires objective/measurement parity deep-dive.
- 2026-02-19: EX-010 upstream legacy parity audit + unified objective/suffix gating patch.
  - Sources reviewed: upstream `README.md`, `utility.py`, `DATA1_model_demo.ipynb`, `DATA2_model_demo.ipynb`, `DATA2_visualization.ipynb`.
  - Unified updates in `unified_codebase_library.py`:
    - persisted MAT gating metadata (`n_v0`, `n_extra`, `n`) on `ExperimentalData`
    - applied legacy vial gates to mass/cF/cV measurement objective terms
    - applied matching gates to ParmEst suffix labeling
    - added zero-value skip logic and scalar-or-series handling for `cV_avg`
  - Validation run: `20260219-004414` (DATA1 only smoke run) completed end-to-end.
  - Outcome: execution parity improved structurally; numerical parity still requires additional alignment (likely interpolation-vs-discrete labeling and DATA1 legacy model-form details).
- 2026-02-19: EX-011 interpolation-consistent measurement mapping implemented and tested.
  - Unified updates in `unified_codebase_library.py`:
    - added linear interpolation helper over discretized `tau` nodes for arbitrary measurement times
    - added measurement-record builder with legacy gating/weighting rules
    - added measured-output proxy variables and linking constraints (`_obs_pred == interpolated state`)
    - switched ParmEst/DoE suffix labeling to use proxy outputs (instead of nearest-node state values)
  - Validation runs: `20260219-121031`, `20260219-121117` (DATA1 Stage A only).
  - Outcome: interpolation parity mechanism is active and stable, but DATA1 Stage A fit parity did not improve yet (`theta.B` still at lower bound; objective and `theta.Lp` remain outside tolerance).
- 2026-02-19: EX-012 DATA1 legacy fit alignment implemented.
  - Unified updates in `unified_codebase_library.py`:
    - DATA-mode MAT initial `cH[1,0]` now follows legacy (`0.8 * first-vial cV_avg`)
    - DATA-mode legacy vial-linking constraints for `mV` and `cVmV` restored
    - MAT loader now maps `data_config.ni` (ionic species count) to model osmotic term
    - added direct legacy grouped objective fit path for DATA/MAT/single-B runs
    - added utility-style casadi initialization attempt before discretized legacy fit
  - Script update in `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`:
    - DATA1 canonical run now uses non-logit sigma parameterization
  - Validation run: `20260219-164626` (DATA1 Stage A, `nfe=300`).
  - Outcome vs baseline:
    - `theta.B`: from bound (`1e-6`) to `0.2842` (relative error ~7.54%)
    - `theta.Lp`: improved to `3.2900` (relative error ~10.58%)
    - `objective_legacy`: improved to `9.59e5` (relative error ~30.2%)
- 2026-02-19: EX-013 initialization/solver framework modernization.
  - Unified updates in `unified_codebase_library.py`:
    - added `build_guess_from_experiment_v24` to source default guesses from MAT `theta0/Lp0/B0/sigma0`
    - estimation path now avoids `ef_ipopt` mapping and remains `ipopt`-only in unified code
    - added optional XLSX-only multistart hook (`ModelOptions.use_multistart_for_xlsx`, `multistart_iterations`, `multistart_strategy`) using `pyomo.contrib.multistart`
  - Runner update in `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`:
    - DATA1 canonical case now uses MAT-derived initialization by default when no explicit override is provided
  - Validation run: `20260219-171416` (DATA1 Stage A quick run) completed successfully with MAT-derived initialization.
- 2026-02-20: EX-014 internal restart fallback + common flow cleanup.
  - `unified_codebase_library.py`:
    - replaced runtime `multistart` solver call with internal utility-style randomized restart loop in ipopt fallback
    - retained `ModelOptions` restart controls so this can be refactored to official multistart later
    - kept MAT initial-value usage (`theta0/Lp0/B0/sigma0`) through shared `build_guess_from_experiment_v24`
  - outcome: no regression in DATA1 quick validation (`20260220-221507`), and fallback path no longer depends on external multistart availability.
- 2026-02-20: EX-015 MAT restart sweep enabled.
  - `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py`:
    - DATA1 canonical spec now enables MAT-legacy restart sweep (`use_multistart_for_mat_legacy=True`, `multistart_iterations=30`)
  - `unified_codebase_library.py`:
    - legacy DATA/MAT fit path now optionally uses the shared internal restart routine
  - validation run: `20260220-223530` (DATA1 Stage A quick run).
  - outcome: converged to the same best basin as prior MAT-initialized fit (`Lp~3.294`, `B~0.283`, objective~9.56e5), indicating remaining mismatch is structural/data-treatment, not local-minimum search failure.
- 2026-02-20: EX-016 Stage B DATA1 artifact generation completed.
  - `UnifiedFramework/DATA3/scripts/reproduce/reproduce_data1_data2.py` updates:
    - added `_build_data1_stage_b_outputs` to generate:
      - DATA1 main: Fig.4 (`data1_fig4_*`), Fig.5 (`data1_fig5_*`), Fig.6 (`data1_fig6_*`)
      - DATA1 SI: Fig.S2-S7 (`data1_si_figS*_*.png`)
    - added `data1_stageB_artifact_manifest.csv` with target-to-file mapping
    - wired Stage B generation into `main()` metadata output
  - run: `20260220-223843`.
  - `paper_target_matrix.md` updated to mark DATA1 Stage B targets as `Generated (20260220-223843)`.
- 2026-02-21: EX-017 deterministic DATA1 notebook regeneration and panel checklist completed.
  - added script: `UnifiedFramework/DATA3/scripts/reproduce/regenerate_data1_notebook_artifacts.py`
  - run: `20260221-data1-notebook`
  - output manifest: `UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/figures/data1_notebook_regen/manifest_data1_notebook_regen.json`
  - output checklist: `UnifiedFramework/DATA3/docs/validation/data1_panel_checklist.md`
  - generated artifacts: 44 PNGs (main + SI-oriented notebook outputs)
- 2026-02-21: EX-018 DATA1 figure mapping lock completed.
  - locked main/SI mapping via PDF image extracts under:
    - `UnifiedFramework/DATA3/results/reproduction/20260221-data1-notebook/pdf_extract/data1_si/`
  - updated figure mapping statuses:
    - `UnifiedFramework/DATA3/docs/validation/data1_panel_checklist.md`
    - `UnifiedFramework/DATA3/docs/validation/paper_target_matrix.md`
- 2026-02-21: EX-019 unified module rename + shared pipeline API baseline completed.
  - module rename:
    - `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_library.py`
    - `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/unified_codebase_runfile.py`
  - compatibility wrappers retained:
    - `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/archives/experiment_dataload_OOP_v24.py`
    - `UnifiedFramework/DATA3/ExperimentalDataAnalysis/UnifiedCode/experiment_load_v24.py`
  - shared MAT/XLSX pipeline API added:
    - `UnifiedPipelineConfigV24`
    - `run_unified_pipeline_v24(...)`
- 2026-02-21: EX-020 unified orchestration internals refactor completed.
  - `unified_codebase_library.py`:
    - added source-aware model options defaults and resolver for MAT vs XLSX
    - added shared stage helpers for estimation and DoE execution
    - added run metadata builder including stage/target IDs/artifact-root context
    - extended `UnifiedPipelineConfigV24` with reporting fields and callback hook
    - updated `run_unified_pipeline_v24` to route through helper stages with one API contract
- 2026-03-06: EX-023 to EX-026 validation consolidation completed (main-paper priority).
  - extraction + evidence:
    - run `20260306-023356-paper-pdf-extract` produced page/image/text artifacts for DATA1/DATA2 main+SI papers
    - `UnifiedFramework/DATA3/docs/validation/target_pdf_page_index.csv` now maps all targets (`32/32`) to matched paper pages
  - figure/panel lock progress:
    - main paper locked: `D2-M-F8`, `D2-M-F9`
    - SI locked: `D2-S-FS2..FS7` via cross-verification outputs + panel visual checks
  - main-paper tables:
    - generated warning-baseline artifacts for `D2-M-T3..T6` under `UnifiedFramework/DATA3/results/reproduction/20260306-mainpaper-table-baseline/tables/`
    - statuses kept as `MISSING_VALUE` until unified-side values can be produced in a stable runtime
  - semantics correction:
    - `D2-M-VARS` set to `NOT_APPLICABLE` as reference-only metadata integrity target (not reproducible numeric result)
  - current replication snapshot from consolidated status:
    - DATA1: `FAIL=1` (`D1-M-T1`), `NOT_APPLICABLE=12`
    - DATA2: `FAIL=2` (`D2-M-T2`, `D2-M-F6`), `MISSING_VALUE=4` (`D2-M-T3..T6`), remainder `NOT_APPLICABLE`
