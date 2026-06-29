# Prompt: Paper-Comparison Pytest Layer for DATA1 / DATA2

> Authored 2026-06-09 during the diffusivity-coefficient model-update discussion.
> This is the brief that drove the implementation in
> `tests/test_data1_paper_comparison.py` + `tests/test_data2_paper_comparison.py`.
> Kept as a permanent artifact so the scope, decisions, and trade-offs are
> auditable.

## Origin

While planning a model-physics change (adding salt-specific diffusion coefficients D
to the Sherwood mass-transfer correlation, Eq. 4 of Ouimet et al. 2022), the
operator needs a deterministic way to answer:

> "After I change the process model, do the DATA1 and DATA2 figures I would
> reproduce still match the published-paper figures, or did something move?"

Currently the pytest suite catches **structural** drift (Pyomo Var/Constraint
counts, loader scalars). It does NOT catch the case where the model preserves
its structure but the numeric outputs shift — which is exactly the kind of
change a diffusivity coefficient edit produces.

## What to build

A new pytest layer that compares the outputs of the refactored DATA1 / DATA2
workflows to the published paper figures, using the locked figure mapping and
the committed numeric baselines that already exist in the repo.

## Inputs (read-only, all pre-existing)

| Asset | Path | Purpose |
|---|---|---|
| Locked figure mapping | `UnifiedFramework/DATA3/docs/validation/target_notebook_source_map.csv` | 33 rows: `target_id` → paper figure label → exact filename in `paper_artifacts/data{1,2}/notebook_figures/` |
| Numeric baselines | `UnifiedFramework/DATA3/docs/validation/simulation_validation_data_files/` | Per-figure CSV baselines (23 rows DATA1, 33 rows DATA2) generated from the legacy MATLAB stack |
| Numeric manifest | `simulation_validation_manifest.csv` | Per-row paper_id, figure_id, panel_id, baseline_csv path, artifact_kind |
| Paper figure renders | `UnifiedFramework/DATA3/results/reproduction/20260306-023356-paper-pdf-extract/pdf_extract/` | Rasterized pages (`data1_main/pages/page-NN.png`) + figure crops (`data2_{main,si}/images/img-NNN.png`). NOTE: `data1_main` has no crops, only full pages. |
| Source PDFs | `UnifiedFramework/DATA3/published_works/{DATA1_main,DATA1_SI,DATA2_main,DATA2_SI}.pdf` | Human-readable reference; NOT consumed by the pytest. |

## Outputs

1. **`tests/test_data1_paper_comparison.py`** — parametrized over DATA1 rows of the manifests.
2. **`tests/test_data2_paper_comparison.py`** — parametrized over DATA2 rows.
3. **`helpers/paper_comparison.py`** — shared loading + comparison utilities.
4. **`_runs/run-*/paper_comparison_report.csv`** — per-target drift table emitted each run.

## Scope decision: numeric-only first

Three options were on the table:
- (a) Numeric only — CSV-vs-CSV comparison via `simulation_validation_manifest.csv`. ~30 min build.
- (b) Numeric + visual — add SSIM image similarity against `pdf_extract/`. ~1 hr build, heavier dep on Pillow / scikit-image.
- (c) Mapping audit first — no test code, just a report of which targets are covered.

**Decision: ship (a) now; document (b) as a future extension.**

Rationale: the diffusivity edit will move *numbers*, not the *style* of a plot.
A pure-CSV check answers the operator's question with the lowest dependency
footprint and zero PDF reading. SSIM image checks can be added later as a
separate test file without touching this one.

## Three sub-layers within the comparison layer

The comparison work decomposes cleanly:

| Sub-layer | Marker | What it does |
|---|---|---|
| Smoke | `smoke` | Asserts every baseline CSV named in `simulation_validation_manifest.csv` exists on disk and has the expected columns. No produced artifacts read. < 1 s. |
| Light regression | `regression` | If the operator has *already* run `materialize_all(campaign="DATA1"/"DATA2")` (e.g. they ran the deck-build pipeline yesterday), compare the produced CSV/PNG against the baseline WITHOUT re-running the solver. Auto-skipped if produced artifacts are absent. |
| Full regression | `regression + nightly` | Cached fixture runs `materialize_all` once per family, then runs the same comparison. Heavy: 10 min / family. |

The light regression sub-layer is what the operator actually wants — it doesn't
need a fresh solver run, just compares whatever's currently in
`paper_artifacts/data{1,2}/notebook_figures/`. After the diffusivity edit:

```bash
# 1. operator re-runs materialize_all manually (or skip if recent run exists)
python refactored_codes_v1/refactored_ucb_runfile.py  # pick DATA1, then DATA2

# 2. pytest compares to baselines — fast, no solver
MPLCONFIGDIR=/tmp/mplconfig pytest -q -m "regression and not nightly" \
    pytest_refactored_codes_v1/tests/test_data1_paper_comparison.py \
    pytest_refactored_codes_v1/tests/test_data2_paper_comparison.py
```

## Metrics + tolerance

Per-row from `simulation_validation_manifest.csv`:
- `artifact_kind == "timeseries"` → NRMSE between time-aligned series, tolerance 1% default
- `artifact_kind == "scatter"` → MAE on points, tolerance 1% default
- `artifact_kind == "contour_grid"` → max absolute relative diff over the grid, tolerance 2% default
- `artifact_kind == "table"` → cell-by-cell relative diff, tolerance 1% default

Tolerances live in `baselines/paper_comparison_tolerances.csv` (override per target).
Rows without an explicit tolerance default to 1%.

## Drift report format

After each test run, `_runs/run-<id>/paper_comparison_report.csv`:

```
paper,figure_id,panel_id,artifact_kind,target_id,baseline_csv,produced_csv,metric,value,tolerance,verdict,notes
DATA1,Fig. 2,A,timeseries,D1-M-F2_A,data1_main/fig2/data1_main_fig2_a.csv,paper_artifacts/data1/notebook_figures/mass-dat501.1.png.csv,nrmse,0.003,0.01,PASS,
DATA1,Fig. 4,A-mass,timeseries,D1-M-F4_A_mass,data1_main/fig4/data1_main_fig4_a_mass.csv,…,nrmse,0.014,0.01,FAIL,diffusivity edit candidate
```

A `FAIL` here is the answer to "did my edit change DATA1?"

## Non-goals

- **No edits to legacy code paths.** `simulation_validation_data_files/` is consumed read-only.
- **No re-extraction of paper figures.** `pdf_extract/` is reused as-is for (b) if/when it's built.
- **No PDF parsing.** The pytest never reads the source PDFs.
- **No visual layer in v1.** SSIM / image similarity is deferred to a sibling test file.
- **No edits to `target_notebook_source_map.csv`.** It's the locked authority on mapping.
- **No new baselines generated.** Refreshing baselines is a separate manual workflow run from the legacy stack.

## Success criteria

1. `pytest -m smoke pytest_refactored_codes_v1/tests/test_data{1,2}_paper_comparison.py` passes
   in < 2 s, asserting every baseline CSV in the manifest is on disk with the right columns.
2. `pytest -m "regression and not nightly" pytest_refactored_codes_v1/tests/test_data{1,2}_paper_comparison.py`
   runs in < 30 s when `paper_artifacts/data{1,2}/notebook_figures/` is populated.
3. After a deliberate model edit, the drift report names every figure that moved
   AND quantifies the drift.
4. The layer integrates cleanly with the existing 4-layer architecture documented
   in the suite's README.md — same markers, same `_runs/` convention, same
   conftest fixtures.
5. The pytest never tries to run the solver itself. The cached fixture for the
   full-regression sub-layer is opt-in via the `nightly` marker.

## What this catches (and what it doesn't)

CATCHES:
- Diffusivity coefficient enters Eq. 4 with the wrong exponent → mass-transfer
  coefficient shifts → concentration polarization shifts → DATA1 Fig. 2 panel
  drifts. The numeric comparison flags it.
- Salt-specific D added but KCl value not preserved → DATA1 + DATA2 numbers move
  uniformly.
- A loader fix that changes the data_stru scalar (e.g. C_F0_value) → numeric
  fit moves.

MISSES:
- A pure visual change (e.g. axis labels reformatted, color palette change).
  → That's what the visual sub-layer (b) would add. Out of scope for v1.
- Drift that's smaller than the per-target tolerance. Tune tolerances in
  `baselines/paper_comparison_tolerances.csv`.
- New panels added to a figure that weren't in the original paper. The locked
  mapping only knows about published panels.
