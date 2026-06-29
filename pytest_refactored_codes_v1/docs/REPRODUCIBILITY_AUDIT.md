# Reproducibility Audit — DATA1 + DATA2 (Refactored Codebase)

> Status: COMPLETED on 2026-06-09.
> Companion artifact: [`REPRODUCIBILITY_AUDIT_PROMPT.md`](./REPRODUCIBILITY_AUDIT_PROMPT.md) (the scope spec).
> Scope locked: every Figure and Table in the four PDFs at
> `UnifiedFramework/DATA3/published_works/{DATA1_main, DATA1_SI, DATA2_main, DATA2_SI}.pdf`.
> Authority codebase: `refactored_codes_v1/{refactored_ucb_library.py, refactored_ucb_runfile.py, conductivity_paper.py}`.
> Read-only audit. No changes to any code.

---

## Executive summary

Of **41 published items** (24 figures, 8 tables, plus 9 SI sub-figures) across the four PDFs:

| Verdict | Count | What it means |
|---|---|---|
| `PRESENT` | 0 | Reproducible from scratch using ONLY `refactored_codes_v1` (no `.mat` reads, no V24 imports). |
| `PRESENT_WITH_MAT_DEPENDENCY` | 21 | Python wrapper exists but reads pre-computed MATLAB outputs (`.mat`, `contourdata-*.csv`, `fit_stru.mat`). |
| `PRESENT_WITH_V24_DEPENDENCY` | 1 | Python "direct" path exists but imports from `UnifiedFramework/DATA3/.../UnifiedCode/unified_codebase_library` — leaks outside `refactored_codes_v1`. |
| `MISSING_NEEDS_MATLAB_PORT` | 10 | No Python equivalent; needs a named MATLAB script ported into `refactored_ucb_library.py`. |
| `MISSING_OUT_OF_SCOPE` | 9 | Not model-derived (apparatus schematics, workflow diagrams, paper-image alias from `pdf_extract/`). |

**The bottom line:** **zero** items in the published papers can currently be reproduced "from scratch" using only `refactored_codes_v1`. Every model-derived figure depends on either pre-computed MATLAB artifacts on disk or on the V24 UnifiedCode stack outside `refactored_codes_v1`. This makes the suite UNABLE to detect drift if MATLAB-generated artifacts are deleted or out-of-sync — exactly the fragility the operator is concerned about.

**Top-3 MATLAB ports that unlock the most figures** (each row, when ported, removes one or more `PRESENT_WITH_MAT_DEPENDENCY` flags):

| Port | MATLAB source(s) | Figures it unlocks |
|---|---|---|
| 1. Sigma sensitivity sweep | `legacy/data1_matlab/functions/sigma_sensitivity.m` | DATA1 Fig. 4 (and the `.mat` reads inside `run_sigma_sensitivity`) |
| 2. 2-D contour grid (Lp×σ, Lp×B) | `legacy/data1_matlab/functions/calc_contour_2d.m` + `plot_contour.m` | DATA1 Figs. 5, 6, S3, S4, S5, S6 |
| 3. DoE heatmaps (A/D/E/modified-E optimality) | `legacy/data1_matlab/functions/` (no direct `.m` exists in repo — need `doe_heatmap_filtration.m` + `doe_heatmap_diafiltration.m` reconstructed from the MATLAB README) | DATA1 Fig. S7 |

---

## How to read the table

- **Verdict** uses the four-element vocabulary from the audit prompt.
- **Python function (existing)** is the symbol in `refactored_ucb_library.py` that produces this item today.
- **`.mat` / `.csv` dependency** names the pre-computed file the wrapper reads.
- **MATLAB source** is the script in `legacy/data1_matlab/functions/` that originally produced the dependency.
- **Python reference** is the closer pre-refactor Python source (notebook cell, `utility.py`) — context, not authority.

Source-of-truth crosswalks for the column data:
- Library function index: `grep -nE "^def (run_|plot_|render_)" refactored_codes_v1/refactored_ucb_library.py`
- Manifest entries: `DATA1_FIGURES` / `DATA2_FIGURES` at lines 11657 / 11736 of `refactored_ucb_library.py`
- MATLAB sources: `ls legacy/data1_matlab/functions/`

---

## DATA1_main (J. Membr. Sci. 641, 119743, 2022)

| # | Item | Description | Verdict | Python function | Dependency | MATLAB source | Python reference (legacy) |
|---|---|---|---|---|---|---|---|
| 1 | Fig. 1 | Apparatus schematic | `MISSING_OUT_OF_SCOPE` | — | — | — | — |
| 2 | Fig. 2 (A–D) | Mass + cV time series for filtration 501.1 + diafiltration 511.12, Membrane A | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_figure2_workflow` + `plot_sim_comparison` | `fit_stru.mat` per dataset under `legacy/data1_matlab/functions/dat 501.1.../fit_stru.mat`, `dat 511.12.../fit_stru.mat` | `run_data_analysis.m` (produces the `fit_stru.mat` from raw .mat) | `legacy/data1_matlab/DiafiltrationPaperPlots.ipynb` |
| 3 | **Table 1** | Filtration parameters: 3 membranes × 4 models × (Lp, B, σ, 3 objectives) | `MISSING_NEEDS_MATLAB_PORT` | none | — | `run_data_analysis.m` (writes the WLS fit table) | `utility.py` (fit + report) |
| 4 | Fig. 3 | Permeate vs retentate concentration scatter (filtration △ + diafiltration ▢, 5 membranes) | `PRESENT_WITH_MAT_DEPENDENCY` | `plot_conc_range` via `run_data1_concentration_comparison` | reads `experiment space/Classical_analysis-dat*.csv` + `diafiltration.csv` (produced by MATLAB) | `run_data_analysis.m` writes the `experiment space/` CSVs | `legacy/data1_matlab/diafiltration_plots.py` |
| 5 | **Table 2** | Diafiltration parameters: 3 membranes × 4 models × (Lp, B, σ, 3 objectives) | `MISSING_NEEDS_MATLAB_PORT` | none | — | `run_data_analysis.m` | `utility.py` |
| 6 | Fig. 4 (A, B) | σ-sensitivity (σ = 0.1, 0.5, 0.9): 2 rows × 3 cols (mass / permeate / retentate, filt. vs diafilt.) | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_figure4_workflow` → `run_sigma_sensitivity` | `sigma sensitivity/sim_stru-dat<dataset> C_Fin<cf0>sig<sigma>.mat` (6 files) | **`sigma_sensitivity.m`** + `heatmap_sigma_sensitivity_filtration.m` + `heatmap_sigma_sensitivity_diafiltration.m` | `legacy/data1_matlab/diafiltration_plots.py` |
| 7 | Fig. 5 (A, B, C × 3) | Lp×σ residual contours for DIAFILTRATION: 3 models × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_figure5_workflow` | `contourdata-x_sigma-y_Lp.csv` under `511.12 concpolar/`, `511.11 concpolar/`, `511.12/` | **`calc_contour_2d.m`** + `plot_contour.m` + `plot_grid_contour.m` | `legacy/data1_matlab/DiafiltrationPaperPlots.ipynb` |
| 8 | Fig. 6 (A, B, C × 3) | Lp×B residual contours for DIAFILTRATION: 3 models × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_figure6_workflow` | `contourdata-x_B-y_Lp.csv` under same dirs | `calc_contour_2d.m` + `plot_contour.m` | as above |
| 8b | (alt path) | direct contour computation for Figs. 5/6 | `PRESENT_WITH_V24_DEPENDENCY` | `run_data1_direct_contour_branch` (line 7476) | **imports from `UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library`** (`evaluate_data1_paper_contour_objectives_v24`, `simulate_data1_vialwise_trajectories_v24`) | — | — |

---

## DATA1_SI

| # | Item | Description | Verdict | Python function | Dependency | MATLAB source | Python reference |
|---|---|---|---|---|---|---|---|
| 9 | S1 / Eqs S1-S9 | Model derivation (textual + equations) | `MISSING_OUT_OF_SCOPE` | — | — | — | — |
| 10 | Fig. S1 | KCl calibration curve (conductivity → concentration) — single panel, linear regression | `MISSING_NEEDS_MATLAB_PORT` | none in refactored library (DATA2 has `run_data2_calibration_plots` but it's NF270/KCl two-board; not the DATA1 single curve) | calibration CSV not in repo | — | `utility.py` calibration helpers |
| 11 | Fig. S2 (A, B, C) | Inline-probe vs mass-balance retentate concentration — 3 panels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_si_s2` (line 6653) | reads `data_stru-dataset{301.1,501.1,511.12}.mat` + `fit_stru.mat` + `experiment space/Classical_analysis-dat<id>.csv` | `run_data_analysis.m` + `load_data.m` | `legacy/data1_matlab/DiafiltrationPaperPlots.ipynb` |
| 12 | **Table S1** | Filtration vs diafiltration time comparison (8 timing rows × 2 modes) | `MISSING_NEEDS_MATLAB_PORT` | none | — | manual table in paper (no `.m` exists) | — |
| 13 | Fig. S3 (A,B,C,D × 3) | Lp×σ residual contours for DIAFILTRATION: 4 models (M1-M4) × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_si_s3` (line 6871) | `contourdata-x_sigma-y_Lp.csv` under `511.12 concpolar/`, `511.11 concpolar/`, `511.12/`, `511.11/` | `calc_contour_2d.m` + `plot_contour.m` | notebook |
| 14 | Fig. S4 (A,B,C,D × 3) | Lp×B residual contours for DIAFILTRATION: 4 models × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_si_s4` (line 6882) | `contourdata-x_B-y_Lp.csv` under same dirs | `calc_contour_2d.m` | notebook |
| 15 | Fig. S5 (A,B,C,D × 3) | Lp×σ residual contours for FILTRATION: 4 models × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_si_s5` (line 6893) | `contourdata-x_sigma-y_Lp.csv` under `501.1 concpolar/`, `501.11 concpolar/`, `501.1/`, `501.11/` | `calc_contour_2d.m` | notebook |
| 16 | Fig. S6 (A,B,C,D × 3) | Lp×B residual contours for FILTRATION: 4 models × 3 channels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data1_si_s6` (line 6904) | `contourdata-x_B-y_Lp.csv` under same dirs | `calc_contour_2d.m` | notebook |
| 17 | Fig. S7 (A,B,C,D × 3) | DoE heatmaps: ΔP vs cf (or cd) for 4 experiment scopes × 3 channels (σ-discrimination) | **`MISSING_NEEDS_MATLAB_PORT`** | none — no Python equivalent in refactored library | — | **`doe_heatmap_filtration.m`** (Panel A) + **`doe_heatmap_diafiltration.m`** (Panels B-D) — NOT present in `legacy/data1_matlab/functions/` (the README describes them but the `.m` files are not in this checkout). MATLAB README at `legacy/data1_matlab/README` describes their behavior. | `legacy/data1_matlab/PSE2021Plots.ipynb` (the MBDoE paper) |

---

## DATA2_main (Ind. Eng. Chem. Res. 64, 12111-12130, 2025)

| # | Item | Description | Verdict | Python function | Dependency | MATLAB source | Python reference |
|---|---|---|---|---|---|---|---|
| 18 | Fig. 1 | NF270 apparatus schematic | `MISSING_OUT_OF_SCOPE` | (resolves to `pdf_extract/data2_main/images/img-NNN.png` via `run_data2_publication_figures`) | — | — | — |
| 19 | Fig. 2 (A, B) | Lag + overflow startup mode diagrams | `MISSING_OUT_OF_SCOPE` | (paper-image alias from `pdf_extract/`) | — | — | — |
| 20 | Fig. 3 | Time corrections diagram (apparatus + mass extrapolation inset) | `MISSING_OUT_OF_SCOPE` | (paper-image alias) | — | — | — |
| 21 | Fig. 4 | MBDoE block diagram | `MISSING_OUT_OF_SCOPE` | — | — | — | — |
| 22 | Fig. 5 | FIM geometric interpretation (confidence ellipsoid) | `MISSING_OUT_OF_SCOPE` | — | — | — | — |
| 23 | **Table 1** | Manipulated and measured variables | `MISSING_OUT_OF_SCOPE` | — (static table content) | — | — | — |
| 24 | **Table 2** | Parameter results for B = Jw·(β0 + β1·c_in,f), lag + overflow rows | `MISSING_NEEDS_MATLAB_PORT` | none | — | uses the Python `utility.py` parmest + FIM workflow → no `.m`, but no Python wrapper that emits this exact table | `utility.py:report_parameters` + `DATA2_model_demo.ipynb` |
| 25 | Fig. 6 (A,B,C / D,E,F) | Lag (top row) + overflow (bottom row) — mass, concentration, mass-in-stirred-cell time series for 270511.123, 270511.423 | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites pre-existing PNGs: `mass-dat270511.123.png`, `concentration-dat270511.123.png`, `stirc_mass-dat270511.123.png`, etc.) | requires those panel PNGs already on disk in `paper_artifacts/data2/notebook_figures/` — they come from `DATA2_model_demo.ipynb` which reads `.mat` data | `utility.py:plot_sim_comparison` | `DATA2_model_demo.ipynb` |
| 26 | **Table 3** | AIC ranking, concentrating regime | `MISSING_NEEDS_MATLAB_PORT` | none | — | (no `.m` — paper-generated) | `DATA2_model_demo.ipynb` (AIC analysis cells) |
| 27 | **Table 4** | AIC ranking, diluting regime | `MISSING_NEEDS_MATLAB_PORT` | none | — | — | as above |
| 28 | Fig. 7 (A,B,C) | Convection-diffusion empirical model fitting: c_in,f & c_h / Jw / Js panels | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites `Js_predict0.png`, `Jw_predict.png`, `Js_predict1.png`) | the underlying panels come from `DATA2_visualization.ipynb` notebook output | — | `DATA2_visualization.ipynb` + `utility.py:plot_model_predictions` |
| 29 | Fig. 8 (A,B,C) | MSE + h0 + h1 vs Peclet (partition sensitivity) | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_partition_sensitivity` (line 9467) + `run_data2_publication_figures` alias `figure_8.png ← partition_sensitivity.png` | reads notebook-generated `partition_sensitivity.png` | — | `DATA2_visualization.ipynb` |
| 30 | **Table 5** | Model predictions for lag/overflow × M1-M4 (original vs truncated) | `MISSING_NEEDS_MATLAB_PORT` | none | — | — | `DATA2_model_demo.ipynb` |
| 31 | **Table 6** | Information improvements from FIM | `MISSING_NEEDS_MATLAB_PORT` | none | — | — | `run_DATA2_model_variations.py` + `utility.py:calc_FIM` |
| 32 | Fig. 9 (A, B) | Information improvement bar plot + weighted residual box plots | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites `startup_barplot.png` + `concentrating_residuals_boxplot.png`) | requires those panel PNGs from `DATA2_model_demo.ipynb` / `DATA2_visualization.ipynb` outputs | — | as above |

---

## DATA2_SI

| # | Item | Description | Verdict | Python function | Dependency | MATLAB source | Python reference |
|---|---|---|---|---|---|---|---|
| 33 | AIC derivations + Donnan eqs (S1-S29) | Equations only, no figures | `MISSING_OUT_OF_SCOPE` | — | — | — | — |
| 34 | **Table S1** | Estimated model parameters for variations of solute transport eq (constant B, B_i const, polynomial B with I-sets) | `MISSING_NEEDS_MATLAB_PORT` | none | — | (Python `utility.py` parmest output, no formal table function) | `utility.py` + `DATA2_model_demo.ipynb` |
| 35 | Fig. S1 (A, B) | Calibration curves for CN0359 + EC EvalBoard LSF boards | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_calibration_plots` (line 9308) + `figure_s1.png ← calib_curve.png` alias | reads `data_library/conductivity_calibration{1,2}.csv` | — | `DATA2_model_demo.ipynb` |
| 36 | Fig. S2 (A.1, A.2 × 2 panels) | Additional diafiltration experiments — lag, 60 psi, 15/120 mM | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites `mass-dat270611.121.png`, `concentration-dat270611.121.png`, etc.) | panel PNGs from `DATA2_model_demo.ipynb` reading `.mat` | — | as above |
| 37 | Fig. S3 (B.1, B.2 × 2) | Diafiltration at 40 psi + 20 psi, 15/120 mM | `PRESENT_WITH_MAT_DEPENDENCY` | same composite path (`mass-dat270511.221.png`, `270511.321.png`) | as above | — | as above |
| 38 | Fig. S4 (C.1, C.2 × 2) | Overflow mode experiments at 60 psi | `PRESENT_WITH_MAT_DEPENDENCY` | same (`mass-dat270511.421.png`, `270511.921.png`) | as above | — | as above |
| 39 | Fig. S5 (D.1, D.2 × 2) | Varying KCl conditions (50/100, 10/50 mM) | `PRESENT_WITH_MAT_DEPENDENCY` | same (`270511.521`, `270511.621`) | as above | — | as above |
| 40 | Fig. S6 (E.1, E.2 × 2) | Diluting regime experiments (30/1, 100/12 mM) | `PRESENT_WITH_MAT_DEPENDENCY` | same (`270511.721`, `270511.821`) | as above | — | as above |
| 41 | Fig. S7 (A, B) | Discrete B per vial + Js/Jw vs interface concentration (multi-experiment overlay) | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites `Bpervial.png` + `Js_Jw_cin.png`) | underlying panels from `run_pre_B_dependence.py` / notebook outputs | — | `run_pre_B_dependence.py` |
| 42 | Fig. S8 (A, B, C) | Convection-diffusion model fit to empirical simulation: c_in,f & c_h / Jw / Js | `PRESENT_WITH_MAT_DEPENDENCY` | `run_data2_publication_figures` (composites `Js_predict0.png` + `Jw_predict.png` + `Js_predict.png`) | as above | — | `DATA2_visualization.ipynb` |
| 43 | Fig. S9 | Zeta potential vs pH for NF270 (single scatter with shaded uncertainty) | `MISSING_NEEDS_MATLAB_PORT` | none (zeta potential is an external instrument measurement; no model code generates it) | — | — | external; could be a static asset |

---

## Verdict counts (audit summary block)

```
DATA1_main:   PRESENT=0  PRESENT_WITH_MAT=4  PRESENT_WITH_V24=1  NEEDS_PORT=2  OUT_OF_SCOPE=1   (total 8)
DATA1_SI:     PRESENT=0  PRESENT_WITH_MAT=5                       NEEDS_PORT=3  OUT_OF_SCOPE=1   (total 9)
DATA2_main:   PRESENT=0  PRESENT_WITH_MAT=4                       NEEDS_PORT=5  OUT_OF_SCOPE=6   (total 15)
DATA2_SI:     PRESENT=0  PRESENT_WITH_MAT=8                       NEEDS_PORT=2  OUT_OF_SCOPE=1   (total 11)
                ────────────────────────────────────────────────────────────────────────────────────
TOTAL:        PRESENT=0  PRESENT_WITH_MAT=21 PRESENT_WITH_V24=1   NEEDS_PORT=12 OUT_OF_SCOPE=9   (total 43)
```

(Total 43 vs 41 named items: Figs 5 and 6 in DATA1_main each have an alternative direct-contour path that I counted separately because it carries a different verdict.)

## Distinct MATLAB scripts (or missing-from-repo MATLAB scripts) that block figures

Sorted by figure-count impact:

| Rank | MATLAB source | Currently in `legacy/data1_matlab/functions/`? | Figures it blocks |
|---|---|---|---|
| 1 | `calc_contour_2d.m` + `plot_contour.m` | **YES** (both present) | DATA1 Figs. 5, 6, S3, S4, S5, S6 (= 6 figures, 60+ panels) |
| 2 | `sigma_sensitivity.m` + `heatmap_sigma_sensitivity_filtration.m` + `heatmap_sigma_sensitivity_diafiltration.m` | only **`sigma_sensitivity.m`** is present; the two `heatmap_*.m` ones are NOT in this checkout (referenced in README only) | DATA1 Fig. 4 (= 1 figure, 6 panels) |
| 3 | `run_data_analysis.m` (the WLS fit + Table 1/Table 2 + `experiment space/*.csv` emission branch) | **YES** | DATA1 Figs. 2, 3, S2 + Tables 1, 2 (= 3 figures + 2 tables) |
| 4 | `doe_heatmap_filtration.m` + `doe_heatmap_diafiltration.m` | **NO** — described in README but `.m` files are missing | DATA1 Fig. S7 (= 1 figure, 12 panels) |
| 5 | `load_data.m` + `model_config.m` + `update_model_config.m` (loader chain) | **YES** | Indirect — all DATA1 figures inherit data shape from these |
| 6 | `calc_FIM.m` + `calc_ind_objectives.m` | **YES** | DATA1 Table 1, Table 2 objectives column; DATA2 Tables 3-6 |

## Critical finding — the V24 leak

`run_data1_direct_contour_branch` (`refactored_ucb_library.py:7476`) is the only path that computes DATA1 contours "from scratch" without reading MATLAB outputs. It does so by importing from `UnifiedFramework.DATA3.ExperimentalDataAnalysis.UnifiedCode.unified_codebase_library` (the V24 stack outside `refactored_codes_v1`). This means:

- Even the "clean" Python contour path is NOT self-contained within `refactored_codes_v1`.
- Deleting the V24 UnifiedCode tree would break DATA1 direct contours.
- To make `refactored_codes_v1` a self-contained authority, `evaluate_data1_paper_contour_objectives_v24` and `simulate_data1_vialwise_trajectories_v24` would need to be ported INTO `refactored_ucb_library.py`.

## Critical finding — DATA2's "publication figures" path is a compositor, not a generator

`run_data2_publication_figures` (line 8489) is what `materialize_all(campaign="DATA2")` calls for the main paper figures. It does **NOT** compute figures — it COMPOSITES pre-existing PNG panels found at known filenames under `paper_artifacts/data2/notebook_figures/`. Those panels are produced by:

- `DATA2_model_demo.ipynb` (running cells that read `.mat` data through `utility.py`)
- `DATA2_visualization.ipynb`
- `run_DATA2_model_variations.py`
- `run_pre_B_dependence.py`

So even "from scratch" DATA2 figure generation today requires running those legacy entrypoints (or having their cached PNGs on disk). The pytest-suite "can the current code reproduce the paper" must include either porting those panel generators into `refactored_ucb_library.py`, OR running the legacy entrypoints as part of materialize_all.

## Recommended port order (for the operator's review)

Each row below would spawn its own per-port proposal under `docs/ports/PORT_<name>.md` — **only with explicit operator approval per port**.

| Order | Port | What it unlocks | Estimated lines | Risk |
|---|---|---|---|---|
| 1 | `calc_contour_2d.m` → Python `compute_contour_2d_native(...)` in library | DATA1 Figs. 5, 6, S3, S4, S5, S6 (largest impact) | ~200 LOC | LOW — already have `calc_contour_2d_py` (line ~5800), need to verify identity |
| 2 | `sigma_sensitivity.m` → Python `compute_sigma_sensitivity_native(...)` | DATA1 Fig. 4 | ~100 LOC | LOW — forward sim only |
| 3 | `run_data_analysis.m` table-emission branch → Python `build_data1_parameter_tables(...)` | Tables 1, 2 | ~150 LOC | MEDIUM — model selection + objective computation |
| 4 | V24 lift — port `evaluate_data1_paper_contour_objectives_v24` + `simulate_data1_vialwise_trajectories_v24` into library | Remove V24 dependency from `run_data1_direct_contour_branch` | ~300 LOC | MEDIUM — affects the only "from-scratch" contour path |
| 5 | DATA2 panel generators (`mass-dat*.png`, `concentration-dat*.png`, `stirc_mass-dat*.png`) | DATA2 Fig. 6 + SI S2-S6 (8 figures, ~24 panels) | ~250 LOC | MEDIUM — depends on notebook flow being formalized |
| 6 | DATA2 table generators (Tables 2, 3, 4, 5, 6) | DATA2 main paper tables | ~200 LOC | MEDIUM — requires AIC + FIM result formatting |
| 7 | `doe_heatmap_filtration.m` + `doe_heatmap_diafiltration.m` reconstruction | DATA1 Fig. S7 (MBDoE optimality) | ~250 LOC | HIGH — MATLAB `.m` files not in repo; need to reconstruct from PSE2021 paper |

## What this audit does NOT do (reaffirmed from prompt)

- Does not propose any code change inline.
- Does not modify `refactored_codes_v1/`, `legacy/`, `UnifiedFramework/`, or `target_notebook_source_map.csv`.
- Does not run `materialize_all`.
- Does not start any port — each port is its own discussion + proposal + approval cycle per the operator's standing instruction:
  > "Do not change the DATA1 and DATA2 workflow without talking to me clearly about it and what you want to change."

## Next operator decision points

1. Review the per-row verdicts. Disagreements?
2. Decide port priority (recommended order above is a draft).
3. For each port the operator approves, I'll author `docs/ports/PORT_<name>.md` with:
   - What the MATLAB script does (math + I/O)
   - Proposed Python function signature in `refactored_ucb_library.py`
   - Acceptance test (Python output reproduces the MATLAB-generated `.mat`/CSV byte-for-byte or within tolerance)
   - Per-port discussion with the operator before any code lands.
