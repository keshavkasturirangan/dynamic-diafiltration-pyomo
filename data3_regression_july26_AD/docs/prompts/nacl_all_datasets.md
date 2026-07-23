# Prompt for Sonnet — regression across all NaCl datasets (report §2.4)

You are extending the DATA3 membrane-regression analysis to **all available
single-salt NaCl datasets**. Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`.
Use the documented conda env: `conda activate data3-regression` (Python 3.11,
numpy/scipy/pandas/matplotlib, pyomo 6.10, ipopt). Run everything in that env and
verify scripts execute cleanly.

> **Sequencing / coordination (important).** A separate session is concurrently
> adding a bibliography to `docs/reports/main.tex` (creating `refs.bib`). **Run
> this task only after that bibliography edit is committed.** In this task you
> **must not edit `docs/reports/main.tex` or `refs.bib`, and must not recompile
> the report** (`latexmk`) — you only (re)generate image files under
> `docs/reports/figures/`. The report will be recompiled separately once both
> streams of work are merged.

## Goal
Fit the NaCl (1:1) Donnan model to every processed NaCl experiment and assess
**reproducibility** of the membrane parameters `|χ|` (fixed charge, mM) and `δ*`
(partition factor) across (a) different membrane cuts (MC2 vs MC5) and (b) a
repeat run on the same cut (MC5 `NaCl` vs `NaCl (2)`). This populates §2.4 ("Workflow
applied to all available NaCl datasets") of `docs/reports/main.tex`.

## Datasets (all processed, standard schema, ready now)
All three live in `data/Rejection_Analysis.xlsx` (and identically in
`data/BoE Analysis.xlsx`), 15-column schema, header row 1, data from row 2.
Columns (1-indexed): **B** = `c_int` (interfacial conc, mM), **H** = `c_p`
(permeate interfacial conc, mM), **J** = `J_w` (m³/m²/s).

| Membrane | Sheet name | Valid rows |
|---|---|---|
| MC2 | `NF270_MC2 05.07.24_NaCl`     | 482 |
| MC5 | `NF270_MC5 05.27.26_NaCl`     | 747 |
| MC5 | `NF270_MC5 05.27.26_NaCl (2)` | 747 (analyzed in §2.1–2.2) |

(The 6 raw-only NaCl runs in `data/NF270_MC{2,3,4,5}.xlsx` are **out of scope
here** — they need the interfacial-conc/Jw preprocessing first. Note them as
future work only.)

## Model and code reuse
- Reuse the model and fit machinery already written:
  `analysis/step2_nacl/regress_nacl.py` (functions `build_deltaC_m`, `f_matlab`,
  `f_corrected`, `fit_model`) and the confidence-region / profile utilities in
  `analysis/step3_nacl_uncertainty/scipy_uncertainty.py` (F-test threshold, SSE
  surface, profile likelihood). **Do not duplicate that logic** — import it.
- `step2.load_sheet()` currently hardcodes the sheet name; generalize it to accept
  an optional `sheet_name` argument **defaulting to the current value** (a
  backward-compatible change so `step3` still works), or add a thin loader in your
  new script. Positional columns B/H/J apply to all three sheets.
- **Primary model = `f_corrected`** (the physically correct ½-form), consistent
  with §2.1–2.2. Also report the uncorrected `f` estimates in a secondary column
  for continuity with Bill.
- **Windowing:** for comparability across datasets (which lack Bill's ad hoc
  57–748 window), use the **full valid data range** for each as the primary fit,
  and additionally report a **startup-trimmed** fit (drop the low-concentration
  transient, e.g. first ~7% of rows) as a sensitivity check. State the rule you use.

## Analysis per dataset
For each of the 3 datasets, using `f_corrected`:
1. Fit `(|χ|, δ*)` with `scipy.optimize.least_squares`; report estimate, Wald SE,
   SSE, n.
2. Confirm optimizer reliability with a short multi-start (e.g. 100 starts over
   `|χ|∈[1,200]`, `δ*∈[1e-3,10]` log-uniform), using method `lm` or `trf`
   (recall `dogbox` can get trapped at `|χ|→0`); report % reaching the global
   optimum.
3. Compute the 95% nonlinear (F-test) confidence region and profile-likelihood CIs.

Then **compare across datasets**: are `|χ|` and `δ*` consistent within their
confidence regions? Do the two MC5 repeats agree? Does MC2 differ from MC5?
Compare `|χ|` to the independently-used `χ≈44 mM`.

## Publication-quality figures + consistent style across ALL figures (required)
Follow the guidelines at
<https://ndcbe.github.io/data-and-computing/notebooks/01/Publication-Quality-Figures.html>.
**Create a shared style helper `analysis/plot_style.py`** (a function like
`apply_style()` that sets `matplotlib.rcParams`, plus a fixed colorblind-friendly
palette and a consistent color-per-dataset mapping). **Every figure in the project
must use it** — not just the new ones. Concretely, the helper sets:
- Figure size 4×4 in (or 4×6.4 for tall/multi-panel); **dpi 300** (PNG), vector PDF too.
- Line width **3**, marker size **8**.
- Tick label font size **15**, ticks `direction="in"`, `top=True, right=True`.
- Axis labels font size **16**, `fontweight='bold'`, **units in brackets**, e.g.
  `$|\chi|$ [mM]`, `$\delta^*$ [–]`.
- Legend font size 10–14, `frameon=False`; external placement via
  `bbox_to_anchor` with `borderaxespad=0` when it declutters.
- **Colorblind-friendly palette** (e.g. seaborn colorblind or Okabe–Ito hex list);
  one consistent color per dataset across all figures.
- Save **both** `.png` (dpi 300) and `.pdf`, always `bbox_inches='tight'`.
- Multi-panel: use `label_outer()` to drop redundant labels; `fill_between` for
  any uncertainty bands; `set_aspect` for square axes where sensible.

Produce (into `docs/reports/figures/`, PNG+PDF):
- `nacl_all_estimates.{png,pdf}` — point estimates of `|χ|` and `δ*` with 95% CI
  error bars across the 3 datasets (e.g. two stacked panels or a forest-style plot).
- `nacl_all_confidence_regions.{png,pdf}` — the three 95% nonlinear confidence
  regions overlaid on the `(|χ|, δ*)` plane, one color per dataset, optima marked.
- `nacl_all_fits.{png,pdf}` — multi-panel `f_corrected` fit vs. `Δc_m` data, one
  panel per dataset (shared axes, `label_outer()`).

### Restyle the existing Step 2 & Step 3 figures for consistency
So the whole report is visually uniform, regenerate the earlier figures through
the same `plot_style`:
- Add `import plot_style; plot_style.apply_style()` (styling only — **no logic
  change**) near the top of `analysis/step2_nacl/regress_nacl.py`,
  `analysis/step3_nacl_uncertainty/scipy_uncertainty.py`, and
  `analysis/step3_nacl_uncertainty/parmest_uncertainty.py`, and make each also
  save a `.pdf` beside every `.png` (via `bbox_inches='tight'`). Ensure axis
  labels carry units in brackets (e.g. `$|\chi|$ [mM]`, `$\delta^*$ [–]`).
- Re-run all three scripts so every PNG/PDF in `docs/reports/figures/` is
  regenerated with the shared style. (Note: `parmest_uncertainty.py` takes
  ~6 min due to the confidence-region grid + bootstrap; that is expected.)
- The existing figure **filenames must stay the same** (the report references them
  by name); you are only restyling their contents and adding PDF siblings.

## Deliverables
- `analysis/nacl_all_datasets/regress_all_nacl.py` — loops the 3 datasets, reuses
  step2/step3 code, prints a summary table, writes the figures.
- `analysis/plot_style.py` — the shared rcParams helper.
- `analysis/nacl_all_datasets/README.md` — how to run + results summary table.
- Updated `analysis/step2_nacl/regress_nacl.py` and both
  `analysis/step3_nacl_uncertainty/*.py` — with the `plot_style` import and PDF
  output added (styling only), plus the backward-compatible `load_sheet(sheet_name=...)`
  change in step2.
- All figures under `docs/reports/figures/` regenerated (PNG **and** PDF) with the
  shared style — the new `nacl_all_*` figures and the restyled Step 2/3 figures.
- **Allowed edits:** new files under `analysis/`; the `plot_style` import + PDF
  saving + `load_sheet` argument in step2/step3 (no change to their numerical
  logic). **Do NOT modify** `prior_analysis/`, `data/`, `docs/reports/main.tex`,
  or `docs/reports/refs.bib`, and **do NOT run `latexmk`** (see the sequencing note
  at the top). Verify numerical outputs are unchanged from the committed values
  after your edits (the estimates must not move — this is a styling pass).

## Report back (paste-ready for §2.4)
1. A table: per dataset — membrane, n, `|χ| ± 95% CI`, `δ* ± 95% CI`, SSE, and
   multi-start % global (corrected model; full-range primary).
2. The startup-trimmed sensitivity (how much estimates move).
3. A 2–3 sentence reproducibility verdict: do the parameters agree across cuts and
   across the MC5 repeat? Is `|χ|` near 44 mM?
4. Confirmation all scripts ran in the `data3-regression` env and all figures
   (PNG+PDF) were written; note the filenames. Confirm the restyled Step 2/3
   scripts reproduce the **previously committed estimates unchanged** (this is a
   styling-only pass) and that you did not touch `main.tex`/`refs.bib` or run
   `latexmk`.
