# Prompt for Sonnet — CaCl2 (2:1) Donnan regression + adsorption-hypothesis diagnostics (report §3)

You are extending the DATA3 analysis to **single-salt \ce{CaCl2}** (a 2:1 salt) on
NF270. Work in `/Users/adowling/DowlingLab/Membranes/data3_regression`, in the
`data3-regression` conda env (`conda activate data3-regression`). Run everything
there and verify it executes cleanly.

> **Run this after the NaCl figure-cleanup task** (`docs/prompts/nacl_figure_cleanup.md`),
> which finalizes `analysis/plot_style.py`. Import and use that shared style so
> the CaCl2 figures match the rest of the report automatically.

## Why this is different from NaCl (read first)
The PI's **working hypothesis** is that **\ce{Ca^{2+}} adsorbs onto the membrane
and changes its effective fixed charge $|\chi|$ in a concentration-dependent
way.** The current model treats $|\chi|$ as a single constant and does **not**
represent adsorption. So we *expect the CaCl2 fits to be worse than NaCl*, and the
**main deliverable is a diagnostic test of that hypothesis**, not just a fit. The
signatures to look for: poor overall fit, and **systematic residual structure vs.
concentration** / an **apparent $|\chi|$ that drifts with concentration**.

## Datasets (all processed; use `data/BoE Analysis.xlsx`)
> **Important:** read from **`data/BoE Analysis.xlsx`** — the 27-sheet superset.
> Do **not** use `data/Rejection_Analysis.xlsx`; it is missing the 06.27.26 run.

All three use the standard 15-column schema (headers say "CaCl2"), header row 1,
numeric data in cols **B**=`c_int` (CaCl2 interfacial conc, mM), **H**=`c_p`
(CaCl2 permeate interfacial conc, mM), **J**=`J_w` (m³/m²/s) — position-readable
exactly like NaCl.

| Membrane | Sheet | Valid rows | c_int range (mM) |
|---|---|---|---|
| MC3 | `NF270_MC3 07.11.24_SCaCl2`  | 614  | 1.0 – 47   |
| MC5 | `NF270_MC5 05.27.26_CaCl2`   | 1645 | 1.6 – 108  |
| MC5 | `NF270_MC5 06.27.26_CaCl2`   | 1885 | 1.3 – 103  |

(The MC5 pair 05.27.26 vs 06.27.26 is a near-repeat — useful reproducibility
check. Raw-only CaCl2 runs exist but are out of scope here.)

## Model (2:1 co-ion cubic) — WORKING ASSUMPTIONS (Opus)
These choices mirror the NaCl treatment and are sufficient for the hypothesis
test; they connect to the open items in report §1.2 ("Modeling questions") and can
be refined later. State them in your README.

- **Co-ion = \ce{Cl-}**; its solution concentration is set by stoichiometry:
  `c_co_s = 2 * c_CaCl2` (2 Cl⁻ per CaCl2). Use this for both the feed-side
  (`c_int`) and permeate (`c_p`) interfacial values.
- **Donnan co-ion cubic** (report Eq. 21co), lumped partition constant `K` (≥0)
  absorbing `δ* δ°_co` and constant stoichiometric factors:
  `(c_co_m)^3 + |χ|(c_co_m)^2 − K*(c_co_s)^3 = 0`.
  There is exactly **one positive real root** (Descartes: sign pattern +,+,0,−);
  select it via `numpy.roots([1, |χ|, 0, -K*c_co_s**3])` → the unique real, >0 root.
- **Transformed response** (Cl⁻ flux form, analogous to corrected NaCl):
  `deltaC_m = J_s_Cl * l / D_m`, with `J_s_Cl = 2 * c_p * J_w` (Cl⁻ flux = 2× the
  CaCl2 molar flux), `l = 80e-9` m, and the **2:1 ambipolar** membrane diffusivity
  `D_m = 3*D_Ca_m*D_Cl_m/(2*D_Ca_m + D_Cl_m)` using membrane values
  `D_Ca_m = 0.79e-12`, `D_Cl_m = 2.03e-12` m²/s (the same ~1000× solution scaling
  Bill used for NaCl; solution values are `D_Ca=0.79e-9`, `D_Cl=2.03e-9`).
- **Prediction** = `c_co_m(c_int) − c_co_m(c_p)` (solve the cubic at each face).
  Fit `(|χ|, K)` to `deltaC_m` by `scipy.optimize.least_squares`. Suggested start
  `(|χ|, K) = (44, 0.3)`; a lumped constant prefactor mismatch would bias absolute
  `|χ|`/`K` but **not** the residual-shape diagnostic, which is the point here.

## Per-dataset analysis
For each of the 3 datasets:
1. Fit `(|χ|, K)`; report estimate, Wald SE, SSE, and a scale-free fit metric
   (e.g. RMSE / range, or pseudo-$R^2$). Multi-start (100 starts, method `lm`/`trf`
   — avoid `dogbox`, which trapped at the boundary for NaCl) to confirm the global
   optimum and report % global.
2. 95% nonlinear (F-test) confidence region + profile CIs (the 2:1 fit may be less
   linear/identifiable than NaCl — note if the region departs from the ellipse).

## Adsorption-hypothesis diagnostics (the core deliverable)
3. **Residuals vs. concentration:** plot fit residuals (`deltaC_m − prediction`)
   against `c_int`. Systematic (non-random) structure supports concentration-
   dependent physics the model misses. Quantify (e.g. sign runs, or a smoothed
   trend).
4. **Apparent $|\chi|$ vs. concentration:** in a sliding window over `c_int` (or by
   binning), refit `|χ|` locally (K fixed at global, or refit both) and plot the
   local effective `|χ|` vs. concentration. A monotonic drift is the direct
   signature of concentration-dependent charge (adsorption).
5. Compare across datasets (MC3 vs MC5; the MC5 05.27 vs 06.27 near-repeat) and
   contrast fit quality with NaCl (which fit well).

## Figures — use the finalized shared style `analysis/plot_style.py`
The NaCl figure-cleanup task finalized `plot_style.py`; **follow its conventions
exactly** so CaCl2 figures match the report with no rework (see its docstring):
- Call `plot_style.apply_style()`; set each figure's `figsize` **width** via
  `plot_style.fig_width("full")` (≈6.50 in) or `fig_width("half")` (≈3.19 in) so
  text renders at its true ~11pt size on the page (do **not** use ad hoc
  `figsize=(4,4)`).
- **Legends go below the axes** (`loc="upper center"`, `bbox_to_anchor` just below,
  `frameon=False`), never `loc="best"`/inside. For **multi-panel figures with the
  same legend in every panel, use ONE shared legend** below the whole figure (as
  `nacl_multistart`/`nacl_all_fits` now do) — do not repeat per-panel legends.
- One consistent colorblind-palette color per dataset; units in bracketed axis
  labels; save **PNG + PDF** with `bbox_inches='tight'`.

Produce into `docs/reports/figures/` (you will embed these into report §3 and the
CaCl2 appendix — see "Report presentation" below; save at the widths below):
- `cacl2_fits.{png,pdf}` — **full-width** multi-panel (one panel per dataset),
  shared legend below (like `nacl_all_fits`).
- `cacl2_residuals_vs_conc.{png,pdf}` — **half-width**; **the key diagnostic**.
- `cacl2_effective_chi_vs_conc.{png,pdf}` — **half-width**; local `|χ|` vs. `c_int`.
  (These two half-width diagnostics will sit side-by-side in the report.)
- `cacl2_estimates.{png,pdf}` — **half-width**; `(|χ|, K)` with CIs across datasets.

## Report presentation — match the NaCl standard (important to the PI)
This report is a **comprehensive preliminary handoff/discussion document, not a
summarized paper.** Favor completeness and transparency over brevity — it is fine
(desired) for it to be long and to **show every dataset individually.** Write up
CaCl2 to the *same standard as NaCl*: populate the CaCl2 section (`\S`~`sec:cacl2`,
currently a TBD stub) and add a CaCl2 per-dataset appendix, mirroring the NaCl
structure (`sec:naclAll` §2.4 + `sec:naclappendix` §6) so the two salts read
consistently. Specifically:
1. **What you did** — a short intro: the 2:1 cubic model + working assumptions, the
   datasets, and why CaCl2 is the adsorption-hypothesis test.
2. **Inventory / provenance (kept current)** — state exactly which CaCl2 runs were
   analyzed. If only the 3 processed runs, say so explicitly and list the raw-only
   CaCl2 runs (MC2 05.07.24, MC3 07.11.24\_CaCl2, MC3 07.12.24\_S2CaCl2, and the
   mislabeled MC5 05.26.26) as not-yet-processed — do not leave a stale/ambiguous
   status.
3. **Unified results table** — every CaCl2 dataset regressed in ONE table with an
   explicit count and a provenance column, matching `tab:naclall`'s style
   (`|χ|`, `K`, 95% CI, SSE, fit metric, multi-start %).
4. **What you found** — the adsorption-hypothesis verdict in prose (is the fit
   materially worse than NaCl? systematic residual-vs-concentration structure?
   effective-`|χ|` drift, and in which direction? consistent across datasets?),
   plus a reproducibility note on the MC5 05.27 vs. 06.27 near-repeat.
5. **Show the results — per-dataset appendix** — one block per CaCl2 dataset with
   its fit-vs-data, the **residual-vs-`c_int`** and **effective-`|χ|`-vs-`c_int`**
   diagnostics, and the parameter/CI summary, so a collaborator can dig into each
   experiment (mirrors `sec:naclappendix`).
The section must let a reader answer, for CaCl2, *what did you do, what did you
find, and show me the results so I can dig in myself* — without opening the code.
Figure quality: **no vertical compression**; if per-dataset confidence regions span
very different scales, use the **zoom-inset** approach of
`nacl_all_confidence_regions`; distinct color per dataset from the extended
`plot_style.PALETTE`; and **rasterize (`pdftoppm`) and eyeball every figure**
(rendering issues have bitten us before).

## Deliverables & guardrails
- `analysis/cacl2/regress_cacl2.py`, `analysis/cacl2/README.md`.
- Reuse `analysis/step2_nacl/regress_nacl.py` helpers where sensible (data extent,
  fit utilities, `plot_style`); the model function is new (cubic root). You may
  make a **backward-compatible** change to `step2.load_sheet` to accept an optional
  `data_file=` (and existing `sheet_name=`) argument, defaulting to the current
  values — verify the NaCl scripts still reproduce their committed numbers after.
- Keep new analysis files under `analysis/`; figures under `docs/reports/figures/`.
- You **may** edit `docs/reports/main.tex` **only** within the CaCl2 section
  (`sec:cacl2`) and by adding a CaCl2 per-dataset appendix subsection; compile with
  `latexmk -pdf` (confirm **0 undefined citations, 0 overfull hboxes**), then
  `latexmk -c`. **Do NOT modify** any other part of `main.tex` (the §1 model math /
  notation, §2 NaCl prose/tables/numbers, `\author`/preamble), nor `refs.bib`,
  `data/`, or `prior_analysis/`. Do not alter any NaCl results.

## Report back (paste-ready for §3)
1. Per-dataset table: `|χ|`, `K` (95% CI), SSE, scale-free fit metric, multi-start
   % global.
2. **Verdict on the adsorption hypothesis:** is the fit materially worse than NaCl?
   Do residuals show systematic concentration structure? Does the local effective
   `|χ|` drift with concentration, and in which direction? Is the effect consistent
   across the 3 datasets?
3. A section-by-section summary of the CaCl2 write-up you added (`sec:cacl2` +
   appendix): the tables/figures created and where.
4. Confirmation: everything ran in `data3-regression`; all figures written
   (PNG+PDF) and **rasterize-verified**; the report compiles clean (0 undefined,
   0 overfull) with the new page count; NaCl numbers unchanged if `load_sheet` was
   touched; and only `sec:cacl2` + the new appendix were edited in `main.tex`
   (§1/§2/`refs.bib` untouched).
