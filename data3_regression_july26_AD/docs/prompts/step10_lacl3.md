# Prompt for Sonnet — Step 10: LaCl3 (3:1) Donnan regression via nonlinear constraints + adsorption diagnostics (report §4)

You are extending the DATA3 analysis to **single-salt \ce{LaCl3}** (a 3:1 salt) on
NF270 — the last of the single-salt trilogy. Work in
`/Users/adowling/DowlingLab/Membranes/data3_regression`, env `data3-regression`.
Populate the report's **`sec:lacl3` (§4, currently a TBD stub)** and add a LaCl3
per-dataset appendix, to the **same comprehensive standard as CaCl2** (§3 /
`sec:cacl2appendix`). Explain the math inline (exploratory report). Do not change
NaCl/CaCl2 results or the model math.

> **Run after Step 9 has been committed** (both edit `main.tex`; avoid concurrent
> edits). Use the finalized `analysis/plot_style.py` (boxed legends, bold math via
> `\boldsymbol{}`) — the figures must comply.

## What's distinctive about Step 10 (the methodological point)
The 3:1 co-ion Donnan equation is a **quartic with no compact closed form**. The
plan is to solve it as a **nonlinear constraint in Pyomo/ParmEst** rather than
coding a closed-form root — this is the approach that will generalize to the
multicomponent case (Step 14), where no closed form exists. So do the fit **two
ways and cross-check** (as in Step 3's scipy-vs-ParmEst validation):
- **(A) scipy workhorse:** solve the quartic numerically per point (e.g.\
  `numpy.roots` → unique real positive root, or a bracketed Newton) and fit by
  `scipy.optimize.least_squares` — fast, and drives the diagnostics.
- **(B) ParmEst / Pyomo with the quartic as a nonlinear constraint:** introduce the
  membrane co-ion concentration as a Pyomo `Var` constrained by the quartic
  (positivity-bounded to select the physical root), estimate $(|\chi|,K)$ with
  ParmEst/Ipopt. Confirm it matches (A). Note this is the formulation that scales
  to Step 14.

## Dataset
- **One processed LaCl3 run:** `NF270_MC2 05.21.24_LaCl3` in `data/BoE Analysis.xlsx`
  (also in `Rejection_Analysis.xlsx`), 623 valid rows, standard 15-col schema
  (B=`c_int`, H=`c_p`, J=`J_w`). **Narrow concentration range, c_int ≈ 2–28 mM.**
- A raw-only **MC4 `07.11.24_SLaCl3`** exists (raw campaign workbook) but is **out
  of scope** here (would need the preprocessing pipeline; note it as available for
  future work, like the raw-only NaCl runs).

**Anticipate two headwinds (state them):** (i) the narrow 2–28 mM range gives
**weak |χ| identifiability** (the same range effect flagged for MC2-NaCl and MC3
07.09.24); (ii) La³⁺ rejection is low, so the fit may be ill-conditioned — handle
gracefully (as the NaCl multi-dataset code does) and report honestly if so.

## Model (3:1 co-ion quartic) — WORKING ASSUMPTIONS (Opus; mirror the CaCl2 treatment)
- **Co-ion = \ce{Cl-}**, stoichiometry `c_co_s = 3 * c_LaCl3` (3 Cl⁻ per LaCl3), for
  both `c_int` and `c_p`.
- **Donnan co-ion quartic** (report Eq. 31co), lumped constant `K` (≥0) absorbing
  `δ*(δ°_co)^2` and stoichiometric factors:
  `(c_co_m)^4 + |χ|(c_co_m)^3 − K*(c_co_s)^4 = 0`.
  Descartes ⇒ exactly one positive real root; select it.
- **Transformed response:** `deltaC_m = J_s_Cl * l / D_m`, with `J_s_Cl = 3*c_p*J_w`
  (Cl⁻ flux = 3× the LaCl3 molar flux), `l = 80e-9` m, and the **3:1 ambipolar**
  membrane diffusivity `D_m = 4*D_La_m*D_Cl_m/(3*D_La_m + D_Cl_m)` with
  `D_La_m = 0.626e-12`, `D_Cl_m = 2.03e-12` m²/s (the ~1000× solution scaling used
  for NaCl/CaCl2; solution `D_La = 6.26e-10`).
- **Prediction** = `c_co_m(c_int) − c_co_m(c_p)`. Fit `(|χ|, K)` **bounded (≥0)**
  (as for CaCl2 — χ enters the quartic asymmetrically, so an unconstrained fit can
  wander negative). Suggested start `(|χ|, K) = (44, 0.01)`; a lumped prefactor
  mismatch biases absolute `|χ|`/`K` but not the residual-shape diagnostic.

## Analysis (mirror CaCl2 §3)
1. Fit `(|χ|, K)` [Method A]; report estimate, Wald + profile CI, SSE, a scale-free
   fit metric (pseudo-R²), and multi-start % global (`lm`/`trf`, not `dogbox`).
2. Cross-check with Method B (ParmEst nonlinear-constraint) — confirm agreement.
3. **Adsorption diagnostics** (the core deliverable, as for CaCl2 — La³⁺ is
   expected to adsorb even more strongly than Ca²⁺): residuals vs. `c_int`
   (sign-runs test) and sliding-window effective `|χ|` vs. `c_int`. Over the narrow
   range with few reliable windows, the effective-χ diagnostic may be weak — flag
   unreliable windows and don't over-interpret.
4. Compare qualitatively to NaCl (well-identified) and CaCl2 (χ→0, adsorption): does
   LaCl3 also drive |χ|→0 / show residual structure? Is the 3:1 case even more
   extreme?

## Figures (finalized `plot_style`, PNG+PDF, rasterize-verify)
Into `docs/reports/figures/` — `lacl3_fit.{png,pdf}` (fit vs. `deltaC_m`),
`lacl3_residuals_vs_conc.{png,pdf}`, `lacl3_effective_chi_vs_conc.{png,pdf}`. With a
single dataset, a compact multi-panel appendix figure `lacl3_appendix.{png,pdf}`
(fit, residual, effective-χ, parameter summary) is appropriate. Boxed legends, bold
math, no compression.

## Report write-up (populate `sec:lacl3` + a LaCl3 appendix), matching CaCl2
- **What you did**: the 3:1 quartic model + working assumptions, the
  **nonlinear-constraint (ParmEst) approach and why** (generalizes to multicomponent),
  the single dataset + its narrow-range/low-rejection caveats.
- **Results**: a results table (one row — state clearly there is a single processed
  LaCl3 dataset; note the raw-only MC4 run as future work), Method-A-vs-B agreement,
  fit quality.
- **What you found**: the adsorption-hypothesis verdict for La³⁺ (fit quality vs.
  NaCl/CaCl2; residual structure; effective-χ; whether |χ|→0), honestly caveated by
  the weak identifiability.
- **Show the results**: the per-dataset appendix block.
Keep working assumptions explicitly flagged; note that both the closed-form-free
constraint approach and the diagnostics carry directly into Step 14 (multicomponent).

## Guardrails
- New code under `analysis/lacl3/` (reuse `step2`/`cacl2` helpers + `plot_style`;
  a backward-compatible `load_sheet(data_file=, sheet_name=)` reuse is fine — verify
  NaCl/CaCl2 numbers unchanged after).
- Edit `main.tex` **only** within `sec:lacl3` and by adding a LaCl3 appendix
  subsection; do not touch §1/§2/§3, the model math, `refs.bib`, `data/`,
  `prior_analysis/`, or any NaCl/CaCl2 numbers.
- Compile clean (`latexmk -pdf`: 0 undefined, 0 overfull); rasterize + eyeball every
  figure; `latexmk -c` after.

## Report back
The LaCl3 `(|χ|, K)` ± CI, SSE, fit metric, multi-start %; Method-A-vs-B agreement;
the adsorption-diagnostic verdict (with the identifiability caveat); the `sec:lacl3`
+ appendix write-up added; and confirmation of a clean compile with NaCl/CaCl2
numbers unchanged.
