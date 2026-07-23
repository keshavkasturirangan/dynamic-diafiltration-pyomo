# Step 3 prompt for Sonnet — nonlinear confidence regions, multi-start, and convergence diagnostics

You are implementing **Step 3** of the DATA3 membrane-regression project. Work in
the repo `/Users/adowling/DowlingLab/Membranes/data3_regression`. Use the
documented conda environment: `conda activate data3-regression` (Python 3.11 with
numpy, scipy, pandas, openpyxl, matplotlib, pyomo 6.10, parmest, ipopt). Run
everything in that env and verify scripts execute cleanly.

## Goal
For the NaCl (1:1) Donnan regression, quantify parameter uncertainty with
**nonlinear confidence regions** and stress-test **optimizer reliability**. Do
this **first in scipy, then in Pyomo/ParmEst**, and compare. Two motivations from
the PI:
1. Bill observed the MATLAB curve-fitting optimizer **did not converge reliably** —
   we want to understand *why* and whether scipy / Ipopt do better.
2. We want nonlinear (not just linearized-ellipse) confidence regions, and to
   exercise ParmEst's confidence-region and multi-start features, which are not
   yet well tested by us.

## Background: the model and data (reuse Step 2)
The existing Step-2 code is `analysis/step2_nacl/regress_nacl.py` — **reuse its
data loading and model functions** (import or copy; do not duplicate logic
needlessly). Key facts:

- **Data:** `data/Rejection_Analysis.xlsx`, sheet
  `NF270_MC5 05.27.26_NaCl (2)`. Columns (1-indexed): B = `c_int` (NaCl
  interfacial conc, mM), H = `c_p` (permeate interfacial conc, mM), J = `J_w`
  (m³/m²/s). Use **rows 57–748** as the primary window (matches Bill; skips the
  low-concentration startup) and also report **full valid range rows 2–748**.
- **Transformed response:** `deltaC_m = (c_p*J_w) * l / D_NaCl_m` [mM], with
  `l = 80e-9` m and ambipolar `D_NaCl_m = 2*D_Na_m*D_Cl_m/(D_Na_m+D_Cl_m)`,
  `D_Na_m=1.33e-12`, `D_Cl_m=2.03e-12` m²/s.
- **Model:** `f(β; c_int, c_p) = sqrt(β1²+4β2 c_int²) − sqrt(β1²+4β2 c_p²)` with
  `β1 = |χ|` [mM], `β2 = δ*` [-]. The **physically correct** form restores a
  factor of ½: `f_corr = 0.5 * f`. Use **`f_corr` as the primary model** (it is
  the true Donnan co-ion concentration difference); also report the uncorrected
  `f` for comparison to Bill.
- **Step-2 point estimates (for verification):** uncorrected `f`, rows 57–748:
  `|χ| ≈ 25.18`, `δ* ≈ 0.0816`, SSE ≈ 56.52, n=692. Corrected `f_corr`:
  `|χ| ≈ 50.36`, `δ* ≈ 0.326`, SSE identical. (Note the *exact* homogeneity
  `f(2β1,4β2)=2f(β1,β2)`, so the two forms are related by |χ|→2×, δ*→4×.)

## Part A — scipy

Create `analysis/step3_nacl_uncertainty/scipy_uncertainty.py` implementing:

1. **Multi-start optimization (convergence diagnostic).** Draw a large set (e.g.
   200–500) of initial guesses over a physically plausible box — `|χ| ∈ [1, 200]`
   mM (uniform) and `δ* ∈ [1e-3, 10]` (log-uniform). From each start run
   `scipy.optimize.least_squares` on the residuals. Collect: converged `(|χ|, δ*)`,
   final SSE, success flag, iteration count. Report:
   - How many starts reach the global optimum vs. get stuck / fail.
   - Whether distinct local minima exist or it's one flat valley.
   - Compare optimizer **methods** (`lm`, `trf`, `dogbox`) and the effect of
     **parameter scaling**: |χ|~25 and δ*~0.08 differ by ~2–3 orders of
     magnitude, so `JᵀJ` is likely ill-conditioned. Test a rescaled/`x_scale`
     parameterization and/or **fitting `log10(δ*)`** instead of `δ*`, and report
     the condition number of `JᵀJ` before/after. Hypothesis to test: Bill's
     unreliable MATLAB convergence is a **conditioning/scaling** problem, not a
     multimodality problem.
   - Figure `docs/reports/figures/nacl_multistart.png`: scatter of converged
     starts in (|χ|, log10 δ*) space colored by final SSE, and/or a histogram of
     final SSE showing the basin structure.

2. **Nonlinear (likelihood-ratio / F-test) confidence region.** For Gaussian
   nonlinear least squares the exact joint region is
   `{β : SSE(β) ≤ SSE_min · (1 + (p/(n−p))·F_{p, n−p, 1−α})}` with p=2. Compute
   and **plot the 95% (and 50%, 90%) region** as SSE contours over the (|χ|, δ*)
   plane, overlaid on the `log10(SSE)` surface, with the optimum marked. Overlay
   the **linearized elliptical** region from the covariance
   `σ̂²(JᵀJ)⁻¹`. Comment on how different the nonlinear region is from the ellipse
   (the Step-2 contour showed a long curved valley → expect meaningful
   asymmetry / non-ellipticity). Figure:
   `docs/reports/figures/nacl_confidence_scipy.png`.

3. **Profile-likelihood intervals.** For each parameter, fix it on a grid and
   re-optimize the other; the profile-SSE crossing
   `SSE_min·(1 + F_{1,n−p,1−α}/(n−p))` gives an individual nonlinear CI. Report
   the profile CIs for |χ| and δ* and compare to the ± std-error (Wald) intervals.

4. *(Optional if time)* **Bootstrap** (resample residuals or cases, ~1000 reps) to
   get empirical parameter distributions; compare spread to the analytical region.

**Statistical caveat to state (not fix):** the ~692 points are a dense,
autocorrelated time series, not independent samples, so nominal confidence levels
are optimistic. Note this; do not attempt to correct it in Step 3.

## Part B — Pyomo / ParmEst

Create `analysis/step3_nacl_uncertainty/parmest_uncertainty.py` implementing the
same regression in Pyomo and using ParmEst. **Check the installed ParmEst API
before coding** (`python -c "import pyomo.contrib.parmest as p; print(p.__file__)"`
and read the current version's docstrings/examples) — Pyomo 6.10 uses the newer
**`Experiment`-class interface** (`parmest.Estimator([experiment_objects], ...)`),
not the old callback/dict API. Then:

1. Build a Pyomo model whose parameters are `|χ|`, `δ*` and whose objective is the
   SSE of `f_corr` vs. `deltaC_m` over the data; solve with **Ipopt**. Confirm the
   estimates match scipy (and Step 2).
2. **Confidence regions in ParmEst:** exercise and compare
   - `theta_est()` (point estimate + covariance),
   - the **likelihood-ratio / confidence-region test** (look for
     `likelihood_ratio_test` / `confidence_region_test`), and
   - **bootstrap** (`theta_est_bootstrap`).
   Plot ParmEst's confidence region over (|χ|, δ*) and overlay scipy's nonlinear
   region from Part A on the same axes:
   `docs/reports/figures/nacl_confidence_parmest.png`. Do they agree?
3. **Multi-start in ParmEst:** ParmEst recently added a multi-start capability —
   find it (check the `Estimator`/`theta_est` options and the parmest source for
   a `multistart`/`initialize_parmest`/scenario-sweep feature), use it, and
   compare its convergence reliability to the scipy multi-start from Part A. Does
   Ipopt (with/without scaling) converge more reliably than scipy `lm`?
4. Note any friction with the ParmEst API/features (this is partly an evaluation
   of those features for us).

## Deliverables
- `analysis/step3_nacl_uncertainty/scipy_uncertainty.py`
- `analysis/step3_nacl_uncertainty/parmest_uncertainty.py`
- `analysis/step3_nacl_uncertainty/README.md` — how to run + a results summary.
- Figures listed above under `docs/reports/figures/`.
- Keep all new files under `analysis/`; figures under `docs/reports/figures/`.
  **Do not modify** `prior_analysis/`, `docs/reports/main.tex`, or Step-2 code
  (importing from Step 2 is fine).

## Report back (as plain values/observations I will paste into the LaTeX writeup)
1. **Convergence diagnosis:** multi-start success rate (scipy) by method and with
   vs. without scaling / log-δ*; condition number of `JᵀJ`; conclusion on whether
   the unreliability is conditioning-driven and the recommended fix.
2. **Nonlinear confidence region (scipy):** the 95% region description and how it
   compares to the linearized ellipse; profile-likelihood CIs for |χ| and δ* vs.
   Wald intervals.
3. **ParmEst:** point estimates (match?), its confidence-region and bootstrap
   results, whether they agree with scipy, and how ParmEst multi-start convergence
   compares.
4. Any API friction or caveats worth telling the group.
5. Confirm every script runs cleanly in the `data3-regression` env and all figures
   were written.
