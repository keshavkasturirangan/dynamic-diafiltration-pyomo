# Step 8: AR(1) autocorrelation-corrected confidence regions (report §2.3)

Corrects the NaCl confidence regions for the autocorrelation that the i.i.d.
error assumption (used throughout Steps 2-3) ignores. The ~692 points of the
primary NaCl dataset are a dense single-run time series, so residuals are
strongly autocorrelated and every i.i.d.-based interval (Wald, F-test,
profile, bootstrap) is too narrow.

## How to run

```bash
conda activate data3-regression
python analysis/step8_autocorrelation/ar1_correction.py
```

Takes a few seconds. Reuses `analysis/step2_nacl/regress_nacl.py` (data,
model, `FitResult`) and `analysis/step3_nacl_uncertainty/scipy_uncertainty.py`
(`load_primary_data`, `f_test_threshold`, `sse_surface`,
`profile_likelihood`) unmodified; the AR(1)/whitened analogs (parameterized
by `phi`, which s3's i.i.d. functions don't take) are new in this script.

## Dataset

Primary NaCl dataset, `NF270_MC5 05.27.26_NaCl (2)`, MATLAB row window
57-748 (`n=692`), corrected model `f_corrected`, matching report §2.1-2.2 and
Step 3 (`scipy_uncertainty.py`) exactly — same data, same OLS point estimate.

## Method (Prais-Winsten / Cochrane-Orcutt)

1. **AR(1) residual model.** `r_i = phi*r_{i-1} + a_i`, `a_i` iid
   `N(0, sigma_a^2)`. `phi_hat = sum(r_i r_{i-1}) / sum(r_i^2)` on the
   time-ordered OLS residuals.
2. **Effective sample size** `n_eff = n*(1-phi)/(1+phi)`; naive inflation
   factor `sqrt((1+phi)/(1-phi))`.
3. **GLS via Prais-Winsten whitening.** Because whitening is a *linear*
   operator on the time-ordered sequence, `whiten(y - f(theta)) ==
   whiten(y) - whiten(f(theta))` exactly, so the whitened residual is simply
   `whiten(r(theta))` with `tilde_r_1 = sqrt(1-phi^2)*r_1`,
   `tilde_r_i = r_i - phi*r_{i-1}` (`i>=2`). This means the GLS refit is a
   straightforward nonlinear least-squares fit of `whiten(r(theta; phi
   fixed))` — no separate whitening of `y` and `f` is needed, and scipy's
   `least_squares` computes the correct (whitened) Jacobian automatically
   via finite differences on this composed residual function. Cochrane-Orcutt
   iterates `phi_hat <-> theta_hat` to convergence (5 iterations here,
   tolerance `1e-10`).
4. **Corrected uncertainty — TWO CONSISTENT METHODS, not complementary.**
   Whitening and the `n_eff` deflation are *alternative* ways to account for
   the same autocorrelation; combining them (whitened SSE **and**
   `n_eff - p` df) double-counts the correction — this is exactly the bug in
   the first version of this script, corrected here (Step 8b). Reported
   both, side by side:
   - **Method A (GLS/whitening).** Wald covariance from the whitened
     Jacobian, `V = sigma_a^2 (Jtilde^T Jtilde)^-1`,
     `sigma_a^2 = SSE_whitened/(n-p)`; nonlinear F-test/profile region on the
     whitened SSE surface with the **actual** `n-p` degrees of freedom.
   - **Method B (`n_eff` heuristic).** No whitening: keep the OLS fit and
     OLS SSE surface, but replace `n` by `n_eff` everywhere — Wald
     `SE_B = SE_OLS * sqrt(n/n_eff)`; F-test/profile region on the OLS SSE
     surface with `n_eff - p` degrees of freedom.
5. **Sufficiency of AR(1).** Compare the ACF of the OLS residuals vs. the
   whitened (GLS) residuals against a `+/-1.96/sqrt(n)` white-noise band.

## Results

```
n = 692
OLS fit:  |chi| = 50.35551 mM,  delta* = 0.326342,  SSE = 56.5191
Lag-1 autocorrelation of OLS residuals: phi_hat = 0.7734
Cochrane-Orcutt converged in 5 iterations: phi_hat = 0.7758
GLS fit:  |chi| = 49.86961 mM,  delta* = 0.325167,  SSE_whitened = 21.4271
n_eff = 87.35  (n_eff/n = 0.1262, i.e. ~87 "independent" observations out of 692)
Naive CI-inflation factor sqrt((1+phi)/(1-phi)) = 2.815
Point-estimate shift OLS -> GLS: |chi| -0.965%, delta* -0.360%
```

**|χ| and δ\* barely move** (both shift by <1%) — the point estimate is
robust to the autocorrelation correction; only the *uncertainty* changes.

### i.i.d. vs. Method A vs. Method B 95% CIs

| Param | point | Wald (i.i.d.) | Wald (A) | Wald (B) | fold A / B |
|---|---|---|---|---|---|
| `\|χ\|` (mM) | 50.3555 | (49.825, 50.886) | (48.430, 51.309) | (48.862, 51.849) | **2.71× / 2.81×** |
| `δ*` | 0.3263 | (0.3253, 0.3274) | (0.3224, 0.3279) | (0.3235, 0.3292) | **2.70× / 2.81×** |

| Param | point | Profile (i.i.d.) | Profile (A) | Profile (B) | fold A / B |
|---|---|---|---|---|---|
| `\|χ\|` (mM) | 50.3555 | (49.679, 51.038) | (48.048, 51.734) | (48.417, 52.344) | **2.71× / 2.89×** |
| `δ*` | 0.3263 | (0.3250, 0.3277) | (0.3217, 0.3287) | (0.3226, 0.3302) | **2.70× / 2.89×** |

**Methods A and B agree closely** (2.7–2.9× for both Wald and profile), and
**profile tracks Wald under each method** — profile/Wald ratio for `|χ|` is
1.281 (i.i.d.), 1.280 (A), 1.315 (B): all close to the i.i.d.\ case's own
ratio, i.e. mild, comparable nonlinearity throughout. **This confirms the two
corrections are consistent alternatives, and that the region does NOT become
dramatically more nonlinear once the correction is applied without
double-counting.**

An earlier version of this script combined Method A's whitened SSE surface
with Method B's `n_eff - p` degrees of freedom, which double-counts the
correlation and produced a spurious ~7.9× profile fold-widening with a
visibly more elongated region — that number should not be used; see git
history for the intermediate (incorrect) version.

**Open question for the PI**: Methods A and B are both defensible and agree
closely here — which should be the report's headline convention? Separately,
is the `n_eff` heuristic (Method B) even well justified given AR(1) is not
fully sufficient (below) — its derivation assumes the AR(1) model is exact,
whereas Method A's whitened residuals still show some structure. Method A (a
full GLS refit) is arguably more rigorous on that basis, but both are
reported here rather than silently picking one.

### ACF sufficiency check

```
+/-1.96/sqrt(n) band = +/-0.0745
OLS residual ACF:      lag-1 = 0.7734;  30/30 lags exceed the band
Whitened residual ACF: lag-1 = -0.2971; 12/30 lags exceed the band
```

AR(1) removes most but **not all** of the structure: the whitened ACF is
far smaller in magnitude than the raw ACF and no longer shows the raw
series' slow monotonic decay, but a substantial fraction of lags (12/30,
including a non-trivial negative lag-1 of -0.30) still exceed the white-noise
band. **Verdict: AR(1) is a large first-order improvement but is not fully
sufficient** — the negative lag-1 in the whitened series suggests
overcorrection at short lags (or that a moving-average/short-memory
component remains), and higher-order AR(p) or ARMA models are warranted as
future work (as already flagged in the report's roadmap).

## Figures (in `docs/reports/figures/`, `plot_style.py` conventions)

- `nacl_ar1_acf.png` — half-width; OLS vs. whitened residual ACF with the
  `±1.96/√n` significance band.
- `nacl_ar1_regions.png` — full-width, two panels (zoom + full view) showing
  the i.i.d. vs. Method-A vs. Method-B 95% F-test regions overlaid on the
  `(|χ|, δ*)` plane. All three are visible at comparable scale in both
  panels since Methods A and B give similar (~2.7–2.9×) widening once
  applied consistently (no zoom/full disparity like the earlier
  double-counted ~8× version).

## Guardrails honored

- `data/`, `prior_analysis/`, `refs.bib` untouched.
- `analysis/step2_nacl/regress_nacl.py` and
  `analysis/step3_nacl_uncertainty/scipy_uncertainty.py` untouched (imported,
  not modified) — the OLS point estimate/SSE printed above match Step 3's
  reported values exactly (`|χ|=50.3555`, `δ*=0.326342`, `SSE=56.519`).
- `main.tex` edited only by adding `sec:naclAR1` and updating the one
  "Planned correction" pointer note in `sec:naclUQ`.
