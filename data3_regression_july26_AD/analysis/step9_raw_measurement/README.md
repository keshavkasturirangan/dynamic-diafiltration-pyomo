# Step 9: regression on the raw measured time series (report §2.4)

Redoes the NaCl regression posed directly on the **raw measured permeate
conductivity** `kappa_p(t)`, with no transformation of the response, and
compares it against the transformed-`deltaC_m` fit of report §2.1-2.2 —
an errors-in-variables-free (in the response) alternative formulation.

## How to run

```bash
conda activate data3-regression
python analysis/step9_raw_measurement/raw_regression.py
```

Takes a few seconds. Reuses `analysis/step2_nacl/regress_nacl.py` (model
constants, `FitResult`, the primary transformed fit for comparison),
`analysis/step3_nacl_uncertainty/scipy_uncertainty.py` (`load_primary_data`,
`f_test_threshold`, i.i.d. `profile_likelihood`), `analysis/step8_autocorrelation/
ar1_correction.py` (`whiten`, `lag1_autocorr`, `profile_ci_from_curve` — all
generic utilities, reused unmodified), and `analysis/preprocessing/
preprocess_nacl.py` (raw sheet loading, `compute_jw`, calibration) — none of
these are modified.

## Why (the errors-in-variables problem)

The §2.1-2.2 fit forms `deltaC_m = c_p*Jw*l/Dm` and regresses that — so the
noisy *processed* quantities `c_p` and `Jw` appear inside the response, and
`c_p` appears on both sides of the estimation (both as a derived predictor,
via the calibration, and as part of what's predicted). This mis-states the
uncertainty. This script instead places the error on the **raw measured**
signal (permeate conductivity) and predicts it from the model.

## Dataset

Raw sheet `NF270_MC5.xlsx / 05.27.26_NaCl` (the raw twin of the "(2)" sheet
used in §2.1-2.2), restricted to the **same row window** Bill's MATLAB /
Steps 2-3-8 use (spreadsheet rows 57-748, `n=692`) for a direct,
apples-to-apples comparison. The raw-ts-row ↔ "(2)"-sheet-row correspondence
(`raw ts.iloc[i]` ↔ `"(2)"` sheet row `i+2`) was validated in
`analysis/preprocessing/README.md`.

**Bug caught during development**: the first version of `load_primary_raw_window()`
called `preprocess_nacl.compute_jw()` with its generic default
`window=15`, not the *validated/tuned* `window=51` (the value
`analysis/preprocessing/README.md` found minimizes median relative error,
~3.8%, against the "(2)" sheet's own `Jw` column). At `window=15`, the
model's predicted `c_p` was biased low by ~30-50% relative to the measured
`c_p` even at the known-good `theta` from §2.1-2.2 — a forward-model sanity
check (see git history) that caught this before it silently confounded the
comparison. Fixed to `window=51`.

## Working forward model (exploratory — a reduced model, not the fully
rigorous DATA 2.0-style dynamic DAE on all raw channels; see "Future work" below)

1. **Known inputs per time `t`.** Retentate conductivity `kappa_r(t)` →
   feed-bulk concentration `c_f(t)` via the linear calibration
   `c=(kappa-a)/s` (embedded MC5 values `a=-63.706`, `s=76.685`, the same
   "given" calibration `preprocess_nacl.py` uses for the "(2)" sheet's own
   raw source). Permeate mass `m(t)` → water flux `Jw(t)` via
   `preprocess_nacl.compute_jw` (window=51). Feed interface
   `c_int(t) = c_f(t)` directly (non-CP, matching the "(2)" sheet convention
   used throughout report §2).
2. **Model prediction of the measured `kappa_p`.** The co-ion (Cl⁻) mass
   balance that *defines* `deltaC_m` in §2.1 (`deltaC_m = J_s*l/Dm`,
   `J_s = c_p*Jw`) is, before rearranging into "`deltaC_m` as the response",
   the *implicit* equation for `c_p` at fixed `theta=(|chi|, delta*)`:
   ```
   Jw*c_p = (Dm/l) * [c_co_m(c_int; theta) - c_co_m(c_p; theta)]
   ```
   with the Donnan co-ion root (report Eq. 11co)
   `c_co_m(c_s; theta) = 0.5*(-|chi| + sqrt(|chi|^2 + 4*delta*c_s^2))`. This
   has a unique root in `[0, c_int]` (bracketed: `g(0) <= 0`, `g(c_int) > 0`
   for `|chi|,delta* >= 0` and `Jw,c_int > 0`), solved via
   `scipy.optimize.brentq`. Then `kappa_p_hat(t) = a + s*c_p(theta,t)`
   (invert the calibration).
3. **Response = the raw measured `kappa_p(t)`**, unsmoothed; residual =
   `kappa_p_meas(t) - kappa_p_hat(t)`.

## Measurement-error model (working assumption, flagged)

No replicate conductivity measurements exist in this dataset to calibrate
`sigma_kappa` directly, so a generic constant-floor + proportional-term form
is assumed (DATA 2.0-style) and stated explicitly:
```
sigma_kappa(kappa) = 5.0 [uS/cm] + 0.01 * |kappa|
```
Weighted residual = `(kappa_meas - kappa_hat(theta)) / sigma_kappa(kappa_meas)`;
`theta` fit by nonlinear least squares on these weighted residuals, bounded
`|chi|, delta* >= 0` (unlike the `deltaC_m` *difference* form, the standalone
`c_co_m` is not invariant to the sign of `chi`, so bounding matters here —
same class of fix as the CaCl2 script).

**Robustness check**: refitting with uniform (unweighted) residuals gives
`|chi|=140.3` mM, `delta*=1.114` — the same qualitative (large) shift from
the transformed fit persists regardless of the weighting scheme, confirming
this is a genuine effect of the raw-measurement formulation, not an artifact
of the ad hoc error-model weighting choice.

## AR(1) correction

Method A only (GLS/whitening with actual `n-p` degrees of freedom, per the
Step 8b resolution of the double-counting issue) — reuses Step 8's
`whiten`/`lag1_autocorr`/`profile_ci_from_curve` (generic utilities) with a
bespoke Cochrane-Orcutt loop for this model's residual function.

## Results

```
n = 692 (0 dropped from the 692-row MATLAB window)
Raw-measurement WLS fit (i.i.d.): |chi| = 114.894 mM, delta* = 1.01428, weighted SSE = 1262.19
  i.i.d. Wald:    |chi| CI = (112.788, 117.000), delta* CI = (1.00237, 1.02620)
  i.i.d. Profile: |chi| CI = (112.250, 117.604), delta* CI = (0.99940, 1.02951)
Lag-1 autocorrelation of weighted i.i.d. residuals: phi = 0.9293
Cochrane-Orcutt (Method A) converges in 25 iterations: phi_hat = 0.9989 (!)
  n_eff = 0.38  (n_eff/n = 0.0006)
AR(1)-corrected fit: |chi| = 4.921 mM, delta* = 0.00157, SSE_whitened = 31.19
  AR(1) Wald:    |chi| CI = (-8.638, 18.480), delta* CI = (-0.00087, 0.00402)
  AR(1) Profile: |chi| CI = (1.796, 60.265), delta* CI = (0.00055, 0.01155)
```

### Comparison to the transformed fit

| Formulation | `\|χ\|` (mM) | 95% profile CI | `δ*` | 95% profile CI |
|---|---|---|---|---|
| Transformed (§2.1-2.2) | 50.36 | (49.68, 51.04) | 0.3263 | (0.3250, 0.3277) |
| Raw measurement (i.i.d.) | 114.89 | (112.25, 117.60) | 1.0143 | (0.9994, 1.0295) |
| Raw measurement (AR(1)) | 4.92 | (1.80, 60.27) | 0.00157 | (0.00055, 0.01155) |

**Point-estimate shift, transformed → raw (i.i.d.): `|χ|` +128%, `δ*` +211%** —
a large, genuine effect, not a rounding-level correction.

## Verdict

**The errors-in-variables-free formulation moves the estimates dramatically**
— `|χ|` more than doubles, `δ*` more than triples. This is a materially
different conclusion from §2.1-2.2, not a minor refinement. Panel (a)/(b) of
`nacl_raw_fit.png` show the raw `kappa_p(t)` curve is fit very well in an
absolute sense (predicted tracks measured closely) — but panel (c)'s residual
plot shows clear systematic wave structure, and the extreme lag-1
autocorrelation (`phi=0.93`, worse than the transformed fit's `phi=0.77`)
confirms the absolute conductivity signal carries much more persistent
structure than the differenced/transformed response did.

**AR(1) is not a viable correction for this formulation**: Cochrane-Orcutt
drives `phi_hat` to the edge of the unit circle (`0.9989`, `n_eff` collapses
to `<1`) and the "corrected" point estimate and CI are not physically
meaningful (`δ*` CI includes negative values under Wald; the profile CI spans
nearly two orders of magnitude). This is itself an informative result: the
raw permeate-conductivity signal in a batch/semi-batch diafiltration run
rises smoothly and almost monotonically over the whole run, behaving more
like a slowly-varying trend (near-integrated/non-stationary) than a
stationary process fluctuating around a fixed nonlinear curve — AR(1)
whitening, which assumes stationarity, is not well-posed here. Likely
contributors (not disentangled by this exploratory script): (i) the reduced
forward model omits real physics (no CP correction, no time-varying
calibration/temperature effects) that would otherwise explain some of the
smooth residual trend; (ii) the measurement-error model (constant floor +
1% proportional term) is a working assumption not calibrated against
replicate data and could be misspecified; (iii) a genuinely richer noise
model (e.g. integrated/ARIMA, or a state-space formulation) may be needed
for raw-signal regression, consistent with the report's stated aspiration
that the fully rigorous version of this idea is a DATA 2.0-style dynamic DAE
on all raw channels, not the single-equation reduced form used here.

**Recommendation to the group**: the size of this discrepancy (not just its
existence) is the headline result — it suggests the transformed-`deltaC_m`
fit's `|χ|≈50`mM may not be robust to the errors-in-variables issue it
was designed to sidestep investigating. This is exploratory and should not
be treated as a corrected replacement for §2.1-2.2 without further work:
validating the measurement-error model, adding CP/temperature corrections,
and moving to a properly stationary or state-space noise treatment before
trusting either the point estimate or its uncertainty here.

## Figures (in `docs/reports/figures/`, `plot_style.py` conventions)

- `nacl_raw_fit.png` — full-width, 3 panels: (a) measured vs. predicted
  `kappa_p` vs. `c_int`, (b) `kappa_p(t)` time series, (c) residual vs.
  `c_int` (clear systematic wave structure).
- `nacl_raw_compare.png` — `(|χ|, δ*)` point estimates ± 95% profile CI,
  transformed (i.i.d.) vs. raw-measurement (i.i.d.) vs. raw-measurement
  (AR(1)) — the AR(1) panel's huge error bar visually communicates the
  method's breakdown here.

## Guardrails honored

- `data/`, `prior_analysis/`, `refs.bib` untouched.
- `analysis/step2_nacl/regress_nacl.py`, `analysis/step3_nacl_uncertainty/
  scipy_uncertainty.py`, `analysis/step8_autocorrelation/ar1_correction.py`,
  and `analysis/preprocessing/preprocess_nacl.py` all imported, not modified
  — the §2.1-2.2 transformed-fit point estimate/SSE printed above
  (`|χ|=50.3555`, `δ*=0.32634`, `SSE=56.519`) match the report exactly.
- `main.tex` edited only within `sec:altform`.
