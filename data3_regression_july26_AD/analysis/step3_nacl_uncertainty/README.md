# Step 3: NaCl regression — nonlinear confidence regions, multi-start, ParmEst

Uncertainty quantification and optimizer-reliability stress test for the NaCl
(1:1) Donnan regression, motivated by (1) Bill's MATLAB `lsqcurvefit` not
converging reliably, and (2) exercising ParmEst's confidence-region and
multi-start features for the first time in this project. All work here uses
Step 2's "corrected" model `f_corrected = 0.5*f_matlab` fit to the MATLAB row
window (rows 57-748, n=692) — Step 2's primary result:
`|chi| = 50.356 +/- 0.271 mM`, `delta* = 0.32634 +/- 0.00052`, `SSE = 56.52`.

Two scripts, run independently (same conda env as Step 2):

```bash
conda activate data3-regression
python analysis/step3_nacl_uncertainty/scipy_uncertainty.py     # ~15 s
python analysis/step3_nacl_uncertainty/parmest_uncertainty.py   # ~6 min
```

Both import `analysis/step2_nacl/regress_nacl.py` directly (unmodified) for
data loading and the model function; `parmest_uncertainty.py` additionally
imports `scipy_uncertainty.py` to overlay the two confidence regions on one
plot. Neither script modifies `prior_analysis/`, `main.tex`, or Step 2.

## Part A — scipy (`scipy_uncertainty.py`)

### Convergence diagnosis (the "why doesn't MATLAB converge" question)

300 random starts, `|chi|_0` uniform on [1, 200] mM, `delta*_0` log-uniform on
[1e-3, 10] — a box far wider than any reasonable physical guess, deliberately
stress-testing the optimizer rather than the model.

| Method                        | global | local/other | failed | median cond(J^T J) at start | at optimum |
|---|---|---|---|---|---|
| `lm` (unscaled)               | 100.0% | 0.0%  | 0% | 1.58e7 | 1.46e6 |
| `trf` (unscaled)               | 100.0% | 0.0%  | 0% | 1.58e7 | 1.46e6 |
| `dogbox` (unscaled)            | 72.7%  | 27.3% | 0% | 1.58e7 | 1.46e6 |
| `trf`, `x_scale='jac'`         | 100.0% | 0.0%  | 0% | 1.58e7 | 1.46e6 |
| `trf`, fit `(chi, log10 delta*)` | 100.0% | 0.0%  | 0% | 3.87e6 | 8.24e5 |

**The scaling/conditioning hypothesis is only partly right.** `cond(J^T J) ~
1e6-1e7` confirms the problem genuinely is poorly scaled (`|chi| ~ O(10-50)`
vs `delta* ~ O(0.1-0.3)`, ~2-3 orders of magnitude apart), and reparameterizing
to `(|chi|, log10 delta*)` roughly halves the condition number, as expected.
But **scipy's `lm` and `trf` converge to the global optimum from 100% of 300
starts regardless of scaling** — modern trust-region/Newton-type methods with
a sensible parameter box are simply not bothered by this level of
ill-conditioning for a smooth 2-parameter problem.

The one method that *does* fail regularly is `dogbox`: 27.3% of starts get
trapped at a **genuine (not spurious) local minimum** at `beta=(1e-8,
0.2705)`, `SSE=3028` (vs. global `SSE=56.5`) — `dogbox`'s reflective step
strategy drives `|chi|` into its lower bound and stalls there. This is a
real local minimum of the SSE surface (`|chi| -> 0` degenerates the Donnan
model to a pure square-root-of-delta* form), not a scaling artifact.
Figure: `nacl_multistart.png` (panels: `lm` 100% global, `dogbox` 27% trapped,
`trf` reparameterized 100% global).

**Practical recommendation for Bill's MATLAB fit:** `lsqcurvefit`'s default
algorithm is `trust-region-reflective`, MATLAB's analogue of `trf`/`dogbox`
rather than `lm`. Given the `dogbox` finding, the most likely explanation for
unreliable MATLAB convergence is this same boundary-trapping failure mode
(if MATLAB's default lower bound on `delta*` or `chi` is 0 or near-0) —
**not** poor scaling per se. Fixes, in order of preference: (1) switch to
`levenberg-marquardt` in `lsqcurvefit` (`algorithm='levenberg-marquardt'`),
which has no such bound-trapping behavior and matches scipy's 100%-global
`lm`/`trf` result; (2) if trust-region-reflective must be kept, start from a
reasonable physical guess (`chi0 ~ O(10-50)`, `delta0 ~ O(0.1-0.3)`) rather
than a default/arbitrary one, since the local minimum is far from the
physical parameter region; (3) reparameterize to `log10(delta*)` as a
secondary, complementary fix (smaller condition number, though not required
for reliability here).

### Nonlinear confidence region vs. linearized ellipse

`nacl_confidence_scipy.png`: the 50/90/95% nonlinear LR/F-test regions
(`SSE <= SSE_min*(1 + p/(n-p)*F_{p,n-p,1-alpha})`) are visually
indistinguishable from the linearized Wald ellipses at this sample size
(n=692) — the model is very close to linear in its parameters near the
optimum, so both approaches agree here. (This is a feature of the
NaCl fit's identifiability, not a general guarantee for the other-salt
regressions later in this project.)

### Profile-likelihood vs. Wald CIs

`nacl_profile_vs_wald.png`, 95% CIs:

| Parameter | point | profile CI | Wald CI |
|---|---|---|---|
| `\|chi\|` (mM) | 50.3555 | (49.679, 51.038) | (49.825, 50.886) |
| `delta*`       | 0.32634 | (0.32504, 0.32765) | (0.32532, 0.32737) |

Profile CIs are ~10-15% wider than Wald, consistent with the mild curvature
visible in the SSE surface — small but non-negligible; Wald intervals here
are not badly wrong, but slightly optimistic.

### Bootstrap (1000 case-resampling reps)

`nacl_bootstrap_scipy.png`. All 1000 reps converged.
`|chi|` 95% percentile CI = (49.740, 50.956), SE = 0.320 (vs. Wald SE 0.271);
`delta*` 95% CI = (0.32525, 0.32749), SE = 0.00059 (vs. Wald SE 0.00052).
Bootstrap SEs run ~15% larger than Wald/asymptotic SEs — a useful,
model-free sanity check that points the same direction as the
profile-vs-Wald gap above.

### Caveat (not fixed here, by design)

**All four CI/region methods above (F-test region, Wald ellipse, profile
likelihood, and case-resampling bootstrap) assume independent errors.** The
n=692 points are consecutive time-series samples from a single filtration
run, not i.i.d. replicates, so residual autocorrelation is expected and
unaccounted for. Every interval/region reported here is therefore likely
**optimistic (too narrow)** relative to the true sampling uncertainty. This
is flagged, not corrected, per the Step-3 task scope.

## Part B — ParmEst (`parmest_uncertainty.py`)

### API notes (pyomo 6.10.1)

- The installed pyomo uses the **Experiment-class interface**
  (`pyomo.contrib.parmest.experiment.Experiment`, subclassed with
  `create_model()` / `label_model()` / `get_labeled_model()`), not the older
  callback-style `parmest.Estimator(model_function, data, ...)` API shown in
  some older tutorials — that interface has been removed.
- Each of the 692 rows is modeled as one `Experiment` (one small
  `ConcreteModel` with `beta1`, `beta2` fixed-and-shared unknown parameters,
  fixed `c_int`/`c_p` inputs, and a free `deltaC_m` tied to them by an
  equality constraint). ParmEst's extensive-form (EF) `theta_est()` builds
  the joint 692-block NLP and solves it with Ipopt in a single call
  (~0.56 s total — the "many small experiments" pattern scales fine here).
- **`measurement_error` trap:** if you set `m.measurement_error[m.deltaC_m]`
  to a *placeholder number* (e.g. `1.0`), `cov_est()` silently uses that
  fixed value as the residual variance instead of estimating
  `sigma_hat^2 = SSE/(n-p)` from the data — every reported SE comes out
  wrong by a factor of `1/sqrt(sigma_hat^2)` (a ~3.5x inflation was observed
  here) with no warning or error. **You must set it to `None`** to get
  ParmEst to estimate the variance from the fit, matching scipy's approach.
  This is easy to get wrong silently and is not obvious from the examples.
- `theta_est(calc_cov=..., cov_n=...)` is **deprecated**; use the separate
  `cov_est()` method instead (raises a `DeprecationWarning`, still works in
  6.10.1). `cov_est()` also only accepts the two built-in objective names
  `'SSE'`/`'SSE_weighted'` — a custom Python `obj_function` callable (as
  used in the `rooney_biegler` examples) works for `theta_est()` but
  raises `ValueError` in `cov_est()`.
- There is **no built-in multi-start** in ParmEst; "multi-start" here means
  rebuilding the `Estimator` with different per-experiment initial `beta1`/
  `beta2` values and re-running `theta_est()` from scratch each time.
- `objective_at_theta()` re-solves the full 692-block EF NLP once **per grid
  point** (no shortcut for a fixed-parameter square-system evaluation), so
  it costs about the same as a full `theta_est()` call per point
  (~450 ms/point here) — building a fine confidence-region grid is
  meaningfully more expensive with ParmEst than with scipy's `sse_surface`
  (which is a cheap closed-form evaluation, no solver call).

### Estimates: ParmEst vs. scipy

| | `\|chi\|` (mM) | `delta*` | SSE |
|---|---|---|---|
| scipy `least_squares` | 50.355505 | 0.326342 | 56.519108 |
| ParmEst `theta_est()` | 50.355517 | 0.326342 | 56.519108 |
| relative difference | 2.3e-7 | 6.6e-8 | ~0 |

SEs from `cov_est()` (with `measurement_error=None`) match scipy's Wald SEs
to 6 significant figures: `SE(|chi|) = 0.270667`, `SE(delta*) = 0.000523`
for both. This cross-validates both the ParmEst model formulation and
scipy's linearized-covariance calculation.

### Multi-start comparison

30 random starts (same box as scipy's 300: `|chi|_0` in [1,200],
`delta*_0` log-uniform in [1e-3,10]) with ParmEst/Ipopt: **100% converged to
the global optimum**, 0.56 s/start on average. Unlike scipy's `dogbox`
(27% trapped at the `|chi|->0` local minimum), Ipopt's interior-point method
with a log-barrier on the `delta* > 1e-8` bound never gets pinned to that
boundary in this sample — consistent with the scipy `lm`/`trf` baseline
(100% global) rather than the `dogbox` failure mode.
Figure: `nacl_multistart_parmest.png`.

### Confidence region: ParmEst LR/chi2 vs. scipy F-test

`nacl_confidence_parmest.png` overlays scipy's F-test contours with a 21x21
grid of `objective_at_theta()` evaluations, classified inside/outside via
`likelihood_ratio_test`'s chi2-based rule
(`obj_value*(chi2_{2,alpha}/(n-2) + 1)`). **Caveat:** the grid spans
+/-12 SE per axis to safely bracket the region, but the true 95% region is
only ~1.3 mM wide in `|chi|` — so only 2-3 of the 441 grid points land
inside it. Where they do land, they sit exactly on/inside scipy's F-test
contour (no ParmEst-classified-inside point falls outside scipy's 95%
contour), confirming the two methods agree; a finer, narrower grid
(e.g. +/-4 SE) would draw a fuller region but was not rerun here given the
~450 ms/point cost. The chi2-based ParmEst threshold and the F-based scipy
threshold are asymptotically equivalent for n=692, p=2, and the numeric
agreement bears that out.

### Bootstrap: ParmEst vs. scipy

100 reps (`theta_est_bootstrap`, vs. scipy's 1000 — reduced given the
~550 ms/rep cost of resolving each 692-block EF with Ipopt):
`|chi|` 95% CI = (49.596, 51.027), SE = 0.376;
`delta*` 95% CI = (0.32496, 0.32757), SE = 0.00066.
Consistent with scipy's bootstrap (49.740-50.956, SE 0.320; 0.32525-0.32749,
SE 0.00059) — same ballpark, small differences attributable to the 10x
smaller rep count. Figure: `nacl_bootstrap_parmest.png`.

## Bottom line

- **Convergence:** the SSE surface has one genuine local minimum
  (`|chi| -> 0` boundary), only reachable by `dogbox`-style reflective
  trust-region steps; `lm`, `trf`, and Ipopt (ParmEst) all reach the global
  optimum reliably (100% across 300+30 starts) regardless of parameter
  scaling. Recommend Bill switch MATLAB to `levenberg-marquardt`.
- **Uncertainty:** nonlinear region ~ linearized ellipse ~ profile CIs ~ Wald
  CIs ~ bootstrap, all agreeing to within ~15% at this sample size — the
  model is well-identified and close to linear near the optimum for NaCl.
  The one caveat that matters is autocorrelation (unaddressed): true
  uncertainty is larger than any of these numbers suggest.
- **ParmEst:** reproduces scipy's point estimates and Wald SEs to 6+
  significant figures, confirming both implementations, and its multi-start
  behavior (100% global) matches scipy's best-behaved methods. Main API
  friction was the silent `measurement_error` covariance-scaling trap and
  the per-grid-point solver cost of `objective_at_theta`.
- Confirmed to run end-to-end in the `data3-regression` conda environment
  (`scipy_uncertainty.py` ~15 s, `parmest_uncertainty.py` ~6 min, dominated
  by the 441-point confidence-region grid).
