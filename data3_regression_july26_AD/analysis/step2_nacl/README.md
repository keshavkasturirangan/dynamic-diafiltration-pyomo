# Step 2: NaCl single-salt Donnan regression (Python port)

This directory ports Bill's preliminary MATLAB regression
(`prior_analysis/B_regression_concentration_NAClanal.m`) to Python/scipy and
reproduces its results, then reports a corrected version of the model with a
missing factor of 1/2 restored (see "Discrepancies to resolve" in
`docs/reports/main.tex`).

## How to run

```bash
/opt/anaconda3/bin/python3 analysis/step2_nacl/regress_nacl.py
```

Requires `numpy`, `scipy`, `pandas`, `openpyxl`, `matplotlib` (all present in
the anaconda3 environment referenced above). The script reads
`data/Rejection_Analysis.xlsx` (sheet `NF270_MC5 05.27.26_NaCl (2)`)
and writes two figures to `docs/reports/figures/`:

- `nacl_fit.png` — data vs. fitted model, `deltaC_m` vs. `c_int`.
- `nacl_residual_contour.png` — filled contour of `log10(SSE)` over
  `|chi| in linspace(20, 30, 110)`, `delta_star in logspace(-1.6, -1, 110)`,
  with the fitted optimum marked (mirrors the MATLAB `contourf`).

All numeric results are also printed to stdout.

## Model

Donnan equilibrium co-ion (Cl-) concentration in the membrane phase:

```
c_co^m(c_s) = ( -|chi| + sqrt(|chi|^2 + 4*delta_star*c_s^2) ) / 2
```

Bill's MATLAB fits (to a transformed response `deltaC_m`, built from the
measured solute flux):

```
f(beta; c_int, c_p) = sqrt(beta1^2 + 4*beta2*c_int^2) - sqrt(beta1^2 + 4*beta2*c_p^2)
```

which omits the factor of 1/2 from the Donnan root, i.e. it fits
`2*(c_co^m(c_int) - c_co^m(c_p))`. The "corrected" model restores the 1/2:

```
f_corrected(beta; c_int, c_p) = 0.5 * f(beta; c_int, c_p)
```

Both are fit here by nonlinear least squares
(`scipy.optimize.least_squares`, Levenberg-Marquardt) from the initial guess
`beta0 = (22.73, 0.09)`. Standard errors are computed from the linearized
covariance `cov = sigma^2 * (J^T J)^-1` with `sigma^2 = SSE / (n - p)`, using
the Jacobian at the optimum — the same approximation `curve_fit` uses
internally.

## Data

`data/Rejection_Analysis.xlsx`, sheet `NF270_MC5 05.27.26_NaCl (2)`.
Header on row 1; columns used: B = `c_int` (NaCl interfacial concentration,
mM), H = `c_p` (NaCl permeate interfacial concentration, smoothed, mM),
J = `J_w` (m^3/m^2/s).

- The MATLAB script reads a fixed row window **B57:B748**, which deliberately
  skips the low-concentration startup transient in rows 2-56.
- Valid (non-NaN) data actually spans rows **2-748** (n = 747); nothing valid
  exists beyond row 748, even though the workbook's reported sheet dimension
  goes to row 1590. No NaNs, zeros, or negative values were found for
  `c_int`, `c_p`, or `J_w` in rows 2-748.

## Results summary (see script stdout for full precision)

Physical constants: `D_Na_m = 1.33e-12 m^2/s`, `D_Cl_m = 2.03e-12 m^2/s`,
`l = 80e-9 m`, ambipolar `D_NaCl_m = 2*D_Na_m*D_Cl_m/(D_Na_m+D_Cl_m)
≈ 1.607e-12 m^2/s`.

| Fit | rows | n | \|chi\| (mM) | SE(chi) | delta* | SE(delta*) | SSE |
|---|---|---|---|---|---|---|---|
| Bill's MATLAB reproduction (uncorrected f) | 57-748 | 692 | 25.178 | 0.135 | 0.081585 | 0.000131 | 56.52 |
| Full valid range (uncorrected f) | 2-748 | 747 | 24.777 | 0.194 | 0.081267 | 0.000187 | 129.62 |
| Corrected model (f_corrected = 0.5*f) | 57-748 | 692 | 50.356 | 0.271 | 0.326342 | 0.000523 | 56.52 |
| Corrected model (f_corrected = 0.5*f) | 2-748 | 747 | 49.554 | 0.388 | 0.325069 | 0.000748 | 129.62 |

### Effect of the factor-of-2 correction

`f_matlab(beta1, beta2; c)` is **exactly homogeneous** under
`(beta1, beta2) -> (2*beta1, 4*beta2)`: pulling `4` out of each square root
gives `f_matlab(2*b1, 4*b2; c) = 2*f_matlab(b1, b2; c)` for all `c`. Since
`f_corrected = 0.5*f_matlab`, fitting `f_corrected` to the same data is
mathematically equivalent to fitting `f_matlab` to `2*deltaC_m`, which is
solved exactly at `(2*beta1_hat, 4*beta2_hat)` with an unchanged SSE surface.
This is confirmed numerically above: `|chi|` exactly doubles and
`delta_star` exactly quadruples between the uncorrected and corrected fits
(SSE identical to machine precision). So the factor-of-2 bug does **not**
leave `|chi|` "roughly unchanged" — it changes it by a clean, exact factor
of 2 — while `delta_star` rescales by exactly 4.

### Effect of using the full valid row range instead of MATLAB's window

Including the low-concentration startup transient (rows 2-56) shifts
`|chi|` down by about 1.6% (25.18 -> 24.78 mM) and `delta_star` down by about
0.4% (0.0816 -> 0.0813), with SSE increasing (more data, and the startup
points are not perfectly described by the same two-parameter model, e.g. see
the low-`c_int` end of `nacl_fit.png` where the fitted line runs slightly low
of the extrapolated data trend). The parameter estimates are not very
sensitive to this choice for this data set.

### Data anomalies

None found in the valid range (rows 2-748): no NaNs, zeros, or negative
values in `c_int`, `c_p`, or `J_w`. No valid data exists beyond row 748
despite the sheet's nominal dimension extending to row 1590.

### Units sanity check

`J_s = c_p * J_w` has units `[mM] * [m/s] = [mol/m^3] * [m/s] = [mol/m^2/s]`
(since mM = mmol/L = mol/m^3). Then
`deltaC_m = J_s * l / D_NaCl_m` has units
`[mol/m^2/s] * [m] / [m^2/s] = [mol/m^3] = [mM]`. Confirmed numerically:
example `deltaC_m` values are O(5-40) mM, consistent with the plotted range
in `nacl_fit.png`.
