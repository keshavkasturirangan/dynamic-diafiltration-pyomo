"""
Step 9: NaCl regression posed directly on the RAW measured time series
(permeate conductivity kappa_p(t)), with NO transformation of the response,
compared against the transformed-deltaC_m fit of report Sec. 2.1-2.2.

Motivation (errors-in-variables)
---------------------------------
The Sec. 2.1-2.2 fit forms deltaC_m = c_p*Jw*l/Dm and regresses THAT -- so the
noisy processed quantities c_p and Jw appear inside the response, and c_p
appears on both sides of the estimation (it is both a derived "predictor" via
the calibration and part of the thing being predicted). This mis-states the
uncertainty. This script instead places the error on the RAW MEASURED signal
(permeate conductivity) and predicts it from the model -- an
errors-in-variables-free (in the response) formulation.

Working forward model (exploratory; a reduced model, not the fully rigorous
DATA 2.0-style dynamic DAE on all raw channels -- see module docstring notes
below and docs/reports/main.tex sec:altform for the full writeup)
-----------------------------------------------------------------------------
Primary dataset: raw sheet NF270_MC5.xlsx / 05.27.26_NaCl (the raw twin of the
"(2)" sheet used in Sec 2.1-2.2), restricted to the SAME row window Bill's
MATLAB / Step 2-3-8 use (spreadsheet rows 57-748, n=692) for a direct,
apples-to-apples comparison -- verified via preprocess_nacl.py's raw-ts-row
<-> "(2)"-sheet-row correspondence (raw ts.iloc[i] <-> "(2)" sheet row i+2).

1. Known inputs per time t: retentate conductivity kappa_r(t) -> feed-bulk
   concentration c_f(t) via the linear calibration c=(kappa-a)/s (embedded
   MC5 values a=-63.706, s=76.685, the same "given" calibration
   preprocess_nacl.py uses for the "(2)" sheet's own raw source). Permeate
   mass m(t) -> water flux Jw(t) via preprocess_nacl.compute_jw (smoothed
   slope), unmodified. CP Phase 3 (docs/prompts/cp_phase3_ar1_rawmeas.md):
   feed interface c_int(t) is now the thin-film CP correction on c_f(t)
   (APPLY_CP=True, the default -- matching Phase 2's standing convention),
   with k imported from analysis/cp_authors_bulk.py (Phase 1's single source
   of truth). APPLY_CP=False reproduces the original c_int(t)=c_f(t) shortcut
   for comparison.
2. Model prediction of the measured kappa_p: the co-ion (Cl-) mass balance
   that defines deltaC_m in Sec. 2.1 (deltaC_m = J_s*l/Dm, J_s = c_p*Jw) is,
   before rearranging into "deltaC_m as the response", the implicit equation
   for c_p at fixed theta=(|chi|, delta*):
       Jw*c_p = (Dm/l) * [c_co_m(c_int; theta) - c_co_m(c_p; theta)]
   with the Donnan co-ion root (report Eq. 11co)
       c_co_m(c_s; theta) = 0.5*(-|chi| + sqrt(|chi|^2 + 4*delta*_c_s^2)).
   Non-CP: c_int=c_f is fixed, so this is a single equation in c_p, solved by
   scipy.optimize.brentq exactly as before (Descartes'-rule-style bracketing,
   g(0)<=0, g(c_int)>0, unique root in [0, c_int]).
   CP (the real subtlety -- see solve_cp_and_cint_CP): c_int is no longer
   fixed but itself a function of c_p via the CP equation
   c_int = c_p + (c_f-c_p)*exp(Jw/k), so (c_int, c_p) must be solved
   JOINTLY. Because the CP equation is affine in c_p, substituting
   c_int(c_p) into the flux balance above collapses the 2-equation system
   to a single equation in c_p alone -- an EXACT reduction, not an
   approximation or a seeded/decoupled shortcut -- solved by the same
   single brentq call (same computational cost as the non-CP path).
   Then kappa_p_hat(t) = a + s*c_p(theta,t) (invert the calibration).
3. Response = the RAW measured kappa_p(t) (permeate conductivity), no
   smoothing/transformation; residual = kappa_p_meas(t) - kappa_p_hat(t).

Measurement-error model (heteroscedastic, DATA 2.0-style; WORKING ASSUMPTION,
flagged -- no replicate conductivity measurements are available in this
dataset to calibrate sigma_kappa directly, so a generic constant-floor +
proportional-term form is assumed and stated explicitly rather than left
implicit):
    sigma_kappa(kappa) = SIGMA0_USCM + PROP_TERM * |kappa|
Weighted residual = (kappa_meas - kappa_hat(theta)) / sigma_kappa(kappa_meas);
theta fit by nonlinear least squares on these weighted residuals (bounded
|chi|,delta* >= 0, matching the CaCl2 fix for the same asymmetric-in-chi
issue: unlike the deltaC_m difference form, c_co_m alone is not chi-sign
invariant).

AR(1) correction: Method A only (GLS/whitening with actual n-p df, per the
Step 8b resolution of the double-counting issue -- reuses
analysis/step8_autocorrelation/ar1_correction.py's whiten/lag1_autocorr/
profile_ci_from_curve utilities, which are generic; a bespoke
Cochrane-Orcutt loop is used here since Step 8's is hardwired to the Sec 2.1
deltaC_m model).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/step9_raw_measurement/raw_regression.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq, least_squares
from scipy.stats import norm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step3_nacl_uncertainty"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step8_autocorrelation"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "preprocessing"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import scipy_uncertainty as s3  # noqa: E402
import ar1_correction as s8  # noqa: E402
import preprocess_nacl as pp  # noqa: E402
import cp_authors_bulk  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)
DATA_DIR = REPO_ROOT / "data"

ALPHA = 0.05
THETA0 = np.array([44.0, 0.3])
FIT_BOUNDS = ([0.0, 0.0], [np.inf, np.inf])

# CP Phase 3 (docs/prompts/cp_phase3_ar1_rawmeas.md): insert the same
# thin-film CP used in Phase 2 into the c_int(t) reconstruction here, instead
# of the previous non-CP c_int(t)=c_f(t) (bulk feed) shortcut. k imported
# from Phase 1's single source of truth (analysis/cp_authors_bulk.py), NOT
# recomputed.
APPLY_CP = True
K_NACL = cp_authors_bulk.K_SALT["NaCl"]

# Measurement-error model (working assumption, see module docstring)
SIGMA0_USCM = 5.0     # constant floor [uS/cm]
PROP_TERM = 0.01      # proportional term [-]

GIVEN_CALIB = pp.Calibration(slope=76.685, intercept=-63.706, n_points=0,
                              r2=float("nan"), source="given (embedded formula, matches (2) sheet)")


# ---------------------------------------------------------------------------
# Data loading (raw sheet, MATLAB row window for direct comparability)
# ---------------------------------------------------------------------------
def load_primary_raw_window():
    """Returns (time_s, Jw, c_f, kappa_p_meas, n_dropped). `c_f` is the
    retentate BULK feed concentration from the calibration -- CP Phase 3:
    this is no longer used directly as c_int (see predict_kappa_p / the
    coupled solve below); it is the CP mechanism's input bulk value, exactly
    as MC5(2)'s column B is in Phase 2."""
    raw = pp.load_raw_sheet(DATA_DIR / "NF270_MC5.xlsx", "05.27.26_NaCl")
    # window=51 is the validated/tuned value (analysis/preprocessing/README.md:
    # minimizes median relative error, ~3.8%, against the "(2)" sheet's own Jw
    # column) -- NOT the module's generic default of 15, which is markedly
    # noisier and was the root cause of an early (~1.6x) c_p bias in this
    # script's forward-model sanity check.
    Jw_full = pp.compute_jw(raw.ts["time_s"], raw.ts["mass_g"], window=51)
    c_f_full = pp.apply_calibration(raw.ts["ret_cond"], GIVEN_CALIB)
    kappa_p_full = raw.ts["perm_cond"].to_numpy()
    time_full = raw.ts["time_s"].to_numpy()

    # raw ts row i (0-indexed) <-> "(2)" sheet spreadsheet row i+2 (validated
    # in analysis/preprocessing/README.md); MATLAB window is spreadsheet rows
    # 57-748 inclusive -> raw ts iloc[55:747] (692 rows).
    lo = step2.MATLAB_FIRST_ROW - 2
    hi = step2.MATLAB_LAST_ROW - 2 + 1
    time_s = time_full[lo:hi]
    Jw = Jw_full[lo:hi]
    c_f = c_f_full[lo:hi]
    kappa_p = kappa_p_full[lo:hi]

    valid = np.isfinite(Jw) & np.isfinite(c_f) & np.isfinite(kappa_p) & (Jw > 0) & (c_f > 0)
    n_dropped = int((~valid).sum())
    return time_s[valid], Jw[valid], c_f[valid], kappa_p[valid], n_dropped


# ---------------------------------------------------------------------------
# Forward model: Donnan co-ion root (absolute, not the difference form) and
# the implicit c_p(theta, t) solve
# ---------------------------------------------------------------------------
def c_co_m(c_s, beta):
    """Report Eq. 11co: c_co^m(c^s;theta) = 0.5*(-|chi| + sqrt(chi^2+4*delta*(c^s)^2))."""
    chi, delta = beta
    return 0.5 * (-chi + np.sqrt(chi ** 2 + 4.0 * delta * c_s ** 2))


def solve_cp_noCP(beta, Jw_i, cint_i):
    """Non-CP path (c_int = c_f directly): unique root in [0, cint_i] of
    g(cp) = Jw*cp - (Dm/l)*(c_co_m(cint)-c_co_m(cp))."""
    def g(cp):
        return Jw_i * cp - (step2.D_NACL_M / step2.L_MEMBRANE) * (c_co_m(cint_i, beta) - c_co_m(cp, beta))

    lo, hi = 0.0, cint_i
    g_lo, g_hi = g(lo), g(hi)
    if g_lo > 0 or g_hi < 0:
        # Degenerate/boundary case (e.g. chi,delta ~ 0): clip to endpoint.
        cp = lo if abs(g_lo) < abs(g_hi) else hi
    else:
        cp = brentq(g, lo, hi, xtol=1e-10, rtol=1e-10)
    return cp, cint_i


def solve_cp_and_cint_CP(beta, Jw_i, c_f_i, k=K_NACL):
    """CP Phase 3 coupled solve. Two simultaneous equations in (c_int, c_p):
        (1) CP (thin-film, permeate unpolarized):
              c_int = c_p + (c_f - c_p) * exp(Jw/k)
        (2) Donnan flux balance (same as the non-CP path):
              Jw*c_p = (Dm/l) * (c_co_m(c_int; theta) - c_co_m(c_p; theta))
    Equation (1) is AFFINE in c_p for fixed (c_f, Jw, k), so c_int(c_p) can be
    substituted directly into (2), collapsing the 2x2 system to a single
    root-find in c_p alone -- an exact reduction, not an approximation.
    Bracket: c_int(c_p=0) = c_f*exp(Jw/k) (CP-elevated bulk) and
    c_int(c_p=c_f) = c_f (CP collapses to bulk when c_p=c_f, regardless of
    the exp factor), so g(0) = -(Dm/l)*c_co_m(c_int(0)) <= 0 and
    g(c_f) = Jw*c_f > 0 for Jw,c_f>0, chi,delta>=0 -- root guaranteed in
    [0, c_f] by the same Descartes'-rule-style argument as the non-CP path."""
    exp_factor = np.exp(Jw_i / k)

    def cint_of_cp(cp):
        return cp + (c_f_i - cp) * exp_factor

    def g(cp):
        cint = cint_of_cp(cp)
        return Jw_i * cp - (step2.D_NACL_M / step2.L_MEMBRANE) * (c_co_m(cint, beta) - c_co_m(cp, beta))

    lo, hi = 0.0, c_f_i
    g_lo, g_hi = g(lo), g(hi)
    if g_lo > 0 or g_hi < 0:
        cp = lo if abs(g_lo) < abs(g_hi) else hi
    else:
        cp = brentq(g, lo, hi, xtol=1e-10, rtol=1e-10)
    return cp, cint_of_cp(cp)


def predict_kappa_p(beta, Jw, c_f, calib=GIVEN_CALIB, apply_cp=None):
    """Returns (kappa_p_hat, c_p, c_int). Non-CP: c_int=c_f, single brentq
    per point (unchanged from the pre-Phase-3 script). CP: the coupled
    (c_int, c_p) solve above, also a single brentq per point (same cost)."""
    if apply_cp is None:
        apply_cp = APPLY_CP
    solver = solve_cp_and_cint_CP if apply_cp else solve_cp_noCP
    results = [solver(beta, jw_i, cf_i) for jw_i, cf_i in zip(Jw, c_f)]
    cp = np.array([r[0] for r in results])
    cint = np.array([r[1] for r in results])
    return calib.intercept + calib.slope * cp, cp, cint


# ---------------------------------------------------------------------------
# Measurement-error weighting and residuals
# ---------------------------------------------------------------------------
def sigma_kappa(kappa_meas):
    return SIGMA0_USCM + PROP_TERM * np.abs(kappa_meas)


def make_weighted_resid_fn(Jw, c_f, kappa_p_meas):
    w = 1.0 / sigma_kappa(kappa_p_meas)

    def resid(beta):
        kappa_hat, _, _ = predict_kappa_p(beta, Jw, c_f)
        return w * (kappa_p_meas - kappa_hat)

    return resid


# ---------------------------------------------------------------------------
# Fitting / uncertainty (generic analogs of step2.fit_model / s3's
# F-test/profile, parameterized by an arbitrary residual function)
# ---------------------------------------------------------------------------
def fit_weighted(resid_fn, beta0, n, label="raw kappa_p (WLS)"):
    res = least_squares(resid_fn, beta0, method="trf", bounds=FIT_BOUNDS)
    beta_hat = res.x
    sse = float(np.sum(res.fun ** 2))
    p = len(beta_hat)
    dof = max(n - p, 1)
    J = res.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma_scale2 = sse / dof
    cov = sigma_scale2 * JTJ_inv
    se = np.sqrt(np.diag(cov))
    return step2.FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


def profile_likelihood_generic(resid_fn, beta_hat, se, sse_min, n, p, param_idx,
                                alpha=0.05, n_grid=41, span_se=8.0):
    lo_bound = max(beta_hat[param_idx] - span_se * se[param_idx], 1e-6)
    grid = np.linspace(lo_bound, beta_hat[param_idx] + span_se * se[param_idx], n_grid)
    other_idx = 1 - param_idx
    profile_sse = np.zeros(len(grid))
    for i, val in enumerate(grid):
        def resid_fixed(other_val):
            beta = np.zeros(2)
            beta[param_idx] = val
            beta[other_idx] = other_val[0]
            return resid_fn(beta)

        x0 = np.array([max(beta_hat[other_idx], 1e-7)])
        res = least_squares(resid_fixed, x0, method="trf", bounds=([1e-8], [1e5]))
        profile_sse[i] = float(np.sum(res.fun ** 2))

    thresh = s3.f_test_threshold(sse_min, n, p, alpha)
    return grid, profile_sse, thresh


def cochrane_orcutt_generic(theta0, resid_fn, max_iter=50, tol=1e-10):
    """Method-A (GLS/whitening) Cochrane-Orcutt for an arbitrary time-ordered
    weighted residual function; reuses step8's whiten()/lag1_autocorr()
    (generic utilities) but a bespoke iteration loop, since step8's own
    cochrane_orcutt() is hardwired to the Sec. 2.1 deltaC_m model."""
    theta = np.asarray(theta0, dtype=float)
    r = resid_fn(theta)
    phi = s8.lag1_autocorr(r)

    for it in range(max_iter):
        def resid_whitened(th):
            return s8.whiten(resid_fn(th), phi)

        res = least_squares(resid_whitened, theta, method="trf", bounds=FIT_BOUNDS)
        theta_new = res.x
        r_new = resid_fn(theta_new)
        phi_new = s8.lag1_autocorr(r_new)
        d_theta = np.max(np.abs(theta_new - theta) / np.maximum(np.abs(theta), 1e-12))
        d_phi = abs(phi_new - phi)
        theta, phi = theta_new, phi_new
        if d_theta < tol and d_phi < tol:
            break

    return theta, phi, it + 1


def whitened_fit_result_generic(theta_hat, phi_hat, resid_fn, label="raw kappa_p GLS/AR(1)"):
    def resid(th):
        return s8.whiten(resid_fn(th), phi_hat)

    res = least_squares(resid, theta_hat, method="trf", bounds=FIT_BOUNDS)
    beta_hat = res.x
    sse = float(np.sum(res.fun ** 2))
    n = len(res.fun)
    p = len(beta_hat)
    dof = max(n - p, 1)
    J = res.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma_a2 = sse / dof
    cov = sigma_a2 * JTJ_inv
    se = np.sqrt(np.diag(cov))
    return step2.FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_raw_fit(time_s, c_int, kappa_p_meas, kappa_p_hat, cp_hat, outfile):
    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 3, figsize=(w, w * 0.32))

    order = np.argsort(c_int)
    ax = axes[0]
    ax.plot(c_int, kappa_p_meas, "o", color="tab:blue", markersize=3, alpha=0.4, label="measured")
    ax.plot(c_int[order], kappa_p_hat[order], "k--", lw=1.5, label="predicted")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\kappa_p}$ [$\mu$S/cm]")
    ax.set_title("(a) measured vs. predicted", fontsize=9)
    ax.legend(loc="best", fontsize=7)

    ax = axes[1]
    ax.plot(time_s, kappa_p_meas, "o", color="tab:blue", markersize=3, alpha=0.4, label="measured")
    ax.plot(time_s, kappa_p_hat, "-", color="k", lw=1.2, label="predicted")
    ax.set_xlabel(r"$\boldsymbol{t}$ [s]")
    ax.set_ylabel(r"$\boldsymbol{\kappa_p}$ [$\mu$S/cm]")
    ax.set_title(r"(b) $\kappa_p(t)$", fontsize=9)
    ax.legend(loc="best", fontsize=7)

    ax = axes[2]
    resid = kappa_p_meas - kappa_p_hat
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.plot(c_int, resid, "o", color="tab:red", markersize=3, alpha=0.5)
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"residual [$\mu$S/cm]")
    ax.set_title("(c) residual vs. $\\boldsymbol{c_{int}}$", fontsize=9)

    fig.suptitle("Raw-measurement fit: NaCl primary dataset (rows 57-748)", fontsize=10)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_compare(fr_transformed, wald_t, profile_t, fr_raw, wald_r, profile_r,
                  fr_raw_ar1, wald_r_ar1, profile_r_ar1, outfile):
    w = plot_style.fig_width("full") * 0.75
    fig, axes = plt.subplots(2, 1, figsize=(w, w * 1.1), sharex=True)

    labels = ["transformed\n(i.i.d.)", "raw meas.\n(i.i.d.)", "raw meas.\n(AR(1))"]
    x = np.arange(3)
    colors = [plot_style.OKABE_ITO["blue"], plot_style.OKABE_ITO["vermillion"],
              plot_style.OKABE_ITO["bluish_green"]]

    ax = axes[0]
    points = [fr_transformed.beta[0], fr_raw.beta[0], fr_raw_ar1.beta[0]]
    profiles = [profile_t[0], profile_r[0], profile_r_ar1[0]]
    for xi, pt, ci, c in zip(x, points, profiles, colors):
        lo, hi = max(pt - ci[0], 0), max(ci[1] - pt, 0)
        ax.errorbar(xi, pt, yerr=[[lo], [hi]], fmt="o", color=c, capsize=5,
                    markersize=8, elinewidth=2)
    ax.set_ylabel(r"$\boldsymbol{|\chi|}$ [mM]")

    ax = axes[1]
    points = [fr_transformed.beta[1], fr_raw.beta[1], fr_raw_ar1.beta[1]]
    profiles = [profile_t[1], profile_r[1], profile_r_ar1[1]]
    for xi, pt, ci, c in zip(x, points, profiles, colors):
        lo, hi = max(pt - ci[0], 0), max(ci[1] - pt, 0)
        ax.errorbar(xi, pt, yerr=[[lo], [hi]], fmt="o", color=c, capsize=5,
                    markersize=8, elinewidth=2)
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)

    for ax in axes:
        ax.label_outer()
    fig.suptitle("Transformed vs. raw-measurement point estimates\n"
                 "$\\pm$ 95% profile-likelihood CI", fontsize=10)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 100)
    print("Step 9: NaCl regression on the RAW measured time series (kappa_p)")
    print("=" * 100)

    print(f"APPLY_CP = {APPLY_CP}")
    time_s, Jw, c_f, kappa_p_meas, n_dropped = load_primary_raw_window()
    n = len(time_s)
    print(f"n = {n} (dropped {n_dropped} invalid rows from the 692-row MATLAB window)")

    resid_fn = make_weighted_resid_fn(Jw, c_f, kappa_p_meas)
    fr_raw = fit_weighted(resid_fn, THETA0, n)
    print(f"Raw-measurement WLS fit: |chi| = {fr_raw.beta[0]:.5f} mM, "
          f"delta* = {fr_raw.beta[1]:.6f}, weighted SSE = {fr_raw.sse:.4f}")

    z = norm.ppf(1 - ALPHA / 2)
    wald_raw = [(fr_raw.beta[k] - z * fr_raw.se[k], fr_raw.beta[k] + z * fr_raw.se[k])
                for k in range(2)]
    profile_raw = []
    for idx in range(2):
        grid, prof_sse, thresh = profile_likelihood_generic(
            resid_fn, fr_raw.beta, fr_raw.se, fr_raw.sse, n, fr_raw.p, idx, alpha=ALPHA,
        )
        profile_raw.append(s8.profile_ci_from_curve(grid, prof_sse, thresh))
    print(f"  i.i.d. Wald:    |chi| CI = ({wald_raw[0][0]:.4f}, {wald_raw[0][1]:.4f}), "
          f"delta* CI = ({wald_raw[1][0]:.5f}, {wald_raw[1][1]:.5f})")
    print(f"  i.i.d. Profile: |chi| CI = ({profile_raw[0][0]:.4f}, {profile_raw[0][1]:.4f}), "
          f"delta* CI = ({profile_raw[1][0]:.5f}, {profile_raw[1][1]:.5f})")

    # --- AR(1) correction (Method A, GLS/whitening) ---
    r0 = resid_fn(fr_raw.beta)
    phi_ols = s8.lag1_autocorr(r0)
    theta_ar1, phi_hat, n_iter = cochrane_orcutt_generic(fr_raw.beta, resid_fn)
    fr_raw_ar1 = whitened_fit_result_generic(theta_ar1, phi_hat, resid_fn)
    n_eff = n * (1 - phi_hat) / (1 + phi_hat)
    print(f"Lag-1 autocorrelation of weighted i.i.d. residuals: phi = {phi_ols:.4f}")
    print(f"Cochrane-Orcutt (Method A) converged in {n_iter} iterations: phi_hat = {phi_hat:.4f}")
    print(f"  n_eff = {n_eff:.2f} (n_eff/n = {n_eff/n:.4f})")
    print(f"AR(1)-corrected fit: |chi| = {fr_raw_ar1.beta[0]:.5f} mM, "
          f"delta* = {fr_raw_ar1.beta[1]:.6f}, SSE_whitened = {fr_raw_ar1.sse:.4f}")

    wald_raw_ar1 = [(fr_raw_ar1.beta[k] - z * fr_raw_ar1.se[k], fr_raw_ar1.beta[k] + z * fr_raw_ar1.se[k])
                    for k in range(2)]

    def whitened_resid_fn(th):
        return s8.whiten(resid_fn(th), phi_hat)

    profile_raw_ar1 = []
    for idx in range(2):
        grid, prof_sse, thresh = profile_likelihood_generic(
            whitened_resid_fn, fr_raw_ar1.beta, fr_raw_ar1.se, fr_raw_ar1.sse, n,
            fr_raw_ar1.p, idx, alpha=ALPHA,
        )
        profile_raw_ar1.append(s8.profile_ci_from_curve(grid, prof_sse, thresh))
    print(f"  AR(1) Wald:    |chi| CI = ({wald_raw_ar1[0][0]:.4f}, {wald_raw_ar1[0][1]:.4f}), "
          f"delta* CI = ({wald_raw_ar1[1][0]:.5f}, {wald_raw_ar1[1][1]:.5f})")
    print(f"  AR(1) Profile: |chi| CI = ({profile_raw_ar1[0][0]:.4f}, {profile_raw_ar1[0][1]:.4f}), "
          f"delta* CI = ({profile_raw_ar1[1][0]:.5f}, {profile_raw_ar1[1][1]:.5f})")

    # --- Transformed-fit comparison (Sec. 2.1-2.2 primary fit, unmodified) ---
    c_int_t, c_p_t, y_t = s3.load_primary_data()
    fr_t = step2.fit_model(step2.f_corrected, c_int_t, c_p_t, y_t, beta0=step2.BETA0,
                            label="transformed deltaC_m (Sec. 2.1-2.2)")
    wald_t = [(fr_t.beta[k] - z * fr_t.se[k], fr_t.beta[k] + z * fr_t.se[k]) for k in range(2)]
    profile_t = []
    for idx in range(2):
        grid, prof_sse, thresh = s3.profile_likelihood(
            fr_t.beta, fr_t.se, fr_t.sse, fr_t.n, fr_t.p, c_int_t, c_p_t, y_t, idx, alpha=ALPHA,
        )
        profile_t.append(s8.profile_ci_from_curve(grid, prof_sse, thresh))
    print()
    print("=" * 100)
    print("COMPARISON: transformed (Sec. 2.1-2.2) vs. raw-measurement (this script)")
    print("=" * 100)
    print(f"Transformed:        |chi| = {fr_t.beta[0]:.4f} mM (profile CI "
          f"{profile_t[0][0]:.4f}-{profile_t[0][1]:.4f}), "
          f"delta* = {fr_t.beta[1]:.5f} (profile CI {profile_t[1][0]:.5f}-{profile_t[1][1]:.5f})")
    print(f"Raw meas. (i.i.d.): |chi| = {fr_raw.beta[0]:.4f} mM (profile CI "
          f"{profile_raw[0][0]:.4f}-{profile_raw[0][1]:.4f}), "
          f"delta* = {fr_raw.beta[1]:.5f} (profile CI {profile_raw[1][0]:.5f}-{profile_raw[1][1]:.5f})")
    print(f"Raw meas. (AR(1)):  |chi| = {fr_raw_ar1.beta[0]:.4f} mM (profile CI "
          f"{profile_raw_ar1[0][0]:.4f}-{profile_raw_ar1[0][1]:.4f}), "
          f"delta* = {fr_raw_ar1.beta[1]:.5f} (profile CI {profile_raw_ar1[1][0]:.5f}-{profile_raw_ar1[1][1]:.5f})")
    d_chi_pct = 100 * (fr_raw.beta[0] - fr_t.beta[0]) / fr_t.beta[0]
    d_delta_pct = 100 * (fr_raw.beta[1] - fr_t.beta[1]) / fr_t.beta[1]
    print(f"Point-estimate shift transformed->raw: |chi| {d_chi_pct:+.2f}%, delta* {d_delta_pct:+.2f}%")
    print()

    # --- CP shift diagnostic (c_f -> c_int at the converged fit) ---
    if APPLY_CP:
        _, _, cint_check = predict_kappa_p(fr_raw.beta, Jw, c_f)
        shift_pct = 100 * (cint_check - c_f) / c_f
        print(f"CP shift in c_int (c_f -> c_int, at converged fit): "
              f"median {np.median(shift_pct):+.1f}%, max {np.max(shift_pct):+.1f}%")

    # --- Figures ---
    kappa_p_hat, cp_hat, cint_hat = predict_kappa_p(fr_raw.beta, Jw, c_f)
    plot_raw_fit(time_s, cint_hat, kappa_p_meas, kappa_p_hat, cp_hat,
                 FIG_DIR / "nacl_raw_fit.png")
    plot_compare(fr_t, wald_t, profile_t, fr_raw, wald_raw, profile_raw,
                 fr_raw_ar1, wald_raw_ar1, profile_raw_ar1,
                 FIG_DIR / "nacl_raw_compare.png")

    return dict(
        fr_t=fr_t, fr_raw=fr_raw, fr_raw_ar1=fr_raw_ar1, phi_hat=phi_hat, n_eff=n_eff,
        wald_t=wald_t, profile_t=profile_t, wald_raw=wald_raw, profile_raw=profile_raw,
        wald_raw_ar1=wald_raw_ar1, profile_raw_ar1=profile_raw_ar1,
    )


if __name__ == "__main__":
    main()
