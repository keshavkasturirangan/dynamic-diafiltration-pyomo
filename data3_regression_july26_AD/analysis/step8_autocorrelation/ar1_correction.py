"""
Step 8: autocorrelation-corrected (AR(1)) confidence regions for the primary
NaCl regression.

Motivation
----------
The Step-2/3 regressions minimize SSE(theta) = sum_i r_i^2, r_i = y_i -
f(theta; x_i), assuming the r_i are i.i.d. N(0, sigma^2). The ~692 points are
a dense single-run TIME SERIES, so the residuals are autocorrelated and every
i.i.d.-based interval (Wald, F-test, profile, bootstrap) is too narrow.

IMPORTANT correction (Step 8b): the first version of this script combined
GLS whitening with an n_eff-based F-test threshold, which DOUBLE-COUNTS the
autocorrelation correction (whitening already encodes the information loss
in the whitened Jacobian; n_eff is a separate, alternative heuristic for the
same effect, not a second correction to stack on top). This version computes
the correction two CONSISTENT, non-double-counting ways and reports both
(see docs/reports/main.tex, sec:naclAR1, for the full write-up):

- Method A (GLS/whitening, rigorous): Prais-Winsten + Cochrane-Orcutt, with
  the F-test/profile region on the WHITENED SSE surface using actual n-p
  degrees of freedom (NOT n_eff - p).
- Method B (effective-sample-size heuristic): no whitening; keep the OLS fit
  and OLS SSE surface, but replace n by n_eff everywhere uncertainty is
  computed (Wald SE scaled by sqrt(n/n_eff); F-test/profile region on the
  OLS SSE surface using n_eff - p degrees of freedom).

Methods A and B are ALTERNATIVES, not complementary -- each is a complete,
internally-consistent correction on its own; this script reports both so the
choice (and the double-counting risk of combining them) is visible on paper
rather than picked silently.

Method (Prais-Winsten / Cochrane-Orcutt, shared machinery)
------------------------------------------------------------
1. AR(1) residual model r_i = phi*r_{i-1} + a_i, a_i ~ iid N(0, sigma_a^2).
   Estimate phi_hat = sum(r_i r_{i-1}) / sum(r_i^2) (lag-1 sample
   autocorrelation of the OLS residuals, time-ordered).
2. Effective sample size n_eff = n*(1-phi)/(1+phi); inflation factor
   sqrt(n/n_eff) = sqrt((1+phi)/(1-phi)) -- this is exactly Method B's Wald
   SE inflation.
3. GLS via Prais-Winsten whitening (Method A only): since whitening is a
   LINEAR operator on the time-ordered residual sequence,
   whiten(y - f(theta)) == whiten(y) - whiten(f(theta)) exactly, so the
   whitened residual is simply whiten(r(theta)) with
       tilde_r_1 = sqrt(1-phi^2) * r_1
       tilde_r_i = r_i - phi*r_{i-1},  i >= 2
   Cochrane-Orcutt iterates: fit theta (OLS) -> estimate phi from r(theta)
   -> refit theta by NLS on whiten(r(theta; phi fixed)) -> re-estimate phi
   from the new r(theta) -> repeat to convergence.
4. Corrected uncertainty:
   - Method A: Wald covariance from the whitened Jacobian,
     V = sigma_a^2 (Jtilde^T Jtilde)^-1, sigma_a^2 = SSE_whitened/(n-p); the
     nonlinear F-test/profile region on the whitened SSE surface with n-p df.
   - Method B: Wald SE_B = SE_OLS * sqrt(n/n_eff); the nonlinear
     F-test/profile region on the (unmodified) OLS SSE surface with
     n_eff - p df.
5. Sufficiency of AR(1): compare the ACF of the OLS residuals vs. the
   whitened (GLS) residuals; if the whitened ACF is within the
   +/-1.96/sqrt(n) band beyond lag 0, AR(1) suffices.

Reuses analysis/step2_nacl/regress_nacl.py (data, model, OLS fit) and
analysis/step3_nacl_uncertainty/scipy_uncertainty.py (F-test threshold,
i.i.d. SSE surface, i.i.d. profile likelihood) unmodified -- Method B reuses
s3's i.i.d. functions directly (passing n_eff in place of n is the entire
method), Method A's whitened analogs are new (parameterized by phi, which
s3's functions don't take).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/step8_autocorrelation/ar1_correction.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import norm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step3_nacl_uncertainty"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import scipy_uncertainty as s3  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
MAX_CO_ITER = 50
CO_TOL = 1e-10
ACF_MAX_LAG = 30

# CP Phase 3 (docs/prompts/cp_phase3_ar1_rawmeas.md): this script reuses
# s3.load_primary_data() -> step2.load_sheet(), so it already inherits
# step2's APPLY_CP default (True since Phase 2). This module-level flag is
# threaded through explicitly (rather than relying on the imported default)
# so the CP state is visible here and independently toggleable for the
# regression test.
APPLY_CP = True


# ---------------------------------------------------------------------------
# AR(1) machinery
# ---------------------------------------------------------------------------
def lag1_autocorr(r: np.ndarray) -> float:
    """phi_hat = sum_{i=2}^n r_i r_{i-1} / sum_{i=1}^n r_i^2, time-ordered r."""
    return float(np.sum(r[1:] * r[:-1]) / np.sum(r ** 2))


def whiten(v: np.ndarray, phi: float) -> np.ndarray:
    """Prais-Winsten whitening transform of a time-ordered sequence."""
    out = np.empty_like(v)
    out[0] = np.sqrt(1.0 - phi ** 2) * v[0]
    out[1:] = v[1:] - phi * v[:-1]
    return out


def acf(r: np.ndarray, max_lag: int = ACF_MAX_LAG) -> np.ndarray:
    """Sample autocorrelation function of a (mean-zero) residual series,
    lags 0..max_lag, normalized by sum(r^2) (lag-0 value is exactly 1)."""
    n = len(r)
    denom = np.sum(r ** 2)
    return np.array([np.sum(r[k:] * r[:n - k]) / denom for k in range(max_lag + 1)])


def whitened_resid_fn(theta, phi, c_int, c_p, y):
    r = y - step2.f_corrected(theta, c_int, c_p)
    return whiten(r, phi)


def cochrane_orcutt(theta0, c_int, c_p, y, max_iter=MAX_CO_ITER, tol=CO_TOL):
    """Iterate phi_hat <-> theta_hat to convergence (Cochrane-Orcutt).
    Returns (theta_hat, phi_hat, n_iterations, history)."""
    theta = np.asarray(theta0, dtype=float)
    r = y - step2.f_corrected(theta, c_int, c_p)
    phi = lag1_autocorr(r)
    history = [(theta.copy(), phi)]

    for it in range(max_iter):
        def resid(th):
            return whitened_resid_fn(th, phi, c_int, c_p, y)

        res = least_squares(resid, theta, method="lm")
        theta_new = res.x
        r_new = y - step2.f_corrected(theta_new, c_int, c_p)
        phi_new = lag1_autocorr(r_new)

        history.append((theta_new.copy(), phi_new))
        d_theta = np.max(np.abs(theta_new - theta) / np.maximum(np.abs(theta), 1e-12))
        d_phi = abs(phi_new - phi)
        theta, phi = theta_new, phi_new
        if d_theta < tol and d_phi < tol:
            break

    return theta, phi, it + 1, history


def whitened_fit_result(theta_hat, phi_hat, c_int, c_p, y, label="GLS/AR(1)"):
    """FitResult-like object (reuses step2.FitResult) for the converged GLS
    fit: Jacobian and covariance computed on the WHITENED residual function
    at (theta_hat, phi_hat) -- this is the correct GLS covariance, distinct
    from an OLS covariance evaluated on whitened data."""
    def resid(th):
        return whitened_resid_fn(th, phi_hat, c_int, c_p, y)

    res = least_squares(resid, theta_hat, method="lm")
    beta_hat = res.x
    r_tilde = res.fun
    sse = float(np.sum(r_tilde ** 2))
    n, p = len(y), len(beta_hat)
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
# Whitened analogs of step3's F-test / profile-likelihood machinery
# (parameterized by phi, unlike s3's, which assume i.i.d. errors)
# ---------------------------------------------------------------------------
def sse_surface_whitened(beta1_grid, beta2_grid, phi, c_int, c_p, y):
    SSE = np.zeros((len(beta1_grid), len(beta2_grid)))
    for i, b1 in enumerate(beta1_grid):
        for j, b2 in enumerate(beta2_grid):
            r_tilde = whitened_resid_fn((b1, b2), phi, c_int, c_p, y)
            SSE[i, j] = np.sum(r_tilde ** 2)
    return SSE


def profile_likelihood_whitened(beta_hat, se, sse_min, dof_n, p, phi, c_int, c_p, y,
                                 param_idx, alpha=0.05, n_grid=41, span_se=8.0):
    """Profile likelihood on the WHITENED SSE surface (Method A). `dof_n` is
    the value used as "n" in the F-test threshold's denominator df (n-p) --
    pass the ACTUAL n for Method A (the correct, non-double-counting choice;
    see module docstring)."""
    lo_bound = max(beta_hat[param_idx] - span_se * se[param_idx], 1e-6)
    grid = np.linspace(lo_bound, beta_hat[param_idx] + span_se * se[param_idx], n_grid)
    other_idx = 1 - param_idx
    profile_sse = np.zeros(len(grid))
    for i, val in enumerate(grid):
        def resid_fixed(other_val):
            beta = np.zeros(2)
            beta[param_idx] = val
            beta[other_idx] = other_val[0]
            return whitened_resid_fn(beta, phi, c_int, c_p, y)

        x0 = np.array([max(beta_hat[other_idx], 1e-7)])
        res = least_squares(resid_fixed, x0, method="trf", bounds=([1e-8], [1e5]))
        profile_sse[i] = float(np.sum(res.fun ** 2))

    thresh = s3.f_test_threshold(sse_min, dof_n, p, alpha)
    return grid, profile_sse, thresh


def profile_ci_from_curve(grid, profile_sse, thresh):
    imin = int(np.argmin(profile_sse))
    lo = grid[0]
    for i in range(imin, 0, -1):
        if profile_sse[i - 1] > thresh >= profile_sse[i]:
            lo = np.interp(thresh, [profile_sse[i], profile_sse[i - 1]],
                            [grid[i], grid[i - 1]])
            break
    else:
        lo = grid[0]
    hi = grid[-1]
    for i in range(imin, len(grid) - 1):
        if profile_sse[i] <= thresh < profile_sse[i + 1]:
            hi = np.interp(thresh, [profile_sse[i], profile_sse[i + 1]],
                            [grid[i], grid[i + 1]])
            break
    else:
        hi = grid[-1]
    return lo, hi


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_acf(acf_raw, acf_white, n, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 0.85))
    lags = np.arange(len(acf_raw))
    band = 1.96 / np.sqrt(n)
    ax.axhspan(-band, band, color="gray", alpha=0.2, label=r"$\pm 1.96/\sqrt{n}$ band")
    ax.stem(lags - 0.15, acf_raw, linefmt="C0-", markerfmt="C0o", basefmt=" ",
            label="OLS (i.i.d.) residuals")
    ax.stem(lags + 0.15, acf_white, linefmt="C1-", markerfmt="C1s", basefmt=" ",
            label="whitened (GLS/AR(1)) residuals")
    ax.axhline(0, color="k", lw=1)
    ax.set_xlabel("lag")
    ax.set_ylabel("sample ACF")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.25), ncol=1, fontsize=8,
              borderaxespad=0)
    fig.suptitle("Residual ACF: raw vs. whitened\n(NaCl primary dataset, rows 57-748)",
                  fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def _draw_three_regions(ax, fr_ols, fr_gls, phi_hat, n, n_eff, b1_grid, b2_grid, c_int, c_p, y):
    """Draw the i.i.d., Method-A (GLS/whitening, n-p df), and Method-B
    (n_eff heuristic on the OLS surface) 95% F-test regions on the same
    axes -- all three consistent (non-double-counting) corrections."""
    SSE_ols = s3.sse_surface(b1_grid, b2_grid, c_int, c_p, y)

    thresh_iid = s3.f_test_threshold(fr_ols.sse, fr_ols.n, fr_ols.p, ALPHA)
    ax.contour(b1_grid, b2_grid, SSE_ols.T, levels=[thresh_iid], colors=["tab:blue"],
               linewidths=2.5)
    ax.plot(fr_ols.beta[0], fr_ols.beta[1], marker="*", color="tab:blue", markersize=12,
            markeredgecolor="k", markeredgewidth=0.5, linestyle="none")

    # Method B: OLS SSE surface (already computed above), n_eff - p df.
    thresh_B = s3.f_test_threshold(fr_ols.sse, n_eff, fr_ols.p, ALPHA)
    ax.contour(b1_grid, b2_grid, SSE_ols.T, levels=[thresh_B], colors=["tab:green"],
               linewidths=2.5, linestyles="dotted")

    # Method A: whitened SSE surface, actual n - p df.
    SSE_gls = sse_surface_whitened(b1_grid, b2_grid, phi_hat, c_int, c_p, y)
    thresh_A = s3.f_test_threshold(fr_gls.sse, n, fr_gls.p, ALPHA)
    ax.contour(b1_grid, b2_grid, SSE_gls.T, levels=[thresh_A], colors=["tab:red"],
               linewidths=2.5, linestyles="dashed")
    ax.plot(fr_gls.beta[0], fr_gls.beta[1], marker="*", color="tab:red", markersize=12,
            markeredgecolor="k", markeredgewidth=0.5, linestyle="none")


def plot_regions(fr_ols, fr_gls, phi_hat, n, n_eff, c_int, c_p, y, outfile):
    """Two panels: LEFT a zoom on the i.i.d. region, RIGHT a wider view --
    with the corrected framing, Methods A and B give SIMILAR (~2.7-2.8x)
    widening, so (unlike the double-counted first draft) all three regions
    are visible without an extreme zoom/full split; the zoom panel is kept
    only because the i.i.d. region is still the tightest of the three."""
    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 2, figsize=(w, w * 0.52))

    se_o = fr_ols.se
    b1_zoom = np.linspace(fr_ols.beta[0] - 6 * se_o[0], fr_ols.beta[0] + 6 * se_o[0], 200)
    b2_zoom = np.linspace(fr_ols.beta[1] - 6 * se_o[1], fr_ols.beta[1] + 6 * se_o[1], 200)
    b2_zoom = b2_zoom[b2_zoom > 0]
    _draw_three_regions(axes[0], fr_ols, fr_gls, phi_hat, n, n_eff, b1_zoom, b2_zoom, c_int, c_p, y)
    axes[0].set_title("zoom: i.i.d. region", fontsize=9)

    # Full-view grid sized by the larger (AR(1)) uncertainty, centered to
    # cover both optima, so all three regions render as closed contours.
    se_g = fr_gls.se
    center1 = 0.5 * (fr_ols.beta[0] + fr_gls.beta[0])
    center2 = 0.5 * (fr_ols.beta[1] + fr_gls.beta[1])
    b1_full = np.linspace(center1 - 5 * se_g[0], center1 + 5 * se_g[0], 250)
    b2_full = np.linspace(center2 - 5 * se_g[1], center2 + 5 * se_g[1], 250)
    b2_full = b2_full[b2_full > 0]
    _draw_three_regions(axes[1], fr_ols, fr_gls, phi_hat, n, n_eff, b1_full, b2_full, c_int, c_p, y)
    axes[1].set_title("full view", fontsize=9)

    for ax in axes:
        ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    axes[0].set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")

    shared_handles = [
        plt.Line2D([], [], color="tab:blue", lw=2.5, label="95% i.i.d. F-test region"),
        plt.Line2D([], [], color="tab:red", lw=2.5, ls="--",
                   label="95% Method A (GLS/whitening, $n-p$ df)"),
        plt.Line2D([], [], color="tab:green", lw=2.5, ls=":",
                   label=r"95% Method B ($n_{\mathrm{eff}}$ heuristic)"),
    ]
    fig.legend(handles=shared_handles, loc="lower center", bbox_to_anchor=(0.5, -0.16),
               ncol=1, fontsize=8, borderaxespad=0)
    fig.suptitle("i.i.d. vs. AR(1)-corrected 95% confidence regions\n"
                 "(NaCl primary dataset, rows 57-748; Methods A, B are consistent alternatives)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 90)
    print("Step 8: AR(1) autocorrelation correction, NaCl primary dataset")
    print("=" * 90)

    s3.APPLY_CP = APPLY_CP
    c_int, c_p, y = s3.load_primary_data()
    n = len(y)
    print(f"APPLY_CP = {APPLY_CP}")

    # --- OLS (i.i.d.) baseline ---
    fr_ols = step2.fit_model(step2.f_corrected, c_int, c_p, y, beta0=step2.BETA0,
                              label="OLS (i.i.d.)")
    r_ols = y - step2.f_corrected(fr_ols.beta, c_int, c_p)
    phi_ols = lag1_autocorr(r_ols)
    print(f"n = {n}")
    print(f"OLS fit: |chi| = {fr_ols.beta[0]:.5f}, delta* = {fr_ols.beta[1]:.6f}, "
          f"SSE = {fr_ols.sse:.4f}")
    print(f"Lag-1 autocorrelation of OLS residuals: phi_hat = {phi_ols:.4f}")

    # --- Cochrane-Orcutt GLS refit ---
    theta_gls, phi_hat, n_iter, history = cochrane_orcutt(fr_ols.beta, c_int, c_p, y)
    fr_gls = whitened_fit_result(theta_gls, phi_hat, c_int, c_p, y)
    print(f"Cochrane-Orcutt converged in {n_iter} iterations: "
          f"phi_hat = {phi_hat:.4f}")
    print(f"GLS fit: |chi| = {fr_gls.beta[0]:.5f}, delta* = {fr_gls.beta[1]:.6f}, "
          f"SSE_whitened = {fr_gls.sse:.4f}")

    n_eff = n * (1 - phi_hat) / (1 + phi_hat)
    inflation = np.sqrt((1 + phi_hat) / (1 - phi_hat))
    print(f"n_eff = {n_eff:.2f}  (n_eff/n = {n_eff/n:.4f})")
    print(f"Naive CI-inflation factor sqrt((1+phi)/(1-phi)) = {inflation:.3f}")

    # --- Point-estimate shift OLS -> GLS ---
    d_chi_pct = 100 * (fr_gls.beta[0] - fr_ols.beta[0]) / fr_ols.beta[0]
    d_delta_pct = 100 * (fr_gls.beta[1] - fr_ols.beta[1]) / fr_ols.beta[1]
    print(f"Point-estimate shift OLS->GLS: |chi| {d_chi_pct:+.3f}%, "
          f"delta* {d_delta_pct:+.3f}%")

    # --- i.i.d. Wald & profile CIs (reuse s3 unmodified) ---
    z = norm.ppf(1 - ALPHA / 2)
    wald_ols = [(fr_ols.beta[k] - z * fr_ols.se[k], fr_ols.beta[k] + z * fr_ols.se[k])
                for k in range(2)]
    profile_ols = []
    for idx in range(2):
        grid, prof_sse, thresh = s3.profile_likelihood(
            fr_ols.beta, fr_ols.se, fr_ols.sse, fr_ols.n, fr_ols.p, c_int, c_p, y, idx,
            alpha=ALPHA,
        )
        profile_ols.append(profile_ci_from_curve(grid, prof_sse, thresh))

    # --- Method A (GLS/whitening): Wald & profile CIs, ACTUAL n-p df ---
    wald_A = [(fr_gls.beta[k] - z * fr_gls.se[k], fr_gls.beta[k] + z * fr_gls.se[k])
              for k in range(2)]
    profile_A = []
    for idx in range(2):
        grid, prof_sse, thresh = profile_likelihood_whitened(
            fr_gls.beta, fr_gls.se, fr_gls.sse, n, fr_gls.p, phi_hat, c_int, c_p, y,
            idx, alpha=ALPHA,
        )
        profile_A.append(profile_ci_from_curve(grid, prof_sse, thresh))

    # --- Method B (n_eff heuristic): no whitening, OLS fit/SSE, n_eff-p df ---
    se_B = fr_ols.se * inflation
    beta_B = fr_ols.beta  # Method B does not refit theta -- same OLS point estimate
    wald_B = [(beta_B[k] - z * se_B[k], beta_B[k] + z * se_B[k]) for k in range(2)]
    profile_B = []
    for idx in range(2):
        grid, prof_sse, thresh = s3.profile_likelihood(
            beta_B, se_B, fr_ols.sse, n_eff, fr_ols.p, c_int, c_p, y, idx, alpha=ALPHA,
        )
        profile_B.append(profile_ci_from_curve(grid, prof_sse, thresh))

    names = ["|chi| (mM)", "delta*"]
    print()
    print("=" * 100)
    print("i.i.d. vs. METHOD A (GLS) vs. METHOD B (n_eff) -- 95% CIs")
    print("=" * 100)
    header = (f"{'Param':<12}{'point':>10}"
              f"{'Wald (iid)':>20}{'Wald (A)':>20}{'Wald (B)':>20}{'foldA':>7}{'foldB':>7}")
    print(header)
    fold_widths = {}
    for idx, name in enumerate(names):
        w_o, w_a, w_b = wald_ols[idx], wald_A[idx], wald_B[idx]
        fold_a = (w_a[1] - w_a[0]) / (w_o[1] - w_o[0])
        fold_b = (w_b[1] - w_b[0]) / (w_o[1] - w_o[0])
        print(f"{name:<12}{fr_ols.beta[idx]:>10.4f}"
              f"  ({w_o[0]:6.4f},{w_o[1]:6.4f})"
              f"  ({w_a[0]:6.4f},{w_a[1]:6.4f})"
              f"  ({w_b[0]:6.4f},{w_b[1]:6.4f})"
              f"{fold_a:>7.2f}{fold_b:>7.2f}")
    print()
    header2 = (f"{'Param':<12}{'point':>10}"
               f"{'Profile (iid)':>20}{'Profile (A)':>20}{'Profile (B)':>20}{'foldA':>7}{'foldB':>7}")
    print(header2)
    for idx, name in enumerate(names):
        p_o, p_a, p_b = profile_ols[idx], profile_A[idx], profile_B[idx]
        fold_a = (p_a[1] - p_a[0]) / (p_o[1] - p_o[0])
        fold_b = (p_b[1] - p_b[0]) / (p_o[1] - p_o[0])
        fold_widths[name] = (fold_a, fold_b)
        print(f"{name:<12}{fr_ols.beta[idx]:>10.4f}"
              f"  ({p_o[0]:6.4f},{p_o[1]:6.4f})"
              f"  ({p_a[0]:6.4f},{p_a[1]:6.4f})"
              f"  ({p_b[0]:6.4f},{p_b[1]:6.4f})"
              f"{fold_a:>7.2f}{fold_b:>7.2f}")
    print()
    print("Profile/Wald ratio, per method (checks 'profile tracks Wald' i.e. mild nonlinearity):")
    for idx, name in enumerate(names):
        r_iid = (profile_ols[idx][1] - profile_ols[idx][0]) / (wald_ols[idx][1] - wald_ols[idx][0])
        r_a = (profile_A[idx][1] - profile_A[idx][0]) / (wald_A[idx][1] - wald_A[idx][0])
        r_b = (profile_B[idx][1] - profile_B[idx][0]) / (wald_B[idx][1] - wald_B[idx][0])
        print(f"  {name:<12} iid={r_iid:.3f}  A={r_a:.3f}  B={r_b:.3f}")
    print()

    # --- ACF sufficiency check ---
    r_final = y - step2.f_corrected(fr_gls.beta, c_int, c_p)
    r_white_final = whitened_resid_fn(fr_gls.beta, phi_hat, c_int, c_p, y)
    acf_raw = acf(r_ols)
    acf_white = acf(r_white_final)
    band = 1.96 / np.sqrt(n)
    n_exceed_raw = int(np.sum(np.abs(acf_raw[1:]) > band))
    n_exceed_white = int(np.sum(np.abs(acf_white[1:]) > band))
    print("=" * 90)
    print("ACF SUFFICIENCY CHECK")
    print("=" * 90)
    print(f"+/-1.96/sqrt(n) band = +/-{band:.4f}")
    print(f"OLS residual ACF: lag-1 = {acf_raw[1]:.4f}; "
          f"{n_exceed_raw}/{len(acf_raw)-1} lags exceed the band")
    print(f"Whitened residual ACF: lag-1 = {acf_white[1]:.4f}; "
          f"{n_exceed_white}/{len(acf_white)-1} lags exceed the band")
    print()

    # --- Figures ---
    plot_acf(acf_raw, acf_white, n, FIG_DIR / "nacl_ar1_acf.png")
    plot_regions(fr_ols, fr_gls, phi_hat, n, n_eff, c_int, c_p, y,
                 FIG_DIR / "nacl_ar1_regions.png")

    return dict(
        fr_ols=fr_ols, fr_gls=fr_gls, phi_hat=phi_hat, n=n, n_eff=n_eff,
        inflation=inflation, wald_ols=wald_ols, wald_A=wald_A, wald_B=wald_B,
        profile_ols=profile_ols, profile_A=profile_A, profile_B=profile_B,
        fold_widths=fold_widths,
        acf_raw=acf_raw, acf_white=acf_white, n_exceed_raw=n_exceed_raw,
        n_exceed_white=n_exceed_white,
    )


if __name__ == "__main__":
    main()
