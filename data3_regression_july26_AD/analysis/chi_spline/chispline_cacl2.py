"""
Step 12 (first cut): concentration-dependent chi(c) spline regression, CaCl2.

Motivation (docs/prompts/step12_chispline_firstcut.md, math in
docs/reports/draft_chispline_futurework.tex): the CaCl2 adsorption-hypothesis
diagnostics (analysis/cacl2/regress_cacl2.py, report sec:cacl2diag) show a
sliding-window "effective |chi|" that drifts with concentration -- consistent
with concentration-dependent counter-ion adsorption, but also with an
artifact of forcing a single constant charge onto a mis-specified mean model.
This script asks the narrow question: is that drift a REAL, IDENTIFIABLE
chi(c), or does it evaporate under an honest joint chi(c)+K fit?

Model
-----
Same CaCl2 2:1 co-ion cubic and transformed response as
analysis/cacl2/regress_cacl2.py, with the scalar |chi| replaced by a smooth
chi(c) = B(c)^T beta, a clamped cubic B-spline evaluated separately at each
face (c_int, c_p) -- draft section Eqs. (2)-(3). K stays a single scalar
(same lumped-partition-constant convention as the constant-chi model). Cubic
B-spline basis functions are non-negative and sum to 1 (partition of unity)
on the clamped domain, so:
  - constraining all control points beta_k >= 0 guarantees chi(c) >= 0 with
    no exp() link needed;
  - a spline with ALL beta_k equal to some c is IDENTICAL, pointwise, to the
    constant model chi(c)=c -- so the existing constant-chi fit
    (analysis/cacl2/regress_cacl2.py, reused unmodified) is the exact
    0-interior-knot nested submodel, and F-test/AIC/BIC model comparison
    against it is valid.

Scope (deliberately narrow, per the prompt): CaCl2 MC5 05.27.26 only (the
widest concentration range, ~1645 points to ~108 mM); scipy
least_squares only (no Pyomo/ParmEst); no roughness penalty (<=2 interior
knots keeps this identifiable without one); an optional NaCl null-check
(does the same machinery manufacture spurious drift on a dataset with no
adsorption hypothesis?).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/chi_spline/chispline_cacl2.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import BSpline
from scipy.optimize import least_squares
from scipy.stats import f as f_dist

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "cacl2"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
import plot_style  # noqa: E402
import regress_cacl2 as cacl2  # noqa: E402
import regress_nacl as step2  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
DEGREE = 3  # cubic B-spline
N_INTERIOR_LIST = [0, 1, 2]
RNG_SEED = 20260706
N_MULTISTART = 30  # flexible-fit convergence check (prompt: "confirm
                    # multi-start/convergence robustness on the flexible fits")

CACL2_DATASET = {"key": "MC5-0527", "sheet": "NF270_MC5 05.27.26_CaCl2",
                  "label": "MC5 (05.27.26)"}


# ---------------------------------------------------------------------------
# Clamped cubic B-spline basis: knots at data quantiles, no extrapolation
# ---------------------------------------------------------------------------
def build_knots(c_data, n_interior, degree=DEGREE):
    c_min, c_max = float(np.min(c_data)), float(np.max(c_data))
    if n_interior == 0:
        interior = np.array([])
    else:
        qs = np.linspace(0.0, 1.0, n_interior + 2)[1:-1]
        interior = np.quantile(c_data, qs)
    knots = np.concatenate((
        np.full(degree + 1, c_min),
        interior,
        np.full(degree + 1, c_max),
    ))
    n_basis = len(knots) - degree - 1
    return knots, n_basis, c_min, c_max


def eval_basis(knots, n_basis, c_min, c_max, c, degree=DEGREE):
    """Design matrix B(c), shape (len(c), n_basis). Clamped knots give a
    partition of unity (rows sum to 1) on [c_min, c_max]; nudge points
    landing exactly on c_max inward by a hair to sidestep scipy's half-open
    knot-span convention (BSpline is left-inclusive/right-exclusive per
    span), which would otherwise evaluate the rightmost basis function as 0
    exactly at c_max."""
    c = np.atleast_1d(np.asarray(c, dtype=float))
    span = c_max - c_min
    c_eval = np.clip(c, c_min, c_max - 1e-9 * span if span > 0 else c_max)
    B = np.empty((len(c), n_basis))
    for k in range(n_basis):
        coef = np.zeros(n_basis)
        coef[k] = 1.0
        B[:, k] = BSpline(knots, coef, degree, extrapolate=False)(c_eval)
    B = np.nan_to_num(B, nan=0.0)
    return B


def check_partition_of_unity(knots, n_basis, c_min, c_max, degree=DEGREE):
    c_test = np.linspace(c_min, c_max, 50)
    B = eval_basis(knots, n_basis, c_min, c_max, c_test, degree)
    row_sums = B.sum(axis=1)
    assert np.allclose(row_sums, 1.0, atol=1e-8), \
        f"B-spline basis is not a partition of unity (max |sum-1|={np.max(np.abs(row_sums-1)):.2e})"


# ---------------------------------------------------------------------------
# CaCl2 2:1 co-ion cubic with a per-point (spline-varying) chi
# ---------------------------------------------------------------------------
def co_root_varchi(chi_arr, K, c_s):
    """Same closed-form positive real root as regress_cacl2.co_root
    (c_m^3 + chi*c_m^2 - K*c_s^3 = 0), generalized to accept a per-point
    chi_arr (same shape as c_s) instead of a single scalar chi -- required
    once chi is a function of concentration evaluated at each face. K stays
    a single scalar. Verified against regress_cacl2.co_root to <1e-9
    relative error by passing chi_arr = np.full_like(c_s, chi) (see
    validate_against_scalar_coroot below)."""
    c_s = np.atleast_1d(np.asarray(c_s, dtype=float))
    chi_arr = np.atleast_1d(np.asarray(chi_arr, dtype=float))
    chi_arr = np.broadcast_to(chi_arr, c_s.shape).astype(float).copy()
    rhs = K * c_s ** 3

    m = np.where(chi_arr > 0, rhs / np.maximum(chi_arr, 1e-30) ** 3, np.inf)
    small = m < 1e-6

    out = np.full_like(c_s, np.nan, dtype=float)

    if np.any(small):
        t0 = np.sqrt(rhs[small] / chi_arr[small])
        for _ in range(3):
            g = t0 ** 3 + chi_arr[small] * t0 ** 2 - rhs[small]
            gp = 3 * t0 ** 2 + 2 * chi_arr[small] * t0
            t0 = t0 - g / np.where(gp != 0, gp, 1.0)
        out[small] = t0

    if np.any(~small):
        cs2 = c_s[~small]
        chi2 = chi_arr[~small]
        a = 1.0
        b = chi2.astype(complex)
        c = np.zeros_like(cs2, dtype=complex)
        d = (-K * cs2 ** 3).astype(complex)

        p = (3 * a * c - b ** 2) / (3 * a ** 2)
        q = (2 * b ** 3 - 9 * a * b * c + 27 * a ** 2 * d) / (27 * a ** 3)
        disc = (q / 2) ** 2 + (p / 3) ** 3
        sqrt_disc = np.sqrt(disc)
        u = (-q / 2 + sqrt_disc) ** (1 / 3 + 0j)
        u_safe = np.where(np.abs(u) < 1e-14, 1e-14, u)
        v = np.where(np.abs(u) < 1e-14, 0.0, -p / (3 * u_safe))
        omega = np.exp(2j * np.pi / 3)
        x1, x2, x3 = u + v, omega * u + omega ** 2 * v, omega ** 2 * u + omega * v
        shift = b / (3 * a)
        ts = np.stack([x1 - shift, x2 - shift, x3 - shift], axis=0)
        tol = 1e-6 * (np.abs(cs2).max() + 1.0) if cs2.size else 1e-6
        is_real_pos = (np.abs(ts.imag) < tol) & (ts.real > 1e-10)
        real_part = np.where(is_real_pos, ts.real, -np.inf)
        r = np.max(real_part, axis=0)
        out[~small] = np.where(np.isfinite(r), r, np.nan)

    return out


def validate_against_scalar_coroot(rng, n=2000):
    """Sanity check: co_root_varchi with a constant chi_arr must reproduce
    regress_cacl2.co_root (scalar chi) to high precision."""
    chi = rng.uniform(1e-3, 1e5, n)
    K = rng.uniform(1e-8, 1e3, n)
    c_s = rng.uniform(1e-2, 1e3, n)
    max_rel_err = 0.0
    for chi_i, K_i, c_i in zip(chi[:200], K[:200], c_s[:200]):
        r_scalar = cacl2.co_root((chi_i, K_i), np.array([c_i]))[0]
        r_vec = co_root_varchi(np.array([chi_i]), K_i, np.array([c_i]))[0]
        if np.isfinite(r_scalar) and r_scalar > 0:
            max_rel_err = max(max_rel_err, abs(r_vec - r_scalar) / r_scalar)
    return max_rel_err


def f_chispline(beta_all, c_int, c_p, knots, n_basis, c_min, c_max):
    """beta_all = [beta_1..beta_n_basis, K]. chi(c) evaluated separately at
    each face; response unchanged (same Cl- stoichiometry, 2x conc., as
    regress_cacl2.f_cacl2)."""
    betas = beta_all[:-1]
    K = beta_all[-1]
    B_int = eval_basis(knots, n_basis, c_min, c_max, c_int)
    B_p = eval_basis(knots, n_basis, c_min, c_max, c_p)
    chi_int = B_int @ betas
    chi_p = B_p @ betas
    c_co_int = co_root_varchi(chi_int, K, 2.0 * np.asarray(c_int))
    c_co_p = co_root_varchi(chi_p, K, 2.0 * np.asarray(c_p))
    return c_co_int - c_co_p


# ---------------------------------------------------------------------------
# Fitting
# ---------------------------------------------------------------------------
def fit_chispline(c_int, c_p, y, n_interior, chi0_flat, K0, rng,
                   n_multistart=N_MULTISTART):
    knots, n_basis, c_min, c_max = build_knots(c_int, n_interior)
    check_partition_of_unity(knots, n_basis, c_min, c_max)

    def resid(beta_all):
        return f_chispline(beta_all, c_int, c_p, knots, n_basis, c_min, c_max) - y

    beta0 = np.concatenate([np.full(n_basis, chi0_flat), [K0]])
    lo = np.zeros(n_basis + 1)
    hi = np.full(n_basis + 1, np.inf)
    res = least_squares(resid, beta0, method="trf", bounds=(lo, hi))
    beta_hat = res.x
    sse = float(np.sum(res.fun ** 2))
    n = len(y)
    p = len(beta_hat)
    dof = max(n - p, 1)
    J = res.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma2 = sse / dof
    cov = sigma2 * JTJ_inv

    # Multi-start robustness: perturb the warm start broadly (log-uniform
    # jitter per control point + K), count fraction reaching the same SSE.
    n_global = 0
    for _ in range(n_multistart):
        jitter_beta = beta0[:-1] * np.exp(rng.uniform(-2.0, 2.0, n_basis))
        jitter_K = K0 * np.exp(rng.uniform(-2.0, 2.0))
        beta0_try = np.concatenate([jitter_beta, [jitter_K]])
        try:
            res_try = least_squares(resid, beta0_try, method="trf", bounds=(lo, hi))
            sse_try = float(np.sum(res_try.fun ** 2))
            if res_try.success and sse_try <= sse * (1 + 1e-4) + 1e-6:
                n_global += 1
        except Exception:
            pass
    pct_global = 100.0 * n_global / n_multistart

    return dict(
        n_interior=n_interior, knots=knots, n_basis=n_basis,
        c_min=c_min, c_max=c_max, beta=beta_hat, betas=beta_hat[:-1],
        K=beta_hat[-1], sse=sse, n=n, p=p, cov=cov, pct_global=pct_global,
    )


def pseudo_r2(sse, y):
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - sse / ss_tot


def aic_bic(sse, n, p):
    aic = n * np.log(sse / n) + 2 * p
    bic = n * np.log(sse / n) + p * np.log(n)
    return aic, bic


def nested_f_test(sse_reduced, p_reduced, sse_full, p_full, n, alpha=ALPHA):
    df1 = p_full - p_reduced
    df2 = n - p_full
    if df1 <= 0 or sse_full >= sse_reduced:
        return dict(F=0.0, df1=df1, df2=df2, p_value=1.0, Fcrit=np.nan, significant=False)
    F = ((sse_reduced - sse_full) / df1) / (sse_full / df2)
    p_value = float(f_dist.sf(F, df1, df2))
    Fcrit = float(f_dist.ppf(1 - alpha, df1, df2))
    return dict(F=float(F), df1=df1, df2=df2, p_value=p_value, Fcrit=Fcrit,
                significant=bool(F > Fcrit))


def chi_band(fit, c_grid):
    """Pointwise chi(c) +/- z*SE band from the Gauss-Newton covariance
    (draft section Eq. (6)): Var(chi_hat(c)) = B(c)^T Cov(beta_hat) B(c),
    restricted to the beta (control-point) block of the full (beta, K)
    covariance."""
    B = eval_basis(fit["knots"], fit["n_basis"], fit["c_min"], fit["c_max"], c_grid)
    chi_hat = B @ fit["betas"]
    cov_beta = fit["cov"][:fit["n_basis"], :fit["n_basis"]]
    var_chi = np.einsum("ij,jk,ik->i", B, cov_beta, B)
    var_chi = np.clip(var_chi, 0.0, None)
    se_chi = np.sqrt(var_chi)
    return chi_hat, se_chi


# ---------------------------------------------------------------------------
# NaCl null-check (1:1, algebraic root -- no cubic solve needed)
# ---------------------------------------------------------------------------
def f_nacl_chispline(beta_all, c_int, c_p, knots, n_basis, c_min, c_max):
    """beta_all = [beta_1..beta_n_basis, delta_star] -- delta* is jointly
    refit alongside the spline control points, exactly as K is jointly
    refit in the CaCl2 model (f_chispline above); holding delta* fixed at
    the constant-model value while only refitting beta_k would give the
    spline extra unearned degrees of freedom relative to that fixed
    parameter and bias this null-check toward finding spurious drift."""
    betas = beta_all[:-1]
    delta = beta_all[-1]
    B_int = eval_basis(knots, n_basis, c_min, c_max, c_int)
    B_p = eval_basis(knots, n_basis, c_min, c_max, c_p)
    chi_int = B_int @ betas
    chi_p = B_p @ betas
    return 0.5 * (np.sqrt(chi_int ** 2 + 4.0 * delta * np.asarray(c_int) ** 2)
                  - np.sqrt(chi_p ** 2 + 4.0 * delta * np.asarray(c_p) ** 2))


def run_nacl_null_check():
    print("-" * 100)
    print("Optional null-check: same 1-knot spline machinery on NaCl MC5 05.27.26")
    print("-" * 100)
    df_full = step2.load_sheet()  # default: "(2)" sheet, our-CP (APPLY_CP=True)
    df_full = df_full.loc[(df_full.index >= step2.MATLAB_FIRST_ROW)
                           & (df_full.index <= step2.MATLAB_LAST_ROW)].dropna()
    df_full = step2.build_deltaC_m(df_full)
    c_int = df_full["c_int"].to_numpy()
    c_p = df_full["c_p"].to_numpy()
    y = df_full["deltaC_m"].to_numpy()

    fr_const = step2.fit_model(step2.f_corrected, c_int, c_p, y, beta0=step2.BETA0,
                                label="NaCl MC5 05.27.26 constant chi")
    chi0, delta_hat = fr_const.beta
    sse_const = fr_const.sse
    p_const = fr_const.p

    knots, n_basis, c_min, c_max = build_knots(c_int, n_interior=1)
    check_partition_of_unity(knots, n_basis, c_min, c_max)

    def resid(beta_all):
        return f_nacl_chispline(beta_all, c_int, c_p, knots, n_basis, c_min, c_max) - y

    beta0_arr = np.concatenate([np.full(n_basis, chi0), [delta_hat]])
    lo = np.zeros(n_basis + 1)
    hi = np.full(n_basis + 1, np.inf)
    res = least_squares(resid, beta0_arr, method="trf", bounds=(lo, hi))
    sse_spline = float(np.sum(res.fun ** 2))
    p_spline = n_basis + 1  # betas jointly refit WITH delta*, same convention as CaCl2

    ftest = nested_f_test(sse_const, p_const, sse_spline, p_spline, len(y))
    aic_c, bic_c = aic_bic(sse_const, len(y), p_const)
    aic_s, bic_s = aic_bic(sse_spline, len(y), p_spline)

    betas_hat = res.x[:-1]
    delta_spline = res.x[-1]
    beta_spread = float(np.max(betas_hat) - np.min(betas_hat))
    beta_rel_spread = beta_spread / chi0 if chi0 else np.nan
    print(f"n = {len(y)}, constant chi = {chi0:.4f} mM, delta* = {delta_hat:.5f}")
    print(f"Constant model:  SSE = {sse_const:.4f}, p = {p_const}, "
          f"AIC = {aic_c:.2f}, BIC = {bic_c:.2f}")
    print(f"1-knot spline:   SSE = {sse_spline:.4f}, p (incl. jointly-refit delta*) = {p_spline}, "
          f"AIC = {aic_s:.2f}, BIC = {bic_s:.2f}, delta*_hat = {delta_spline:.5f}")
    print(f"  control points beta_k = {np.array2string(betas_hat, precision=3)}")
    print(f"  spread (max-min)/chi0 = {beta_rel_spread:.4f} "
          f"({'FLAT (as expected: no spurious drift)' if beta_rel_spread < 0.05 else 'NOT flat -- unexpected'})")
    print(f"  F-test vs. constant: F={ftest['F']:.3f} (crit {ftest['Fcrit']:.3f}), "
          f"p={ftest['p_value']:.4f} "
          f"({'SIGNIFICANT' if ftest['significant'] else 'not significant'})")
    print(f"  AIC improvement (spline - constant): {aic_s - aic_c:+.2f} "
          f"({'improves' if aic_s < aic_c else 'does NOT improve'})")
    print()
    return dict(sse_const=sse_const, sse_spline=sse_spline, ftest=ftest,
                aic_c=aic_c, aic_s=aic_s, beta_rel_spread=beta_rel_spread,
                betas=betas_hat, delta_spline=delta_spline)


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------
def make_figure(df, fr_const, fits, window, outfile):
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    order = np.argsort(c_int)

    best_fit = fits[max(N_INTERIOR_LIST)]  # most flexible model for panel (a)/(b)

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 2, figsize=(w, w * 0.75))

    ax = axes[0]
    ax.plot(c_int, y, "o", color=plot_style.OKABE_ITO["sky_blue"], markersize=3,
            alpha=0.35, label="data")
    y_const = cacl2.f_cacl2(fr_const.beta, c_int[order], c_p[order])
    ax.plot(c_int[order], y_const, "--", color="k", lw=2,
            label=f"constant $\\chi$ ({fr_const.beta[0]:.1f} mM)")
    y_spline = f_chispline(best_fit["beta"], c_int[order], c_p[order],
                            best_fit["knots"], best_fit["n_basis"],
                            best_fit["c_min"], best_fit["c_max"])
    ax.plot(c_int[order], y_spline, "-", color=plot_style.OKABE_ITO["vermillion"],
            lw=2, label=f"{max(N_INTERIOR_LIST)}-knot $\\chi(c)$ spline")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")
    ax.set_title("(a) fit vs. data")
    ax.set_box_aspect(1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), ncol=1)

    ax = axes[1]
    c_grid = np.linspace(best_fit["c_min"], best_fit["c_max"], 300)
    ax.axhline(fr_const.beta[0], color="k", ls="--", lw=1.5, label="constant $\\chi$")
    for n_int, fit in fits.items():
        chi_hat, se_chi = chi_band(fit, c_grid)
        color = {1: plot_style.OKABE_ITO["blue"],
                 2: plot_style.OKABE_ITO["vermillion"]}[n_int]
        ax.plot(c_grid, chi_hat, "-", color=color, lw=2,
                label=f"{n_int}-knot spline")
        ax.fill_between(c_grid, chi_hat - 1.96 * se_chi, chi_hat + 1.96 * se_chi,
                         color=color, alpha=0.15)

    centers, locals_, ses, reliable = window
    if np.any(reliable):
        ax.errorbar(centers[reliable], locals_[reliable], yerr=ses[reliable],
                     fmt="o", color=plot_style.OKABE_ITO["bluish_green"],
                     markersize=4, capsize=2, elinewidth=1,
                     label="sliding-window effective $|\\chi|$")
    if np.any(~reliable):
        ax.plot(centers[~reliable], locals_[~reliable], "x",
                color=plot_style.OKABE_ITO["bluish_green"], markersize=5, alpha=0.35)

    ax.set_xlabel(r"$\boldsymbol{c}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\chi(c)}$ [mM]")
    ax.set_title("(b) regressed $\\chi(c)$ vs. sliding-window diagnostic")
    ax.set_box_aspect(1)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2)

    fig.suptitle("CaCl2 $\\boldsymbol{\\chi(c)}$ spline, first cut (MC5 05.27.26, CP-corrected)")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 100)
    print("Step 12 (first cut): CaCl2 chi(c) spline regression, MC5 05.27.26")
    print("=" * 100)
    print(f"cacl2.APPLY_CP = {cacl2.APPLY_CP}")

    rng = np.random.default_rng(RNG_SEED)

    max_rel_err = validate_against_scalar_coroot(rng)
    print(f"co_root_varchi vs. co_root (scalar chi) sanity check: "
          f"max rel. err = {max_rel_err:.2e} (should be ~1e-9 or smaller)")

    df = cacl2.load_dataset(CACL2_DATASET)
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    n = len(y)
    print(f"n = {n}, c_int range = [{c_int.min():.3f}, {c_int.max():.3f}] mM")

    fr_const = cacl2.fit_model_bounded(cacl2.f_cacl2, c_int, c_p, y,
                                        beta0=cacl2.BETA0,
                                        label="MC5 05.27.26 constant chi")
    print(f"Constant-chi fit: |chi| = {fr_const.beta[0]:.4f} mM, K = {fr_const.beta[1]:.5f}, "
          f"SSE = {fr_const.sse:.4f}")

    fits = {}
    for n_interior in [1, 2]:
        fit = fit_chispline(c_int, c_p, y, n_interior,
                             chi0_flat=fr_const.beta[0], K0=fr_const.beta[1],
                             rng=rng)
        fits[n_interior] = fit

    print()
    print("-" * 100)
    print("Model comparison: constant vs. 1-knot vs. 2-knot")
    print("-" * 100)
    header = f"{'model':<16}{'n_basis+K':>10}{'K_hat':>10}{'SSE':>12}{'pseudo-R2':>11}{'AIC':>10}{'BIC':>10}{'MS % global':>13}"
    print(header)
    r2_const = pseudo_r2(fr_const.sse, y)
    aic_const, bic_const = aic_bic(fr_const.sse, n, fr_const.p)
    print(f"{'constant':<16}{fr_const.p:>10}{fr_const.beta[1]:>10.5f}{fr_const.sse:>12.4f}"
          f"{r2_const:>11.5f}{aic_const:>10.2f}{bic_const:>10.2f}{'--':>13}")
    for n_interior in [1, 2]:
        fit = fits[n_interior]
        r2 = pseudo_r2(fit["sse"], y)
        aic, bic = aic_bic(fit["sse"], n, fit["p"])
        label = f"{n_interior}-knot spline"
        print(f"{label:<16}{fit['p']:>10}{fit['K']:>10.5f}{fit['sse']:>12.4f}"
              f"{r2:>11.5f}{aic:>10.2f}{bic:>10.2f}{fit['pct_global']:>12.1f}%")
    print()

    print("Nested F-tests vs. constant-chi baseline:")
    for n_interior in [1, 2]:
        fit = fits[n_interior]
        ftest = nested_f_test(fr_const.sse, fr_const.p, fit["sse"], fit["p"], n)
        print(f"  constant -> {n_interior}-knot: F({ftest['df1']},{ftest['df2']}) = "
              f"{ftest['F']:.3f} (crit {ftest['Fcrit']:.3f} at alpha={ALPHA}), "
              f"p={ftest['p_value']:.5f} -> "
              f"{'flexibility JUSTIFIED' if ftest['significant'] else 'NOT justified (fitting noise)'}")
    print()
    ftest_1v2 = nested_f_test(fits[1]["sse"], fits[1]["p"], fits[2]["sse"], fits[2]["p"], n)
    print(f"  1-knot -> 2-knot: F({ftest_1v2['df1']},{ftest_1v2['df2']}) = "
          f"{ftest_1v2['F']:.3f} (crit {ftest_1v2['Fcrit']:.3f}), "
          f"p={ftest_1v2['p_value']:.5f} -> "
          f"{'2nd knot JUSTIFIED' if ftest_1v2['significant'] else 'NOT justified'}")
    print()

    for n_interior in [1, 2]:
        fit = fits[n_interior]
        print(f"{n_interior}-knot control points beta_k = "
              f"{np.array2string(fit['betas'], precision=3)} mM")
    print()

    # Sliding-window effective-chi diagnostic (reused unmodified)
    centers, chi_locals, chi_ses, reliable = cacl2.effective_chi_vs_conc(
        c_int, c_p, y, K_fixed=fr_const.beta[1], chi0=fr_const.beta[0],
    )
    window = (centers, chi_locals, chi_ses, reliable)

    print("-" * 100)
    print("Cross-check: does the regressed chi(c_int) track the sliding-window diagnostic?")
    print("-" * 100)
    for n_interior in [1, 2]:
        fit = fits[n_interior]
        chi_at_centers, _ = chi_band(fit, centers[reliable])
        if np.any(reliable) and chi_at_centers.size:
            resid_vs_window = chi_at_centers - chi_locals[reliable]
            rmse_vs_window = float(np.sqrt(np.mean(resid_vs_window ** 2)))
            corr = float(np.corrcoef(chi_at_centers, chi_locals[reliable])[0, 1]) \
                if len(chi_at_centers) > 1 else np.nan
            print(f"  {n_interior}-knot spline vs. reliable windows (n={reliable.sum()}): "
                  f"RMSE = {rmse_vs_window:.3f} mM, corr = {corr:.3f}")
    print()

    run_nacl_null_check()

    make_figure(df, fr_const, fits, window, FIG_DIR / "cacl2_chispline_firstcut.png")

    return dict(fr_const=fr_const, fits=fits, window=window)


if __name__ == "__main__":
    main()
