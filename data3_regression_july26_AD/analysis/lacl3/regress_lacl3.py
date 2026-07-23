"""
LaCl3 (3:1) Donnan regression + adsorption-hypothesis diagnostics (report Sec. 4).

Method A (this script): the 3:1 co-ion quartic has NO compact closed form (unlike
NaCl's quadratic or CaCl2's cubic), so it is solved numerically per point via a
bracketed root-find (Brent) -- verified against numpy.roots on a 3000-point
random parameter scan to <1e-6 relative error, ~2.5x faster than a per-point
numpy.roots call at this dataset's scale (n=623). Method B
(analysis/lacl3/parmest_lacl3.py) instead poses the quartic as a Pyomo
nonlinear CONSTRAINT and estimates via ParmEst/Ipopt -- the formulation that
generalizes to the multicomponent case (no closed form there either), used
here as a cross-check on Method A.

Model
-----
Working hypothesis (as for CaCl2, only more so): La3+ is expected to adsorb
onto the membrane even more strongly than Ca2+, so we anticipate an even
poorer fit / stronger residual structure / |chi| pinned near zero.

Co-ion = Cl-, solution concentration set by stoichiometry: c_co_s = 3*c_LaCl3
(3 Cl- per LaCl3), used for both c_int and c_p.

Donnan co-ion quartic (report Eq. 31co), lumped partition constant K (>= 0)
absorbing delta_star*(delta_co_circ)^2 and constant stoichiometric factors:
    c_co_m^4 + |chi|*c_co_m^3 - K*c_co_s^4 = 0
Exactly one positive real root (Descartes: sign pattern +,+,0,0,-).

Transformed response (Cl- flux form, analogous to corrected NaCl/CaCl2):
    deltaC_m = J_s_Cl * l / D_m
    J_s_Cl = 3 * c_p * J_w         (Cl- flux = 3x the LaCl3 molar flux)
    l = 80e-9 m
    D_m = 4*D_La_m*D_Cl_m / (3*D_La_m + D_Cl_m)   (3:1 ambipolar membrane
                                                     diffusivity)
    D_La_m = 0.626e-12 m^2/s, D_Cl_m = 2.03e-12 m^2/s

Prediction = c_co_m(c_int) - c_co_m(c_p). Fit (|chi|, K) to deltaC_m by
scipy.optimize.least_squares, bounded chi,K >= 0 (as for CaCl2 -- chi enters
the quartic unsquared/asymmetrically, so an unconstrained fit can wander
negative with no physical meaning).

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/lacl3/regress_lacl3.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import brentq, least_squares

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import cp_authors_bulk  # noqa: E402

plot_style.apply_style()

DATA_FILE = REPO_ROOT / "data" / "BoE Analysis.xlsx"
SHEET_NAME = "NF270_MC2 05.21.24_LaCl3"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): the LaCl3 sheet is the
# authors' non-CP bulk convention (Phase 1 validation), so CP is applied on
# top of column B. Default primary run is CP=on; APPLY_CP=False reproduces
# the pre-Phase-2 numbers byte-for-byte.
APPLY_CP = True

RNG_SEED = 20260705
N_MULTISTART = 100
ALPHA = 0.05

# Physical constants (membrane-phase diffusivities, m^2/s; report Sec. 4)
D_LA_M = 0.626e-12
D_CL_M = 2.03e-12
L_MEMBRANE = 80e-9  # m
D_LACL3_M = 4.0 * D_LA_M * D_CL_M / (3.0 * D_LA_M + D_CL_M)

BETA0 = np.array([44.0, 0.01])  # initial guess: (|chi| [mM], K [-])
FIT_BOUNDS = ([0.0, 0.0], [np.inf, np.inf])

CHI_LO, CHI_HI = 1.0, 500.0
K_LO, K_HI = 1e-6, 10.0

N_WINDOW_POINTS = 60
WINDOW_STRIDE = 30
WINDOW_CHI_BOUND = 300.0
WINDOW_RELSE_UNRELIABLE = 2.0

DATASET_LABEL = "MC2 (05.21.24)"


# ---------------------------------------------------------------------------
# Data loading (same 15-column schema as NaCl/CaCl2: B=c_int, H=c_p, J=Jw)
# ---------------------------------------------------------------------------
def load_sheet(data_file=DATA_FILE, sheet_name=SHEET_NAME, apply_cp: bool | None = None):
    import pandas as pd
    if apply_cp is None:
        apply_cp = APPLY_CP
    df = pd.read_excel(data_file, sheet_name=sheet_name, header=0)
    col_map = {
        "c_int": df.columns[1],  # column B
        "c_p": df.columns[7],    # column H
        "J_w": df.columns[9],    # column J
    }
    out = df[[col_map["c_int"], col_map["c_p"], col_map["J_w"]]].copy()
    out.columns = ["c_int", "c_p", "J_w"]
    out.index = np.arange(2, 2 + len(out))
    out.index.name = "sheet_row"
    if apply_cp:
        # Column B is the authors' non-CP bulk convention for this sheet
        # (Phase 1 validation) -- apply CP on top of it.
        out["c_int"] = cp_authors_bulk.apply_cp(out["c_int"], out["c_p"], out["J_w"], "LaCl3")
    return out


def build_deltaC_m(df):
    df = df.copy()
    J_s_cl = 3.0 * df["c_p"] * df["J_w"]  # mol/m^2/s (Cl- flux)
    df["J_s_cl"] = J_s_cl
    df["deltaC_m"] = J_s_cl * L_MEMBRANE / D_LACL3_M  # mM
    return df


def load_dataset():
    df_full = load_sheet()
    valid = df_full.notna().all(axis=1) & (df_full["c_int"] > df_full["c_p"]) \
        & (df_full["c_int"] > 0)
    df_valid = df_full[valid].copy()
    return build_deltaC_m(df_valid)


# ---------------------------------------------------------------------------
# Model: 3:1 co-ion quartic, single positive real root via bracketed Brent
# ---------------------------------------------------------------------------
def _quartic_g(c, chi, K, cs):
    return c ** 4 + chi * c ** 3 - K * cs ** 4


def _co_root_scalar(chi, K, cs):
    if cs <= 0:
        return 0.0
    lo, hi = 1e-12, max(cs, 1e-6)
    g_lo = _quartic_g(lo, chi, K, cs)
    if g_lo > 0:
        return lo
    g_hi = _quartic_g(hi, chi, K, cs)
    n_expand = 0
    while g_hi < 0 and n_expand < 40:
        hi *= 2.0
        g_hi = _quartic_g(hi, chi, K, cs)
        n_expand += 1
    if g_hi < 0:
        return np.nan
    try:
        return brentq(_quartic_g, lo, hi, args=(chi, K, cs), xtol=1e-12, rtol=1e-12)
    except ValueError:
        return np.nan


def co_root(beta, c_s):
    """Unique positive real root of c_m^4 + chi*c_m^3 - K*c_s^4 = 0
    (Descartes: sign pattern +,+,0,0,-), via a per-point bracketed Brent
    search (verified against numpy.roots to <1e-6 relative error over a
    3000-point random scan spanning chi in [1e-3,1e3], K in [1e-6,1e2],
    c_s in [1e-2,1e3]; ~2.5x faster than numpy.roots at this dataset's scale,
    n=623 -- see analysis/lacl3/README.md)."""
    chi, K = beta
    c_s = np.atleast_1d(np.asarray(c_s, dtype=float))
    out = np.empty_like(c_s)
    for i, cs in enumerate(c_s):
        out[i] = _co_root_scalar(chi, K, cs)
    return out


def f_lacl3(beta, c_int, c_p):
    """Prediction = c_co_m(c_int) - c_co_m(c_p), co-ion solution conc. =
    3*c_LaCl3 (stoichiometry) at each face."""
    c_co_int = co_root(beta, 3.0 * np.asarray(c_int))
    c_co_p = co_root(beta, 3.0 * np.asarray(c_p))
    return c_co_int - c_co_p


# ---------------------------------------------------------------------------
# Bounded fit (chi, K are both magnitudes; see CaCl2 for why unconstrained
# lm is unsafe for this class of asymmetric-in-chi model)
# ---------------------------------------------------------------------------
def fit_model_bounded(model_fn, c_int, c_p, y, beta0=BETA0, label="fit"):
    def resid(beta):
        return model_fn(beta, c_int, c_p) - y

    result = least_squares(resid, beta0, method="trf", bounds=FIT_BOUNDS)
    beta_hat = result.x
    residuals = result.fun
    sse = float(np.sum(residuals ** 2))
    n = len(y)
    p = len(beta_hat)
    dof = max(n - p, 1)

    J = result.jac
    JTJ = J.T @ J
    try:
        JTJ_inv = np.linalg.inv(JTJ)
    except np.linalg.LinAlgError:
        JTJ_inv = np.linalg.pinv(JTJ)
    sigma2 = sse / dof
    cov = sigma2 * JTJ_inv
    se = np.sqrt(np.diag(cov))
    return step2.FitResult(label=label, beta=beta_hat, se=se, sse=sse, n=n, p=p, cov=cov)


# ---------------------------------------------------------------------------
# Multi-start reliability (method trf/lm; avoid dogbox)
# ---------------------------------------------------------------------------
def multistart_reliability(c_int, c_p, y, sse_global, rng, n_starts=N_MULTISTART):
    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), n_starts)
    log_K = rng.uniform(np.log(K_LO), np.log(K_HI), n_starts)
    chi0s, K0s = np.exp(log_chi), np.exp(log_K)

    n_global = 0
    for chi0, K0 in zip(chi0s, K0s):
        def resid(beta):
            return f_lacl3(beta, c_int, c_p) - y

        try:
            res = least_squares(resid, np.array([chi0, K0]), method="trf",
                                 bounds=([1e-6, 1e-8], [1e5, 1e5]))
            sse_hat = float(np.sum(res.fun ** 2))
            if res.success and sse_hat <= sse_global * (1 + 1e-6) + 1e-8:
                n_global += 1
        except Exception:
            pass
    return 100.0 * n_global / n_starts


# ---------------------------------------------------------------------------
# F-test region / profile likelihood
# ---------------------------------------------------------------------------
def f_test_threshold(sse_min, n, p, alpha):
    from scipy.stats import f as f_dist
    Fcrit = f_dist.ppf(1 - alpha, p, n - p)
    return sse_min * (1 + (p / (n - p)) * Fcrit)


def profile_likelihood(beta_hat, se, sse_min, n, p, c_int, c_p, y, param_idx,
                        alpha=0.05, n_grid=41, span_se=8.0):
    # Floor must never exceed beta_hat + span_se*se (the natural upper end of
    # the span) or the grid direction reverses -- a real bug caught here: a
    # fixed 1e-6 floor is fine for chi/K ~ O(0.01-100), but LaCl3's K can be
    # O(1e-7), where beta_hat+span_se*se itself is below 1e-6.
    hi_natural = beta_hat[param_idx] + span_se * se[param_idx]
    lo_bound = max(min(beta_hat[param_idx] - span_se * se[param_idx], hi_natural * 0.5), 1e-10)
    grid = np.linspace(lo_bound, beta_hat[param_idx] + span_se * se[param_idx], n_grid)
    other_idx = 1 - param_idx
    profile_sse = np.zeros(len(grid))
    for i, val in enumerate(grid):
        def resid_fixed(other_val):
            beta = np.zeros(2)
            beta[param_idx] = val
            beta[other_idx] = other_val[0]
            return f_lacl3(beta, c_int, c_p) - y

        x0 = np.array([max(beta_hat[other_idx], 1e-7)])
        res = least_squares(resid_fixed, x0, method="trf", bounds=([1e-8], [1e5]))
        profile_sse[i] = float(np.sum(res.fun ** 2))

    thresh = f_test_threshold(sse_min, n, p, alpha)
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
# Adsorption-hypothesis diagnostics
# ---------------------------------------------------------------------------
def sign_runs_test(resid):
    signs = np.sign(resid)
    signs = signs[signs != 0]
    n = len(signs)
    n_pos = int(np.sum(signs > 0))
    n_neg = n - n_pos
    runs = 1 + int(np.sum(signs[1:] != signs[:-1]))
    if n_pos == 0 or n_neg == 0:
        return runs, np.nan, np.nan
    mean_runs = 2 * n_pos * n_neg / n + 1
    var_runs = (2 * n_pos * n_neg * (2 * n_pos * n_neg - n)) / (n ** 2 * (n - 1))
    z = (runs - mean_runs) / np.sqrt(var_runs) if var_runs > 0 else np.nan
    return runs, mean_runs, z


def effective_chi_vs_conc(c_int, c_p, y, K_fixed, n_window=N_WINDOW_POINTS,
                           stride=WINDOW_STRIDE, chi0=44.0):
    order = np.argsort(c_int)
    c_int_s, c_p_s, y_s = c_int[order], c_p[order], y[order]
    n = len(c_int_s)

    centers, chi_locals, chi_ses, reliable = [], [], [], []
    start = 0
    while start + n_window <= n:
        sl = slice(start, start + n_window)
        ci, cp, yy = c_int_s[sl], c_p_s[sl], y_s[sl]

        def resid(chi):
            return f_lacl3((chi[0], K_fixed), ci, cp) - yy

        x0_win = min(max(chi0, 1e-3), WINDOW_CHI_BOUND / 2)
        res = least_squares(resid, np.array([x0_win]), method="trf",
                             bounds=([1e-6], [WINDOW_CHI_BOUND]))
        sse = float(np.nansum(res.fun ** 2))
        dof = max(len(yy) - 1, 1)
        sigma2 = sse / dof
        J = res.jac
        JTJ = J.T @ J
        try:
            se = np.sqrt(sigma2 * np.linalg.inv(JTJ)[0, 0])
        except np.linalg.LinAlgError:
            se = np.nan

        chi_hat = float(res.x[0])
        at_bound = chi_hat > 0.995 * WINDOW_CHI_BOUND
        se_ok = np.isfinite(se) and se < WINDOW_RELSE_UNRELIABLE * max(chi_hat, 1e-3)

        centers.append(float(np.mean(ci)))
        chi_locals.append(chi_hat)
        chi_ses.append(float(se))
        reliable.append(bool(se_ok and not at_bound))
        start += stride

    return (np.array(centers), np.array(chi_locals), np.array(chi_ses),
            np.array(reliable))


# ---------------------------------------------------------------------------
# Per-dataset analysis
# ---------------------------------------------------------------------------
def analyze_dataset(rng):
    df_full = load_dataset()
    c_int = df_full["c_int"].to_numpy()
    c_p = df_full["c_p"].to_numpy()
    y = df_full["deltaC_m"].to_numpy()

    fr = fit_model_bounded(f_lacl3, c_int, c_p, y, beta0=BETA0,
                            label=f"{DATASET_LABEL}: LaCl3 3:1 corrected")

    resid = f_lacl3(fr.beta, c_int, c_p) - y
    rmse = float(np.sqrt(np.mean(resid ** 2)))
    y_range = float(y.max() - y.min())
    rmse_over_range = rmse / y_range
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    pseudo_r2 = 1.0 - fr.sse / ss_tot

    pct_global = multistart_reliability(c_int, c_p, y, fr.sse, rng)

    profile_cis = []
    for idx in range(2):
        grid, prof_sse, thresh = profile_likelihood(
            fr.beta, fr.se, fr.sse, fr.n, fr.p, c_int, c_p, y, idx, alpha=ALPHA,
        )
        lo, hi = profile_ci_from_curve(grid, prof_sse, thresh)
        profile_cis.append((lo, hi))

    runs, mean_runs, z_runs = sign_runs_test(resid[np.argsort(c_int)])

    centers, chi_locals, chi_ses, reliable = effective_chi_vs_conc(
        c_int, c_p, y, K_fixed=fr.beta[1], chi0=fr.beta[0],
    )

    return dict(
        df_full=df_full, fr=fr, resid=resid,
        rmse=rmse, rmse_over_range=rmse_over_range, pseudo_r2=pseudo_r2,
        pct_global=pct_global, profile_cis=profile_cis,
        runs=runs, mean_runs=mean_runs, z_runs=z_runs,
        chi_window_centers=centers, chi_window_locals=chi_locals,
        chi_window_ses=chi_ses, chi_window_reliable=reliable,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def print_summary(r):
    fr = r["fr"]
    lo1, hi1 = r["profile_cis"][0]
    lo2, hi2 = r["profile_cis"][1]
    print("=" * 100)
    print(f"SUMMARY: LaCl3 (3:1) Donnan regression, {DATASET_LABEL}")
    print("=" * 100)
    print(f"n = {fr.n}")
    print(f"|chi| = {fr.beta[0]:.4f} mM  (95% CI {lo1:.4f}, {hi1:.4f})")
    print(f"K     = {fr.beta[1]:.4g}   (95% CI {lo2:.4g}, {hi2:.4g})")
    print(f"SSE = {fr.sse:.4f}   RMSE/range = {r['rmse_over_range']:.4f}   "
          f"pseudo-R2 = {r['pseudo_r2']:.4f}")
    print(f"Multi-start: {r['pct_global']:.1f}% global")
    print(f"Residual sign-runs: {r['runs']} (expected {r['mean_runs']:.1f}), "
          f"z = {r['z_runs']:+.2f}")
    rel = r["chi_window_reliable"]
    n_rel = int(rel.sum())
    print(f"Effective-chi windows: {len(rel)} total, {n_rel} reliable")
    if n_rel >= 2:
        locs = r["chi_window_locals"][rel]
        cens = r["chi_window_centers"][rel]
        slope = np.polyfit(cens, locs, 1)[0]
        print(f"  chi range [{locs.min():.2f}, {locs.max():.2f}] mM, "
              f"linear trend slope = {slope:+.4f} mM per mM c_int")
    else:
        print("  too few reliable windows to characterize a trend")
    print()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_fit(r, outfile):
    fr = r["fr"]
    df = r["df_full"]
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    order = np.argsort(c_int)
    y_fit = f_lacl3(fr.beta, c_int[order], c_p[order])

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 0.95))
    ax.plot(c_int, y, "o", color=plot_style.OKABE_ITO["reddish_purple"], markersize=3,
            alpha=0.5, label="data")
    ax.plot(c_int[order], y_fit, "k--", linewidth=2, label=r"$f_{3:1}$ fit")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=1, fontsize=8,
              borderaxespad=0)
    fig.suptitle(f"LaCl3 (3:1) Donnan model fit vs. data\n{DATASET_LABEL}", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_residuals_vs_conc(r, outfile):
    df = r["df_full"]
    c_int = df["c_int"].to_numpy()
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 0.95))
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.plot(c_int, r["resid"], "o", color=plot_style.OKABE_ITO["reddish_purple"],
            markersize=3, alpha=0.5)
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel("residual [mM]")
    fig.suptitle(f"Residuals vs. concentration\n(adsorption-hypothesis diagnostic; {DATASET_LABEL})",
                  fontsize=10)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_effective_chi_vs_conc(r, outfile):
    centers = r["chi_window_centers"]
    locals_ = r["chi_window_locals"]
    ses = r["chi_window_ses"]
    rel = r["chi_window_reliable"]

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.3))
    color = plot_style.OKABE_ITO["reddish_purple"]
    if np.any(rel):
        ax.errorbar(centers[rel], locals_[rel], yerr=ses[rel], fmt="o-", color=color,
                    markersize=4, capsize=2, elinewidth=1, label="reliable window")
    if np.any(~rel):
        ax.plot(centers[~rel], locals_[~rel], "x", color=color, markersize=5, alpha=0.4,
                label="unreliable window")
    ax.axhline(r["fr"].beta[0], color=color, ls=":", lw=1, alpha=0.6)
    ax.set_xlabel(r"window-center $\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"local $\boldsymbol{|\chi|}$ [mM]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1, fontsize=8,
              borderaxespad=0)
    fig.suptitle(f"Effective $|\\chi|$ vs. concentration\n(adsorption drift diagnostic; {DATASET_LABEL})",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_appendix(r, outfile):
    fr = r["fr"]
    df = r["df_full"]
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    order = np.argsort(c_int)
    y_fit = f_lacl3(fr.beta, c_int[order], c_p[order])
    color = plot_style.OKABE_ITO["reddish_purple"]

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 4, figsize=(w, w * 0.30))

    ax = axes[0]
    ax.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4, label="data")
    ax.plot(c_int[order], y_fit, "k--", linewidth=1.5, label=r"$f_{3:1}$ fit")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]", fontsize=8)
    ax.set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]", fontsize=8)
    ax.set_title("(a) fit vs. data", fontsize=8)
    ax.legend(fontsize=6, loc="best")
    ax.tick_params(labelsize=7)

    ax = axes[1]
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.plot(c_int, r["resid"], "o", color=color, markersize=3, alpha=0.5)
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]", fontsize=8)
    ax.set_ylabel("residual [mM]", fontsize=8)
    ax.set_title(f"(b) residual ($\\boldsymbol{{z}}={r['z_runs']:+.1f}$)", fontsize=8)
    ax.tick_params(labelsize=7)

    ax = axes[2]
    centers = r["chi_window_centers"]
    locals_ = r["chi_window_locals"]
    ses = r["chi_window_ses"]
    rel = r["chi_window_reliable"]
    if np.any(rel):
        ax.errorbar(centers[rel], locals_[rel], yerr=ses[rel], fmt="o-", color=color,
                    markersize=3, capsize=2, elinewidth=1)
    if np.any(~rel):
        ax.plot(centers[~rel], locals_[~rel], "x", color=color, markersize=4, alpha=0.35)
    ax.axhline(fr.beta[0], color=color, ls=":", lw=1, alpha=0.6)
    ax.set_xlabel(r"window-center $\boldsymbol{c_{int}}$ [mM]", fontsize=8)
    ax.set_ylabel(r"local $\boldsymbol{|\chi|}$ [mM]", fontsize=8)
    ax.set_title(r"(c) effective $\boldsymbol{|\chi|(c_{int})}$", fontsize=8)
    ax.tick_params(labelsize=7)
    if np.any(rel):
        ax.set_ylim(-0.05 * locals_[rel].max(), 1.15 * locals_[rel].max())

    ax = axes[3]
    ax.axis("off")
    lo1, hi1 = r["profile_cis"][0]
    lo2, hi2 = r["profile_cis"][1]
    info = (f"$n$ = {fr.n}\n\n"
            f"$|\\chi|$ = {fr.beta[0]:.4g} mM\n"
            f"95% CI ({lo1:.3g}, {hi1:.3g})\n\n"
            f"$K$ = {fr.beta[1]:.4g}\n"
            f"95% CI ({lo2:.4g}, {hi2:.4g})\n\n"
            f"SSE = {fr.sse:.4g}\n"
            f"pseudo-$R^2$ = {r['pseudo_r2']:.4f}\n"
            f"RMSE/range = {r['rmse_over_range']:.4f}\n"
            f"multi-start: {r['pct_global']:.0f}% global")
    ax.text(0.05, 0.95, info, transform=ax.transAxes, fontsize=8, va="top",
            ha="left", family="monospace")

    fig.suptitle(DATASET_LABEL, fontsize=11)
    fig.tight_layout(w_pad=2.0)
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 100)
    print("LaCl3 (3:1) Donnan regression: adsorption-hypothesis diagnostics")
    print("=" * 100)
    print(f"APPLY_CP = {APPLY_CP}")

    rng = np.random.default_rng(RNG_SEED)
    r = analyze_dataset(rng)
    print_summary(r)

    plot_fit(r, FIG_DIR / "lacl3_fit.png")
    plot_residuals_vs_conc(r, FIG_DIR / "lacl3_residuals_vs_conc.png")
    plot_effective_chi_vs_conc(r, FIG_DIR / "lacl3_effective_chi_vs_conc.png")
    plot_appendix(r, FIG_DIR / "lacl3_appendix.png")

    return r


if __name__ == "__main__":
    main()
