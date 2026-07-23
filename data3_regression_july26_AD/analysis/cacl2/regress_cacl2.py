"""
CaCl2 (2:1) Donnan regression + adsorption-hypothesis diagnostics (report Sec. 3).

Model
-----
Working hypothesis under test: Ca2+ adsorbs onto the membrane and changes its
effective fixed charge |chi| in a concentration-dependent way. The current
model treats |chi| as a single constant per dataset (no adsorption), so we
expect CaCl2 fits to be worse than NaCl, with the diagnostic signature being
systematic residual structure vs. concentration and/or an apparent |chi| that
drifts with concentration.

Co-ion = Cl-, solution concentration set by stoichiometry: c_co_s = 2*c_CaCl2
(2 Cl- per CaCl2), used for both c_int and c_p.

Donnan co-ion cubic (report Eq. 21co), lumped partition constant K (>= 0)
absorbing delta_star * delta_co_circ and constant stoichiometric factors:
    c_co_m^3 + |chi|*c_co_m^2 - K*c_co_s^3 = 0
Exactly one positive real root (Descartes: sign pattern +,+,0,-), selected via
numpy.roots([1, |chi|, 0, -K*c_co_s**3]) -> the unique real, positive root.

Transformed response (Cl- flux form, analogous to corrected NaCl):
    deltaC_m = J_s_Cl * l / D_m
    J_s_Cl = 2 * c_p * J_w         (Cl- flux = 2x the CaCl2 molar flux)
    l = 80e-9 m
    D_m = 3*D_Ca_m*D_Cl_m / (2*D_Ca_m + D_Cl_m)   (2:1 ambipolar membrane
                                                     diffusivity)
    D_Ca_m = 0.79e-12 m^2/s, D_Cl_m = 2.03e-12 m^2/s

Prediction = c_co_m(c_int) - c_co_m(c_p) (solve the cubic at each face). Fit
(|chi|, K) to deltaC_m by scipy.optimize.least_squares, bounded to chi,K >= 0
(both are magnitudes by construction -- unlike the NaCl model, where beta1
only ever appears squared inside a sqrt so its sign is irrelevant, the CaCl2
cubic's chi enters unsquared/asymmetrically, so an unconstrained fit can drift
to a spurious negative "chi" with a marginally lower SSE that has no physical
meaning; see analysis/cacl2/README.md).

Reuses analysis/step2_nacl/regress_nacl.py's plot_style and FitResult
dataclass; a bounded fit_model_bounded() is used in place of step2.fit_model
(which is hardcoded to unconstrained method="lm") for exactly this reason.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/cacl2/regress_cacl2.py
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.optimize import least_squares
from scipy.stats import f as f_dist, norm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import cp_authors_bulk  # noqa: E402

plot_style.apply_style()

DATA_FILE = REPO_ROOT / "data" / "BoE Analysis.xlsx"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): all 3 CaCl2 sheets are
# the authors' non-CP bulk convention (Phase 1 validation), so CP is applied
# on top of column B for every dataset here. Default primary run is CP=on;
# APPLY_CP=False reproduces the pre-Phase-2 numbers byte-for-byte.
APPLY_CP = True

RNG_SEED = 20260705
N_MULTISTART = 100
ALPHA = 0.05  # 95% CIs/regions throughout

# Physical constants (membrane-phase diffusivities, m^2/s; report Sec. 3)
D_CA_M = 0.79e-12
D_CL_M = 2.03e-12
L_MEMBRANE = 80e-9  # m
D_CACL2_M = 3.0 * D_CA_M * D_CL_M / (2.0 * D_CA_M + D_CL_M)

BETA0 = np.array([44.0, 0.3])  # initial guess: (|chi| [mM], K [-])
FIT_BOUNDS = ([0.0, 0.0], [np.inf, np.inf])  # |chi|, K are both magnitudes

# Global search box for multi-start (mM, dimensionless)
CHI_LO, CHI_HI = 1.0, 500.0
K_LO, K_HI = 1e-4, 50.0


def fit_model_bounded(model_fn, c_int, c_p, y, beta0=BETA0, label="fit"):
    """Bounded analog of step2.fit_model (which is hardcoded to
    unconstrained method="lm"): trust-region-reflective, chi,K >= 0, since
    both are magnitudes here (see module docstring for why this matters for
    the CaCl2 cubic, unlike the NaCl model)."""
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

# Sliding-window effective-|chi| diagnostic
N_WINDOW_POINTS = 80  # points per window
WINDOW_STRIDE = 40    # points between window starts

DATASETS = [
    {"key": "MC3-0711S", "sheet": "NF270_MC3 07.11.24_SCaCl2",
     "label": "MC3 (07.11.24, S)"},
    {"key": "MC5-0527", "sheet": "NF270_MC5 05.27.26_CaCl2",
     "label": "MC5 (05.27.26)"},
    {"key": "MC5-0627", "sheet": "NF270_MC5 06.27.26_CaCl2",
     "label": "MC5 (06.27.26)"},
]
DATASET_COLORS = plot_style.get_dataset_colors([d["label"] for d in DATASETS])


# ---------------------------------------------------------------------------
# Data loading (same 15-column schema/positions as NaCl: B=c_int, H=c_p, J=Jw)
# ---------------------------------------------------------------------------
def load_sheet(sheet_name: str, apply_cp: bool | None = None) -> pd.DataFrame:
    if apply_cp is None:
        apply_cp = APPLY_CP
    df = pd.read_excel(DATA_FILE, sheet_name=sheet_name, header=0)
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
        # Column B is the authors' non-CP bulk convention for all 3 CaCl2
        # sheets (Phase 1 validation) -- apply CP on top of it.
        out["c_int"] = cp_authors_bulk.apply_cp(out["c_int"], out["c_p"], out["J_w"], "CaCl2")
    return out


def build_deltaC_m(df: pd.DataFrame) -> pd.DataFrame:
    """deltaC_m [mM] = J_s_Cl * l / D_CaCl2_m, J_s_Cl = 2 * c_p * J_w."""
    df = df.copy()
    J_s_cl = 2.0 * df["c_p"] * df["J_w"]  # mol/m^2/s (Cl- flux)
    df["J_s_cl"] = J_s_cl
    df["deltaC_m"] = J_s_cl * L_MEMBRANE / D_CACL2_M  # mM
    return df


def load_dataset(ds: dict):
    df_full = load_sheet(ds["sheet"])
    valid = df_full.notna().all(axis=1) & (df_full["c_int"] > df_full["c_p"]) \
        & (df_full["c_int"] > 0)
    df_valid = df_full[valid].copy()
    df_valid = build_deltaC_m(df_valid)
    return df_valid


# ---------------------------------------------------------------------------
# Model: 2:1 co-ion cubic, single positive real root
# ---------------------------------------------------------------------------
def co_root(beta, c_s):
    """Unique positive real root of c_m^3 + chi*c_m^2 - K*c_s^3 = 0
    (Descartes: sign pattern +,+,0,- guarantees exactly one), vectorized over
    c_s (verified against numpy.roots to <1e-9 relative error over a
    5000-point random scan spanning chi in [1e-3,1e5], K in [1e-8,1e3],
    c_s in [1e-2,1e3]; ~1000x faster than a per-element numpy.roots loop for
    the array sizes here, n up to ~1900).

    Two regimes, selected by m = K*c_s^3/chi^3 (how small the true root t is
    relative to chi):
    - m small (t << chi, cubic term negligible): the general closed-form
      Cardano solution suffers catastrophic cancellation (t is the tiny
      difference of two O(chi) terms). Instead use the asymptotic balance
      chi*t^2 ~= K*c_s^3 -> t0 = sqrt(K*c_s^3/chi), polished by a few Newton
      steps on g(t) = t^3 + chi*t^2 - K*c_s^3 (well-conditioned: no
      cancellation, all terms already the right scale). This regime is
      reached during profile-likelihood/multi-start excursions to extreme
      chi for ill-conditioned (near chi->0) datasets -- see MC3 (07.11.24, S)
      in the README.
    - m not small: general closed-form (complex) Cardano solution of the
      depressed cubic, selecting the unique root with negligible imaginary
      part and positive real part.
    """
    chi, K = beta
    c_s = np.atleast_1d(np.asarray(c_s, dtype=float))
    chi_arr = np.full_like(c_s, chi, dtype=float)
    rhs = K * c_s ** 3

    # The small-root asymptotic only applies for chi > 0 (it assumes t << chi
    # with both positive); least_squares can probe chi <= 0 transiently via
    # finite-difference Jacobian steps even though the physical optimum has
    # chi > 0, so route those through the general Cardano branch instead.
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
        a = 1.0
        b = np.full_like(cs2, chi, dtype=complex)
        c = np.zeros_like(cs2, dtype=complex)
        d = (-K * cs2 ** 3).astype(complex)

        # Depressed cubic x^3 + p*x + q = 0 via t = x - b/(3a)
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


def f_cacl2(beta, c_int, c_p):
    """Prediction = c_co_m(c_int) - c_co_m(c_p), co-ion solution conc. =
    2*c_CaCl2 (stoichiometry) at each face."""
    c_co_int = co_root(beta, 2.0 * np.asarray(c_int))
    c_co_p = co_root(beta, 2.0 * np.asarray(c_p))
    return c_co_int - c_co_p


# ---------------------------------------------------------------------------
# Multi-start reliability (method lm/trf; avoid dogbox -- trapped at the
# boundary for NaCl, per Step 3's finding)
# ---------------------------------------------------------------------------
def multistart_reliability(c_int, c_p, y, sse_global, rng, n_starts=N_MULTISTART):
    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), n_starts)
    log_K = rng.uniform(np.log(K_LO), np.log(K_HI), n_starts)
    chi0s, K0s = np.exp(log_chi), np.exp(log_K)

    n_global = 0
    for chi0, K0 in zip(chi0s, K0s):
        def resid(beta):
            return f_cacl2(beta, c_int, c_p) - y

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
# F-test region / profile likelihood (cubic-model analogs of
# step3_nacl_uncertainty.scipy_uncertainty, which is hardwired to
# step2.f_corrected -- re-implemented here, parameterized by f_cacl2, rather
# than modifying that shared NaCl module)
# ---------------------------------------------------------------------------
def f_test_threshold(sse_min, n, p, alpha):
    Fcrit = f_dist.ppf(1 - alpha, p, n - p)
    return sse_min * (1 + (p / (n - p)) * Fcrit)


def sse_surface(chi_grid, K_grid, c_int, c_p, y):
    SSE = np.zeros((len(chi_grid), len(K_grid)))
    for i, chi in enumerate(chi_grid):
        for j, K in enumerate(K_grid):
            r = f_cacl2((chi, K), c_int, c_p) - y
            SSE[i, j] = np.nansum(r ** 2)
    return SSE


def profile_likelihood(beta_hat, se, sse_min, n, p, c_int, c_p, y, param_idx,
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
            return f_cacl2(beta, c_int, c_p) - y

        x0 = np.array([max(beta_hat[other_idx], 1e-7)])  # beta_hat can sit
        # exactly at the bounded fit's chi=0 boundary; clip so x0 stays
        # strictly inside the re-optimization bounds below.
        res = least_squares(resid_fixed, x0, method="trf", bounds=([1e-8], [1e5]))
        profile_sse[i] = float(np.nansum(res.fun ** 2))

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
    """Count sign runs in residuals sorted by c_int; fewer runs than expected
    under randomness indicates systematic (non-random) structure."""
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


WINDOW_CHI_BOUND = 300.0  # mM; physically generous (>NaCl's largest estimate,
                          # ~170 mM) upper bound for the LOCAL windowed fit,
                          # tighter than the global fit's search box -- keeps
                          # small-sample (n_window=80), weakly-identified
                          # windows from diverging to nonphysical values that
                          # would otherwise dominate the plot's y-scale.
WINDOW_RELSE_UNRELIABLE = 2.0  # exclude a window from the drift plot/trend
                               # if its Wald SE exceeds 2x the point estimate
                               # (or the point estimate sits at the bound)


def effective_chi_vs_conc(c_int, c_p, y, K_fixed, n_window=N_WINDOW_POINTS,
                           stride=WINDOW_STRIDE, chi0=44.0):
    """Sliding window over c_int (sorted): refit |chi| locally (K fixed at
    the global estimate) and return window-center c_int, local |chi| point
    estimate + Wald SE, plus a `reliable` mask flagging windows whose local
    fit is well-behaved (not pinned at the search-box bound, SE not wildly
    larger than the point estimate) -- unreliable windows are still returned
    but excluded from the drift-trend estimate and de-emphasized in the plot,
    so a handful of noisy small-sample windows can't dominate the story."""
    order = np.argsort(c_int)
    c_int_s, c_p_s, y_s = c_int[order], c_p[order], y[order]
    n = len(c_int_s)

    centers, chi_locals, chi_ses, reliable = [], [], [], []
    start = 0
    while start + n_window <= n:
        sl = slice(start, start + n_window)
        ci, cp, yy = c_int_s[sl], c_p_s[sl], y_s[sl]

        def resid(chi):
            return f_cacl2((chi[0], K_fixed), ci, cp) - yy

        x0_win = min(max(chi0, 1e-3), WINDOW_CHI_BOUND / 2)
        res = least_squares(resid, np.array([x0_win]),
                             method="trf", bounds=([1e-6], [WINDOW_CHI_BOUND]))
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
def analyze_dataset(ds, rng):
    label = ds["label"]
    df_full = load_dataset(ds)
    c_int = df_full["c_int"].to_numpy()
    c_p = df_full["c_p"].to_numpy()
    y = df_full["deltaC_m"].to_numpy()

    fr = fit_model_bounded(f_cacl2, c_int, c_p, y, beta0=BETA0,
                            label=f"{label}: CaCl2 2:1 corrected")

    # Scale-free fit metric: RMSE / range(y)
    resid = f_cacl2(fr.beta, c_int, c_p) - y
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
        ds=ds, df_full=df_full, fr=fr, resid=resid,
        rmse=rmse, rmse_over_range=rmse_over_range, pseudo_r2=pseudo_r2,
        pct_global=pct_global, profile_cis=profile_cis,
        runs=runs, mean_runs=mean_runs, z_runs=z_runs,
        chi_window_centers=centers, chi_window_locals=chi_locals,
        chi_window_ses=chi_ses, chi_window_reliable=reliable,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def print_summary_table(results):
    print("=" * 110)
    print("SUMMARY: CaCl2 (2:1) Donnan regression, adsorption-hypothesis diagnostics")
    print("=" * 110)
    header = (f"{'Dataset':<20}{'n':>6}{'|chi| (mM)':>14}{'95% CI':>20}"
              f"{'K':>10}{'95% CI':>20}{'SSE':>14}{'RMSE/range':>12}"
              f"{'pseudo-R2':>11}{'MS % global':>13}{'runs (z)':>14}")
    print(header)
    for r in results:
        fr = r["fr"]
        lo1, hi1 = r["profile_cis"][0]
        lo2, hi2 = r["profile_cis"][1]
        print(f"{r['ds']['label']:<20}{fr.n:>6}{fr.beta[0]:>14.4f}"
              f"  ({lo1:6.3f}, {hi1:6.3f})"
              f"{fr.beta[1]:>10.5f}"
              f"  ({lo2:6.4f}, {hi2:6.4f})"
              f"{fr.sse:>14.3f}{r['rmse_over_range']:>12.4f}"
              f"{r['pseudo_r2']:>11.4f}{r['pct_global']:>12.1f}%"
              f"  {r['runs']:>3d} ({r['z_runs']:+.2f})")
    print()

    print("-" * 110)
    print("Effective |chi| drift (sliding window, K fixed at global estimate):")
    print("-" * 110)
    for r in results:
        centers = r["chi_window_centers"]
        locals_ = r["chi_window_locals"]
        rel = r["chi_window_reliable"]
        n_rel = int(rel.sum())
        if n_rel < 2:
            print(f"{r['ds']['label']:<20}insufficient reliable windows "
                  f"({n_rel}/{len(locals_)})")
            continue
        # linear trend of local chi vs. window-center c_int, RELIABLE windows only
        slope = np.polyfit(centers[rel], locals_[rel], 1)[0]
        print(f"{r['ds']['label']:<20}n_windows={len(locals_):>3} "
              f"({n_rel} reliable)  "
              f"chi range [{locals_[rel].min():.2f}, {locals_[rel].max():.2f}] mM  "
              f"linear trend slope = {slope:+.4f} mM per mM c_int "
              f"({'DRIFTS with concentration' if abs(slope) > 0.05 else 'roughly flat'})")
    print()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_all_fits(results, outfile):
    w = plot_style.fig_width("full")
    n = len(results)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, w * nrows / ncols * 0.95),
                              sharey=False)
    axes = np.atleast_1d(axes).flatten()
    for ax in axes[n:]:
        ax.axis("off")
    for ax, r in zip(axes, results):
        fr = r["fr"]
        label = r["ds"]["label"]
        color = DATASET_COLORS[label]
        df = r["df_full"]
        c_int = df["c_int"].to_numpy()
        c_p = df["c_p"].to_numpy()
        y = df["deltaC_m"].to_numpy()
        order = np.argsort(c_int)
        y_fit = f_cacl2(fr.beta, c_int[order], c_p[order])

        ax.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4, label="data")
        ax.plot(c_int[order], y_fit, "--", color="k", linewidth=2,
                label=r"$f_{2:1}$ fit")
        ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
        ax.set_title(label, fontsize=10)
    fig.supylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")

    shared_handles = [
        Line2D([], [], marker="o", linestyle="none", color="gray", alpha=0.6,
               label="data"),
        Line2D([], [], color="k", linestyle="--", linewidth=2,
               label=r"$f_{2:1}$ fit"),
    ]
    fig.legend(handles=shared_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9, borderaxespad=0)
    fig.suptitle("CaCl2 (2:1) Donnan model fit vs. data, by dataset")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_residuals_vs_conc(results, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.3))
    for r in results:
        label = r["ds"]["label"]
        color = DATASET_COLORS[label]
        c_int = r["df_full"]["c_int"].to_numpy()
        ax.plot(c_int, r["resid"], "o", color=color, markersize=3, alpha=0.4,
                label=label)
    ax.axhline(0, color="k", lw=1, ls=":")
    ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"residual [mM]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1,
              fontsize=8, borderaxespad=0)
    fig.suptitle("Residuals vs. concentration\n(adsorption-hypothesis diagnostic)",
                  fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_effective_chi_vs_conc(results, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 1.3))
    for r in results:
        label = r["ds"]["label"]
        color = DATASET_COLORS[label]
        centers = r["chi_window_centers"]
        locals_ = r["chi_window_locals"]
        ses = r["chi_window_ses"]
        rel = r["chi_window_reliable"]
        if len(centers) == 0:
            continue
        # Reliable windows: solid errorbar + line, in the legend. Unreliable
        # (SE >> point estimate, or pinned at the window search-box bound --
        # a handful of small-sample windows near the very-low-c_int end are
        # essentially unidentified) are shown as faint markers only, so they
        # don't dominate the axis scale or misread as a real trend.
        if np.any(rel):
            ax.errorbar(centers[rel], locals_[rel], yerr=ses[rel], fmt="o-",
                        color=color, markersize=4, capsize=2, elinewidth=1,
                        label=label)
        if np.any(~rel):
            ax.plot(centers[~rel], locals_[~rel], "x", color=color,
                    markersize=5, alpha=0.35)
        fr = r["fr"]
        ax.axhline(fr.beta[0], color=color, ls=":", lw=1, alpha=0.6)
    if any(np.any(r["chi_window_reliable"]) for r in results):
        rel_vals = np.concatenate([r["chi_window_locals"][r["chi_window_reliable"]]
                                    for r in results
                                    if np.any(r["chi_window_reliable"])])
        ax.set_ylim(-0.05 * rel_vals.max(), 1.15 * rel_vals.max())
    ax.set_xlabel(r"window-center $\boldsymbol{c_{int}}$ [mM]")
    ax.set_ylabel(r"local $\boldsymbol{|\chi|}$ [mM]")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=1,
              fontsize=8, borderaxespad=0)
    fig.suptitle("Effective $\\boldsymbol{|\\chi|}$ vs. concentration\n"
                 "(adsorption drift diagnostic; $\\boldsymbol{\\times}$ = unreliable window)",
                 fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


def plot_all_estimates(results, outfile):
    # Displayed at 0.6\linewidth in main.tex (not the "half"=0.49\linewidth
    # minipage convention used elsewhere) -- save at that same fraction of
    # the full page width so on-page text renders at its designed ~11pt
    # scale (see figure-compliance audit, docs/prompts/figure_compliance_audit.md).
    w = plot_style.fig_width("full") * 0.6
    fig, axes = plt.subplots(2, 1, figsize=(w, w * 6.4 / 4.0), sharex=True)
    labels = [r["ds"]["label"] for r in results]
    x = np.arange(len(labels))
    colors = [DATASET_COLORS[label] for label in labels]

    # Point estimates sit at (or interpolate near) the chi=0 boundary for all
    # three datasets, so the profile CI's lower bound can round to a hair
    # above/below the point estimate -- clip yerr at 0 (a physical boundary
    # estimate has a one-sided, not symmetric, CI in that direction).
    ax = axes[0]
    for xi, r, c in zip(x, results, colors):
        fr = r["fr"]
        lo = max(fr.beta[0] - r["profile_cis"][0][0], 0.0)
        hi = max(r["profile_cis"][0][1] - fr.beta[0], 0.0)
        ax.errorbar(xi, fr.beta[0], yerr=[[lo], [hi]], fmt="o", color=c,
                    capsize=5, markersize=8, elinewidth=2)
    ax.set_ylabel(r"$\boldsymbol{|\chi|}$ [mM]")

    ax = axes[1]
    for xi, r, c in zip(x, results, colors):
        fr = r["fr"]
        lo = max(fr.beta[1] - r["profile_cis"][1][0], 0.0)
        hi = max(r["profile_cis"][1][1] - fr.beta[1], 0.0)
        ax.errorbar(xi, fr.beta[1], yerr=[[lo], [hi]], fmt="o", color=c,
                    capsize=5, markersize=8, elinewidth=2)
    ax.set_ylabel(r"$\boldsymbol{K}$ [–]")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")

    for ax in axes:
        ax.label_outer()

    fig.suptitle("CaCl2 point estimates $\\boldsymbol{\\pm}$ 95% profile-likelihood CI\n"
                 "(2:1 cubic model, full range)")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 110)
    print("CaCl2 (2:1) Donnan regression across all processed datasets")
    print("=" * 110)
    print(f"APPLY_CP = {APPLY_CP}")

    rng = np.random.default_rng(RNG_SEED)
    results = []
    for ds in DATASETS:
        print(f"--- {ds['label']} ({ds['sheet']}) ---")
        t0 = time.time()
        r = analyze_dataset(ds, rng)
        print(f"  done in {time.time()-t0:.1f} s")
        results.append(r)

    print()
    print_summary_table(results)

    plot_all_fits(results, FIG_DIR / "cacl2_fits.png")
    plot_residuals_vs_conc(results, FIG_DIR / "cacl2_residuals_vs_conc.png")
    plot_effective_chi_vs_conc(results, FIG_DIR / "cacl2_effective_chi_vs_conc.png")
    plot_all_estimates(results, FIG_DIR / "cacl2_estimates.png")

    return results


if __name__ == "__main__":
    main()
