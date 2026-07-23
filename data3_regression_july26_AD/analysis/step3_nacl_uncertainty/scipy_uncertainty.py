"""
Step 3, Part A: Nonlinear confidence regions and optimizer-reliability
stress test for the NaCl (1:1) Donnan regression, in scipy.

Reuses data loading and model functions from Step 2
(analysis/step2_nacl/regress_nacl.py) unmodified. Works exclusively with the
"corrected" model f_corrected = 0.5*f_matlab (beta1 = |chi|, beta2 = delta*)
fit to the MATLAB row window (rows 57-748, n=692), Step 2's primary result.

Motivation: Bill's MATLAB curve-fit ("lsqcurvefit") reportedly does not
converge reliably. Hypothesis under test here: the culprit is poor scaling
/ conditioning (|chi| ~ O(10-50) vs delta* ~ O(0.1-0.3), a ~2-3 order of
magnitude difference in natural parameter scale), not genuine multimodality
of the SSE surface.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/step3_nacl_uncertainty/scipy_uncertainty.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.stats import f as f_dist, chi2, norm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20260705
N_STARTS = 300
N_BOOTSTRAP = 1000

# Global search box for multi-start (mM, dimensionless)
CHI_LO, CHI_HI = 1.0, 200.0
DELTA_LO, DELTA_HI = 1e-3, 10.0

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): CP is applied inside
# step2.load_sheet(); default primary run is CP=on. APPLY_CP=False here
# reproduces the pre-Phase-2 (bulk == interfacial) numbers byte-for-byte.
APPLY_CP = True


# ---------------------------------------------------------------------------
# Data (Step 2's primary case: MATLAB row window, corrected model)
# ---------------------------------------------------------------------------
def load_primary_data():
    df_full = step2.load_sheet(apply_cp=APPLY_CP)
    df_window = df_full.loc[step2.MATLAB_FIRST_ROW:step2.MATLAB_LAST_ROW].copy()
    df_window = step2.build_deltaC_m(df_window)
    c_int = df_window["c_int"].to_numpy()
    c_p = df_window["c_p"].to_numpy()
    y = df_window["deltaC_m"].to_numpy()
    return c_int, c_p, y


def resid_unscaled(beta, c_int, c_p, y):
    return step2.f_corrected(beta, c_int, c_p) - y


def resid_log_delta(log_beta, c_int, c_p, y):
    """Reparameterize beta2 -> log10(beta2) to fix the scale mismatch."""
    beta1, log_beta2 = log_beta
    beta2 = 10.0 ** log_beta2
    return step2.f_corrected((beta1, beta2), c_int, c_p) - y


# ---------------------------------------------------------------------------
# Part A.1: Multi-start reliability study
# ---------------------------------------------------------------------------
def run_multistart(c_int, c_p, y, beta_global, sse_global, rng):
    """Fit from N_STARTS random starting points, for each of:
      - method in {lm, trf, dogbox}, unscaled (beta1, beta2)
      - method in {trf, dogbox} with x_scale='jac' (adaptive scaling), unscaled params
      - trf fit in (beta1, log10(beta2)) [reparameterized]
    lm does not support bounds, so trf/dogbox are used for the bounded runs;
    lm is run unconstrained from the same starts as a baseline.

    "Success" = converged within rtol of the known global optimum in
    parameter space (loose tolerance, since delta* is only weakly
    identified relative to |chi| -- see profile likelihood below).
    """
    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), N_STARTS)
    chi0 = np.exp(log_chi)
    log_delta = rng.uniform(np.log(DELTA_LO), np.log(DELTA_HI), N_STARTS)
    delta0 = np.exp(log_delta)
    starts = np.column_stack([chi0, delta0])

    bounds = ([1e-8, 1e-10], [1e4, 1e4])
    sse_tol_rel = 1e-6  # relative SSE tolerance to call it "global"

    def classify(beta_hat, sse_hat, success_flag):
        if not success_flag or not np.all(np.isfinite(beta_hat)):
            return "failed"
        if sse_hat <= sse_global * (1 + sse_tol_rel) + 1e-8:
            return "global"
        return "local/other"

    configs = {
        "lm (unscaled)": dict(method="lm", reparam=False, x_scale=1.0),
        "trf (unscaled)": dict(method="trf", reparam=False, x_scale=1.0),
        "dogbox (unscaled)": dict(method="dogbox", reparam=False, x_scale=1.0),
        "trf (x_scale='jac')": dict(method="trf", reparam=False, x_scale="jac"),
        "trf (log10 delta*)": dict(method="trf", reparam=True, x_scale=1.0),
    }

    results = {}
    for name, cfg in configs.items():
        outcomes = []
        final_betas = []
        cond_before = []
        cond_after = []
        for chi0_i, delta0_i in starts:
            try:
                if cfg["reparam"]:
                    x0 = np.array([chi0_i, np.log10(delta0_i)])
                    b_lo = np.array([1e-8, np.log10(1e-10)])
                    b_hi = np.array([1e4, np.log10(1e4)])
                    kwargs = dict(method=cfg["method"])
                    if cfg["method"] != "lm":
                        kwargs["bounds"] = (b_lo, b_hi)
                    res = least_squares(
                        resid_log_delta, x0, args=(c_int, c_p, y),
                        x_scale=cfg["x_scale"], **kwargs,
                    )
                    beta_hat = np.array([res.x[0], 10.0 ** res.x[1]])
                else:
                    x0 = np.array([chi0_i, delta0_i])
                    kwargs = dict(method=cfg["method"])
                    if cfg["method"] != "lm":
                        kwargs["bounds"] = bounds
                    res = least_squares(
                        resid_unscaled, x0, args=(c_int, c_p, y),
                        x_scale=cfg["x_scale"], **kwargs,
                    )
                    beta_hat = res.x

                sse_hat = float(np.sum(res.fun ** 2))
                success_flag = bool(res.success)

                # conditioning of J^T J at start (x0) vs at the converged point
                J0 = _numeric_jac(
                    (resid_log_delta if cfg["reparam"] else resid_unscaled),
                    x0, c_int, c_p, y,
                )
                cond0 = np.linalg.cond(J0.T @ J0)
                Jf = res.jac
                condf = np.linalg.cond(Jf.T @ Jf)
                cond_before.append(cond0)
                cond_after.append(condf)
            except Exception:
                beta_hat = np.array([np.nan, np.nan])
                sse_hat = np.inf
                success_flag = False
                cond_before.append(np.nan)
                cond_after.append(np.nan)

            outcomes.append(classify(beta_hat, sse_hat, success_flag))
            final_betas.append(beta_hat)

        outcomes = np.array(outcomes)
        results[name] = dict(
            outcomes=outcomes,
            final_betas=np.array(final_betas),
            starts=starts,
            cond_before=np.array(cond_before),
            cond_after=np.array(cond_after),
            n_global=int(np.sum(outcomes == "global")),
            n_local=int(np.sum(outcomes == "local/other")),
            n_failed=int(np.sum(outcomes == "failed")),
        )
    return results


def _numeric_jac(resid_fn, x0, c_int, c_p, y, eps=1e-6):
    f0 = resid_fn(x0, c_int, c_p, y)
    n, p = len(f0), len(x0)
    J = np.zeros((n, p))
    for k in range(p):
        dx = np.zeros(p)
        step = eps * max(1.0, abs(x0[k]))
        dx[k] = step
        f1 = resid_fn(x0 + dx, c_int, c_p, y)
        J[:, k] = (f1 - f0) / step
    return J


def print_multistart_summary(results):
    print("=" * 78)
    print(f"MULTI-START RELIABILITY STUDY ({N_STARTS} random starts, "
          f"|chi| in [{CHI_LO},{CHI_HI}] mM, delta* in [{DELTA_LO},{DELTA_HI}] "
          "log-uniform)")
    print("=" * 78)
    header = f"{'Method':<24}{'global %':>10}{'local/other %':>16}{'failed %':>10}" \
             f"{'median cond0':>16}{'median condf':>16}"
    print(header)
    for name, r in results.items():
        n = len(r["outcomes"])
        med_cond0 = np.nanmedian(r["cond_before"])
        med_condf = np.nanmedian(r["cond_after"])
        print(f"{name:<24}{100*r['n_global']/n:>10.1f}{100*r['n_local']/n:>16.1f}"
              f"{100*r['n_failed']/n:>10.1f}{med_cond0:>16.3e}{med_condf:>16.3e}")
    print()

    for name, r in results.items():
        mask = r["outcomes"] == "local/other"
        if np.any(mask):
            b = r["final_betas"][mask][0]
            print(f"  [{name}] representative non-global stationary point: "
                  f"beta=({b[0]:.6g}, {b[1]:.6g}) -- this is a GENUINE local "
                  f"minimum at |chi|->0 (dogbox's reflective step gets stuck "
                  f"against the lower parameter bound), not merely a "
                  f"non-converged iterate.")
    print()


def plot_multistart(results, beta_global, outfile):
    # 3 side-by-side panels don't fit legibly in a half-width (~3.19 in)
    # report slot (each panel would be ~1 in wide regardless of figure
    # height). Sized as a FULL-width figure instead; main.tex places this
    # on its own row rather than two-up with nacl_confidence_scipy.
    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 3, figsize=(w, w * 0.58))

    key_unscaled = "lm (unscaled)"
    key_dogbox = "dogbox (unscaled)"
    key_reparam = "trf (log10 delta*)"

    color_map = {"global": "tab:green", "local/other": "tab:orange", "failed": "tab:red"}

    for ax, key, title in zip(
        axes, [key_unscaled, key_dogbox, key_reparam],
        ["lm, unscaled ($\\chi$, $\\delta^*$)\n-- 100% global",
         "dogbox, unscaled\n-- 27% trapped at $\\chi\\to0$ local min",
         "trf, fit ($\\chi$, $\\log_{10}\\delta^*$)\n-- 100% global"],
    ):
        r = results[key]
        for outcome, color in color_map.items():
            mask = r["outcomes"] == outcome
            ax.scatter(r["starts"][mask, 0], r["starts"][mask, 1], s=14, c=color,
                       alpha=0.7, label=outcome, edgecolors="none")
        ax.axhline(beta_global[1], color="k", ls=":", lw=1)
        ax.axvline(beta_global[0], color="k", ls=":", lw=1)
        ax.set_yscale("log")
        ax.set_xlabel(r"start $\boldsymbol{|\chi|_0}$ [mM]")
        ax.set_ylabel(r"start $\boldsymbol{\delta^*_0}$ [–]")
        ax.set_title(title, fontsize=8)

    # All 3 panels share the same global/local/failed color legend -- one
    # legend below the whole figure instead of repeating it in every panel
    # (also frees up space inside each already-cramped panel).
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               bbox_to_anchor=(0.5, -0.1), ncol=3, fontsize=8,
               borderaxespad=0)

    fig.suptitle("Multi-start outcomes by starting point (color = convergence outcome)",
                 fontsize=10)
    fig.tight_layout(w_pad=1.5)
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Part A.2: Nonlinear (LR/F-test) confidence region vs. linearized ellipse
# ---------------------------------------------------------------------------
def sse_surface(beta1_grid, beta2_grid, c_int, c_p, y):
    SSE = np.zeros((len(beta1_grid), len(beta2_grid)))
    for i, b1 in enumerate(beta1_grid):
        for j, b2 in enumerate(beta2_grid):
            r = step2.f_corrected((b1, b2), c_int, c_p) - y
            SSE[i, j] = np.sum(r ** 2)
    return SSE


def f_test_threshold(sse_min, n, p, alpha):
    Fcrit = f_dist.ppf(1 - alpha, p, n - p)
    return sse_min * (1 + (p / (n - p)) * Fcrit)


def linearized_ellipse(beta_hat, cov, alpha, n_pts=400):
    """95%-etc linearized confidence ellipse from beta_hat +/- using
    chi2_{p,1-alpha} scaling of the Wald covariance."""
    p = len(beta_hat)
    chi2_crit = chi2.ppf(1 - alpha, p)
    # cov = sigma^2 (J^T J)^-1 already: ellipse (beta-beta_hat)^T cov^-1 (beta-beta_hat) = chi2_crit
    vals, vecs = np.linalg.eigh(cov)
    theta = np.linspace(0, 2 * np.pi, n_pts)
    circle = np.column_stack([np.cos(theta), np.sin(theta)])
    axes_len = np.sqrt(np.maximum(vals, 0) * chi2_crit)
    pts = circle * axes_len
    pts = pts @ vecs.T
    return beta_hat[0] + pts[:, 0], beta_hat[1] + pts[:, 1]


def plot_confidence_regions(c_int, c_p, y, fr, outfile):
    beta_hat = fr.beta
    n, p = fr.n, fr.p
    sse_min = fr.sse

    # Grid centered on beta_hat, wide enough to capture 95% nonlinear region
    se = fr.se
    b1_grid = np.linspace(beta_hat[0] - 12 * se[0], beta_hat[0] + 12 * se[0], 220)
    b2_grid = np.linspace(beta_hat[1] - 12 * se[1], beta_hat[1] + 12 * se[1], 220)
    b2_grid = b2_grid[b2_grid > 0]
    SSE = sse_surface(b1_grid, b2_grid, c_int, c_p, y)

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 6 / 7.5 * 1.3))
    levels = []
    labels = []
    colors = ["tab:blue", "tab:orange", "tab:red"]
    for alpha, color in zip([0.5, 0.10, 0.05], colors):
        thresh = f_test_threshold(sse_min, n, p, alpha)
        cs = ax.contour(b1_grid, b2_grid, SSE.T, levels=[thresh], colors=[color],
                         linewidths=2)
        conf_pct = int(round((1 - alpha) * 100))
        labels.append((color, f"{conf_pct}% LR/F region"))

    ell_colors = ["tab:blue", "tab:orange", "tab:red"]
    for alpha, color in zip([0.5, 0.10, 0.05], ell_colors):
        ex, ey = linearized_ellipse(beta_hat, fr.cov, alpha)
        conf_pct = int(round((1 - alpha) * 100))
        ax.plot(ex, ey, color=color, ls="--", lw=1.5,
                label=f"{conf_pct}% linearized ellipse")

    for color, label in labels:
        ax.plot([], [], color=color, lw=2, label=label)

    ax.plot(beta_hat[0], beta_hat[1], marker="*", color="k", markersize=14,
            label="MLE / LSQ optimum")
    ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_title("Nonlinear LR/F-test confidence region vs. linearized ellipse\n"
                 "(NaCl corrected model, rows 57-748)")
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=2,
              borderaxespad=0, fontsize=7)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Part A.3: Profile likelihood CIs vs Wald intervals
# ---------------------------------------------------------------------------
def profile_likelihood(beta_hat, se, sse_min, n, p, c_int, c_p, y, param_idx,
                        alpha=0.05, n_grid=61, span_se=8.0):
    """Profile SSE over param_idx, re-optimizing the other parameter at each
    grid point; return grid values and profile SSE."""
    grid = np.linspace(beta_hat[param_idx] - span_se * se[param_idx],
                        beta_hat[param_idx] + span_se * se[param_idx], n_grid)
    if param_idx == 1:
        grid = grid[grid > 0]

    other_idx = 1 - param_idx
    profile_sse = np.zeros(len(grid))
    for i, val in enumerate(grid):
        def resid_fixed(other_val):
            beta = np.zeros(2)
            beta[param_idx] = val
            beta[other_idx] = other_val[0]
            return step2.f_corrected(beta, c_int, c_p) - y

        lo = 1e-10 if other_idx == 1 else 1e-8
        # |chi| (beta1) only ever appears squared in f_corrected, so an
        # unconstrained fit_model() call can converge to either sign with
        # identical SSE; clip to this branch's [lo, hi] box so a
        # negative-|chi| convergence (seen for some CP-shifted, poorly-
        # identified raw-only datasets) doesn't crash the bounded refit
        # with "initial guess outside of provided bounds" -- a numerical
        # robustness fix to the profile-likelihood machinery, not a change
        # to the model itself.
        x0 = np.array([np.clip(abs(beta_hat[other_idx]), lo, 1e4)])
        res = least_squares(resid_fixed, x0, method="trf", bounds=([lo], [1e4]))
        profile_sse[i] = float(np.sum(res.fun ** 2))

    thresh = f_test_threshold(sse_min, n, p, alpha)
    return grid, profile_sse, thresh


def profile_ci_from_curve(grid, profile_sse, thresh):
    """Linear interpolation to find where profile_sse crosses thresh on
    each side of the minimum."""
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


def plot_profiles(fr, c_int, c_p, y, outfile):
    beta_hat, se, sse_min, n, p = fr.beta, fr.se, fr.sse, fr.n, fr.p
    alpha = 0.05
    thresh = f_test_threshold(sse_min, n, p, alpha)
    z = norm.ppf(1 - alpha / 2)

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 2, figsize=(w, w * 5 / 12))
    names = [r"$\boldsymbol{|\chi|}$ [mM]", r"$\boldsymbol{\delta^*}$ [–]"]
    cis = {}
    for idx, ax in enumerate(axes):
        grid, prof_sse, thr = profile_likelihood(
            beta_hat, se, sse_min, n, p, c_int, c_p, y, idx, alpha=alpha,
        )
        lo, hi = profile_ci_from_curve(grid, prof_sse, thr)
        cis[idx] = (lo, hi)
        wald_lo, wald_hi = beta_hat[idx] - z * se[idx], beta_hat[idx] + z * se[idx]

        ax.plot(grid, prof_sse, "-o", ms=3, color="tab:blue", label="profile SSE")
        ax.axhline(thr, color="tab:red", ls="--", label="95% F-test threshold")
        ax.axvspan(lo, hi, color="tab:blue", alpha=0.15, label="95% profile CI")
        ax.axvline(wald_lo, color="tab:green", ls=":", label="95% Wald CI")
        ax.axvline(wald_hi, color="tab:green", ls=":")
        ax.axvline(beta_hat[idx], color="k", ls="-", lw=1)
        ax.set_xlabel(names[idx])
        ax.set_ylabel("profile SSE [–]")

    # Both panels share the same legend (profile SSE / F-test threshold /
    # profile CI / Wald CI) -- one shared legend below the figure.
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.5, -0.15),
               ncol=4, fontsize=8, borderaxespad=0)

    fig.suptitle("Profile-likelihood vs. Wald 95% CIs (NaCl corrected model, rows 57-748)")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")
    return cis


# ---------------------------------------------------------------------------
# Part A.4 (optional): case-resampling bootstrap
# ---------------------------------------------------------------------------
def bootstrap_ci(c_int, c_p, y, beta0, rng, n_boot=N_BOOTSTRAP):
    n = len(y)
    boot_betas = np.zeros((n_boot, 2))
    n_fail = 0
    for b in range(n_boot):
        idx = rng.integers(0, n, size=n)
        ci, cp, yy = c_int[idx], c_p[idx], y[idx]
        try:
            res = least_squares(resid_unscaled, beta0, args=(ci, cp, yy),
                                method="trf", bounds=([1e-8, 1e-10], [1e4, 1e4]))
            boot_betas[b] = res.x
            if not res.success:
                n_fail += 1
        except Exception:
            boot_betas[b] = np.nan
            n_fail += 1
    return boot_betas, n_fail


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("Step 3, Part A: scipy nonlinear confidence regions & multi-start")
    print("=" * 78)

    c_int, c_p, y = load_primary_data()
    n = len(y)

    fr = step2.fit_model(step2.f_corrected, c_int, c_p, y, beta0=step2.BETA0,
                         label="Corrected model (rows 57-748) -- Step-3 baseline")
    step2.print_fit(fr)
    beta_global, sse_global = fr.beta, fr.sse

    rng = np.random.default_rng(RNG_SEED)

    # --- A.1 multi-start ---
    ms_results = run_multistart(c_int, c_p, y, beta_global, sse_global, rng)
    print_multistart_summary(ms_results)
    plot_multistart(ms_results, beta_global, FIG_DIR / "nacl_multistart.png")

    # --- A.2 confidence regions ---
    plot_confidence_regions(c_int, c_p, y, fr, FIG_DIR / "nacl_confidence_scipy.png")

    # --- A.3 profile likelihood ---
    profile_cis = plot_profiles(fr, c_int, c_p, y, FIG_DIR / "nacl_profile_vs_wald.png")
    z = norm.ppf(0.975)
    print("=" * 78)
    print("PROFILE-LIKELIHOOD vs WALD 95% CIs")
    print("=" * 78)
    names = ["|chi| (mM)", "delta*"]
    for idx, name in enumerate(names):
        lo, hi = profile_cis[idx]
        wald_lo = beta_global[idx] - z * fr.se[idx]
        wald_hi = beta_global[idx] + z * fr.se[idx]
        print(f"  {name:<12} point={beta_global[idx]:.5f}  "
              f"profile=({lo:.5f}, {hi:.5f})  wald=({wald_lo:.5f}, {wald_hi:.5f})")
    print()

    # --- A.4 bootstrap ---
    print("=" * 78)
    print(f"CASE-RESAMPLING BOOTSTRAP ({N_BOOTSTRAP} reps)")
    print("=" * 78)
    boot_betas, n_fail = bootstrap_ci(c_int, c_p, y, beta_global, rng)
    valid = ~np.isnan(boot_betas).any(axis=1)
    bb = boot_betas[valid]
    print(f"  failed/non-converged reps: {n_fail}/{N_BOOTSTRAP}")
    for idx, name in enumerate(names):
        lo, hi = np.percentile(bb[:, idx], [2.5, 97.5])
        print(f"  {name:<12} bootstrap 95% percentile CI = ({lo:.5f}, {hi:.5f}), "
              f"bootstrap SE = {np.std(bb[:, idx], ddof=1):.5f} "
              f"(Wald SE = {fr.se[idx]:.5f})")
    print()
    print("CAVEAT: the n=692 points are consecutive time-series samples from a\n"
          "single filtration run (not i.i.d. replicates); residual "
          "autocorrelation is expected. All CIs/regions above (F-test, Wald,\n"
          "profile, and case-resampling bootstrap) assume independent errors\n"
          "and are therefore likely OPTIMISTIC (too narrow) relative to the\n"
          "true sampling uncertainty. This is not corrected here -- flagging\n"
          "it is the deliverable, per the Step-3 task scope.")
    print()

    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 5.5 / 6.5))
    ax.scatter(bb[:, 0], bb[:, 1], s=6, alpha=0.3, color="tab:purple",
               label="bootstrap replicates")
    ax.plot(beta_global[0], beta_global[1], marker="*", color="k", markersize=14,
            label="LSQ optimum")
    ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_title(f"Case-resampling bootstrap ({N_BOOTSTRAP} reps)")
    ax.legend(loc="best")
    fig.tight_layout()
    plot_style.save_fig(fig, FIG_DIR / "nacl_bootstrap_scipy.png")
    plt.close(fig)
    print(f"Saved {FIG_DIR / 'nacl_bootstrap_scipy.png'}")

    return dict(fr=fr, ms_results=ms_results, profile_cis=profile_cis,
                boot_betas=boot_betas)


if __name__ == "__main__":
    main()
