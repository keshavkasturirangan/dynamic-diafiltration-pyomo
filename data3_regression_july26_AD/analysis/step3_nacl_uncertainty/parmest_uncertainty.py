"""
Step 3, Part B: Pyomo/ParmEst parameter estimation, confidence regions,
bootstrap, and multi-start for the NaCl (1:1) Donnan regression.

Uses pyomo.contrib.parmest's Experiment-class API (pyomo 6.10.x), NOT the
older callback-based API (`parmest.Estimator(model_function, data, ...)`)
found in older pyomo tutorials -- that interface was removed. See the
"ParmEst API notes" section of README.md in this directory for what changed
and what tripped us up.

Reuses data loading and the (beta1=|chi|, beta2=delta*) corrected Donnan
model from Step 2 (analysis/step2_nacl/regress_nacl.py) and reuses the
scipy nonlinear-region / F-test utilities from Part A
(scipy_uncertainty.py) so the two confidence regions can be overlaid on one
plot for a direct comparison.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/step3_nacl_uncertainty/parmest_uncertainty.py
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

import pyomo.environ as pyo
from pyomo.contrib.parmest.experiment import Experiment
import pyomo.contrib.parmest.parmest as parmest

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step3_nacl_uncertainty"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import scipy_uncertainty as s3  # noqa: E402  (reuse F-test / ellipse utilities)

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20260705
N_MULTISTART = 30       # fewer than scipy's 300: each ParmEst theta_est()
                        # solves a 692-experiment extensive-form NLP (~1.3 s)
N_BOOTSTRAP = 100       # fewer than scipy's 1000, for the same reason
GRID_N = 21             # 21x21 = 441 objective_at_theta evaluations
                        # (~0.45 s each -> ~3-4 min for the confidence region)

CHI_LO, CHI_HI = 1.0, 200.0
DELTA_LO, DELTA_HI = 1e-3, 10.0

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): CP is applied inside
# step2.load_sheet(); default primary run is CP=on. APPLY_CP=False here
# reproduces the pre-Phase-2 (bulk == interfacial) numbers byte-for-byte.
APPLY_CP = True


# ---------------------------------------------------------------------------
# Experiment class (Pyomo 6.10 Experiment-class API)
# ---------------------------------------------------------------------------
class NaClDonnanExperiment(Experiment):
    """One experiment == one (c_int, c_p, deltaC_m) data row.

    beta1 = |chi| [mM], beta2 = delta* [-] are declared as Vars, fixed at an
    initial guess, and listed in `unknown_parameters`; ParmEst's extensive
    form (EF) then ties every experiment's beta1/beta2 to a single shared,
    unfixed decision variable and solves the joint SSE-minimization NLP with
    Ipopt.
    """

    def __init__(self, data_row, theta0=None):
        self.data = data_row
        self.model = None
        self.theta0 = theta0 or {"beta1": 50.0, "beta2": 0.33}

    def create_model(self):
        m = pyo.ConcreteModel()
        m.beta1 = pyo.Var(initialize=self.theta0["beta1"], bounds=(1e-6, 1e4))
        m.beta2 = pyo.Var(initialize=self.theta0["beta2"], bounds=(1e-8, 1e4))
        m.beta1.fix()
        m.beta2.fix()

        m.c_int = pyo.Var(initialize=self.data["c_int"])
        m.c_int.fix()
        m.c_p = pyo.Var(initialize=self.data["c_p"])
        m.c_p.fix()

        m.deltaC_m = pyo.Var(initialize=self.data["deltaC_m"])

        def donnan_rule(m):
            return m.deltaC_m == 0.5 * (
                pyo.sqrt(m.beta1 ** 2 + 4 * m.beta2 * m.c_int ** 2)
                - pyo.sqrt(m.beta1 ** 2 + 4 * m.beta2 * m.c_p ** 2)
            )

        m.donnan_eq = pyo.Constraint(rule=donnan_rule)
        self.model = m

    def label_model(self):
        m = self.model

        m.experiment_outputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_outputs.update([(m.deltaC_m, self.data["deltaC_m"])])

        m.unknown_parameters = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.unknown_parameters.update((k, pyo.value(k)) for k in [m.beta1, m.beta2])

        m.experiment_inputs = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.experiment_inputs.update(
            [(m.c_int, self.data["c_int"]), (m.c_p, self.data["c_p"])]
        )

        # IMPORTANT (API friction, see README): leave measurement_error as
        # None so parmest's cov_est() ESTIMATES the residual variance from
        # SSE/(n-p), matching scipy's sigma_hat^2 = SSE/(n-p). Supplying a
        # placeholder number here (e.g. 1.0) silently changes cov_est()'s
        # covariance scale by a factor of (sigma_hat^2)^-1 -- easy to miss.
        m.measurement_error = pyo.Suffix(direction=pyo.Suffix.LOCAL)
        m.measurement_error.update([(m.deltaC_m, None)])

    def get_labeled_model(self):
        if self.model is None:
            self.create_model()
            self.label_model()
        return self.model


def build_experiment_list(df, theta0=None):
    return [NaClDonnanExperiment(df.iloc[i], theta0=theta0) for i in range(len(df))]


def load_primary_data():
    df_full = step2.load_sheet(apply_cp=APPLY_CP)
    df_window = df_full.loc[step2.MATLAB_FIRST_ROW:step2.MATLAB_LAST_ROW].copy()
    df_window = step2.build_deltaC_m(df_window)
    return df_window


# ---------------------------------------------------------------------------
# B.1: theta_est + cov_est, compare to scipy
# ---------------------------------------------------------------------------
def run_theta_est(df, scipy_fr):
    print("=" * 78)
    print("PARMEST theta_est() vs. scipy least_squares")
    print("=" * 78)

    exp_list = build_experiment_list(df, theta0={"beta1": 50.0, "beta2": 0.33})
    pest = parmest.Estimator(exp_list, obj_function="SSE")

    t0 = time.time()
    obj, theta = pest.theta_est()
    t_theta = time.time() - t0

    n = len(df)
    sse_parmest = obj * n  # theta_est's obj is the SSE averaged over n experiments
    cov = pest.cov_est()
    se = np.sqrt(np.diag(cov.to_numpy()))

    beta_parmest = np.array([theta["beta1"], theta["beta2"]])
    print(f"  theta_est() wall time: {t_theta:.2f} s ({n} experiments, 1 solve)")
    print(f"  ParmEst  beta = ({beta_parmest[0]:.6f}, {beta_parmest[1]:.6f}), "
          f"SSE = {sse_parmest:.6f}")
    print(f"  scipy    beta = ({scipy_fr.beta[0]:.6f}, {scipy_fr.beta[1]:.6f}), "
          f"SSE = {scipy_fr.sse:.6f}")
    print(f"  ParmEst  SE   = ({se[0]:.6f}, {se[1]:.6f})  [cov_est(), "
          f"measurement_error=None -> variance estimated from SSE/(n-p)]")
    print(f"  scipy    SE   = ({scipy_fr.se[0]:.6f}, {scipy_fr.se[1]:.6f})  "
          f"[sigma_hat^2 (J^T J)^-1]")
    rel_diff_beta = np.abs(beta_parmest - scipy_fr.beta) / scipy_fr.beta
    rel_diff_se = np.abs(se - scipy_fr.se) / scipy_fr.se
    print(f"  relative diff (beta): {rel_diff_beta}")
    print(f"  relative diff (SE):   {rel_diff_se}")
    print()
    return pest, beta_parmest, se, sse_parmest


# ---------------------------------------------------------------------------
# B.2: multi-start reliability with ParmEst/Ipopt, compare to scipy
# ---------------------------------------------------------------------------
def run_parmest_multistart(df, beta_global, sse_global, rng):
    print("=" * 78)
    print(f"PARMEST MULTI-START RELIABILITY STUDY ({N_MULTISTART} random starts, "
          f"same box as scipy)")
    print("=" * 78)

    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), N_MULTISTART)
    chi0 = np.exp(log_chi)
    log_delta = rng.uniform(np.log(DELTA_LO), np.log(DELTA_HI), N_MULTISTART)
    delta0 = np.exp(log_delta)

    outcomes = []
    betas = []
    times = []
    sse_tol_rel = 1e-6
    n = len(df)
    for chi0_i, delta0_i in zip(chi0, delta0):
        theta0 = {"beta1": float(chi0_i), "beta2": float(delta0_i)}
        exp_list = build_experiment_list(df, theta0=theta0)
        pest = parmest.Estimator(exp_list, obj_function="SSE")
        t0 = time.time()
        try:
            obj, theta = pest.theta_est()
            dt = time.time() - t0
            beta_hat = np.array([theta["beta1"], theta["beta2"]])
            sse_hat = obj * n
            if sse_hat <= sse_global * (1 + sse_tol_rel) + 1e-8:
                outcome = "global"
            else:
                outcome = "local/other"
        except Exception as exc:  # Ipopt failure / infeasible / etc.
            dt = time.time() - t0
            beta_hat = np.array([np.nan, np.nan])
            outcome = "failed"
        outcomes.append(outcome)
        betas.append(beta_hat)
        times.append(dt)

    outcomes = np.array(outcomes)
    betas = np.array(betas)
    n_global = int(np.sum(outcomes == "global"))
    n_local = int(np.sum(outcomes == "local/other"))
    n_failed = int(np.sum(outcomes == "failed"))
    print(f"  global: {100*n_global/N_MULTISTART:.1f}%   "
          f"local/other: {100*n_local/N_MULTISTART:.1f}%   "
          f"failed: {100*n_failed/N_MULTISTART:.1f}%")
    print(f"  mean wall time per start: {np.mean(times):.2f} s")
    print("  Comparison to scipy (Part A): lm/trf reached the global optimum "
          "from 100% of 300 random starts across the same box; only "
          "scipy's dogbox method (27% of starts) got trapped at a genuine "
          "|chi|->0 local minimum. Ipopt (interior-point, exact/AD "
          "derivatives, log-barrier handling of the beta2>0 bound) is being "
          "compared here against that baseline.")
    print()
    return dict(chi0=chi0, delta0=delta0, outcomes=outcomes, betas=betas, times=times)


def plot_parmest_multistart(ms, beta_global, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 5.5 / 6.5))
    color_map = {"global": "tab:green", "local/other": "tab:orange", "failed": "tab:red"}
    for outcome, color in color_map.items():
        mask = ms["outcomes"] == outcome
        if np.any(mask):
            ax.scatter(ms["chi0"][mask], ms["delta0"][mask], s=30, c=color,
                       alpha=0.8, label=f"{outcome} (n={mask.sum()})")
    ax.axhline(beta_global[1], color="k", ls=":", lw=1)
    ax.axvline(beta_global[0], color="k", ls=":", lw=1)
    ax.set_yscale("log")
    ax.set_xlabel(r"start $\boldsymbol{|\chi|_0}$ [mM]")
    ax.set_ylabel(r"start $\boldsymbol{\delta^*_0}$ [–]")
    ax.set_title(f"ParmEst/Ipopt multi-start outcomes ({N_MULTISTART} starts)")
    ax.legend(loc="best")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# B.3: LR/chi2 confidence region (objective_at_theta + likelihood_ratio_test)
# ---------------------------------------------------------------------------
def run_confidence_region(pest, beta_hat, se, sse_min, n, p, df):
    print("=" * 78)
    print(f"PARMEST LIKELIHOOD-RATIO CONFIDENCE REGION ({GRID_N}x{GRID_N} grid, "
          f"objective_at_theta)")
    print("=" * 78)

    b1_grid = np.linspace(beta_hat[0] - 12 * se[0], beta_hat[0] + 12 * se[0], GRID_N)
    b2_grid = np.linspace(beta_hat[1] - 12 * se[1], beta_hat[1] + 12 * se[1], GRID_N)
    b2_grid = b2_grid[b2_grid > 0]
    grid_pts = [(b1, b2) for b1 in b1_grid for b2 in b2_grid]
    theta_vals = pd.DataFrame(grid_pts, columns=["beta1", "beta2"])

    t0 = time.time()
    obj_at_theta = pest.objective_at_theta(theta_vals)
    dt = time.time() - t0
    print(f"  objective_at_theta() over {len(theta_vals)} points: {dt:.1f} s "
          f"({dt/len(theta_vals)*1000:.0f} ms/point)")

    alphas = [0.50, 0.90, 0.95]
    LR, thresholds = pest.likelihood_ratio_test(
        obj_at_theta, sse_min / n, alphas, return_thresholds=True
    )
    print(f"  chi2-based LR thresholds (mean-SSE units): {thresholds.to_dict()}")
    print()
    return LR, thresholds, b1_grid, b2_grid


def plot_confidence_comparison(LR, thresholds, b1_grid, b2_grid, beta_hat, cov,
                                sse_scipy_fr, c_int, c_p, y, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 6 / 7.5 * 1.5))

    # --- scipy nonlinear F-test region (reuse Part A machinery) ---
    SSE = s3.sse_surface(b1_grid, b2_grid, c_int, c_p, y)
    n, p = sse_scipy_fr.n, sse_scipy_fr.p
    colors = ["tab:blue", "tab:orange", "tab:red"]
    for alpha, color in zip([0.5, 0.10, 0.05], colors):
        thresh = s3.f_test_threshold(sse_scipy_fr.sse, n, p, alpha)
        conf_pct = int(round((1 - alpha) * 100))
        ax.contour(b1_grid, b2_grid, SSE.T, levels=[thresh], colors=[color],
                   linewidths=2, linestyles="solid")
        ax.plot([], [], color=color, lw=2, label=f"scipy F-test {conf_pct}%")

    # --- ParmEst LR/chi2 region: scatter grid points inside the 95% region ---
    for alpha, color in zip([0.50, 0.90, 0.95], colors):
        inside = LR[LR[alpha]]
        conf_pct = int(round(alpha * 100))
        ax.scatter(inside["beta1"], inside["beta2"], s=8, color=color, alpha=0.35,
                   marker="s")
    ax.scatter([], [], s=20, color="gray", alpha=0.6, marker="s",
               label="ParmEst LR/chi2 region (grid points inside)")

    ax.plot(beta_hat[0], beta_hat[1], marker="*", color="k", markersize=14,
            label="ParmEst theta_est() optimum")
    ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_title("ParmEst LR/chi2 region vs. scipy nonlinear F-test region\n"
                 "(NaCl corrected model, rows 57-748)")
    ax.xaxis.set_major_locator(plt.MaxNLocator(3))
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=1,
              borderaxespad=0, fontsize=7)
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# B.4: bootstrap
# ---------------------------------------------------------------------------
def run_bootstrap(pest, beta_hat, se_wald):
    print("=" * 78)
    print(f"PARMEST BOOTSTRAP ({N_BOOTSTRAP} reps, theta_est_bootstrap)")
    print("=" * 78)
    t0 = time.time()
    boot_theta = pest.theta_est_bootstrap(N_BOOTSTRAP, seed=RNG_SEED)
    dt = time.time() - t0
    print(f"  wall time: {dt:.1f} s ({dt/N_BOOTSTRAP*1000:.0f} ms/rep)")
    for col in ["beta1", "beta2"]:
        lo, hi = np.percentile(boot_theta[col], [2.5, 97.5])
        print(f"  {col:<8} bootstrap 95% percentile CI = ({lo:.5f}, {hi:.5f}), "
              f"bootstrap SE = {boot_theta[col].std(ddof=1):.5f}")
    print()
    return boot_theta


def plot_bootstrap(boot_theta, beta_hat, outfile):
    w = plot_style.fig_width("half")
    fig, ax = plt.subplots(figsize=(w, w * 5.5 / 6.5))
    ax.scatter(boot_theta["beta1"], boot_theta["beta2"], s=14, alpha=0.5,
               color="tab:purple", label=f"ParmEst bootstrap (n={len(boot_theta)})")
    ax.plot(beta_hat[0], beta_hat[1], marker="*", color="k", markersize=14,
            label="ParmEst theta_est() optimum")
    ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_title("ParmEst bootstrap replicates")
    ax.legend(loc="best")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 78)
    print("Step 3, Part B: Pyomo/ParmEst NaCl Donnan regression")
    print("=" * 78)

    df = load_primary_data()
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    n = len(df)

    scipy_fr = step2.fit_model(step2.f_corrected, c_int, c_p, y, beta0=step2.BETA0,
                               label="scipy baseline (rows 57-748, corrected model)")
    step2.print_fit(scipy_fr)

    pest, beta_parmest, se_parmest, sse_parmest = run_theta_est(df, scipy_fr)

    rng = np.random.default_rng(RNG_SEED)
    ms = run_parmest_multistart(df, beta_parmest, sse_parmest, rng)
    plot_parmest_multistart(ms, beta_parmest, FIG_DIR / "nacl_multistart_parmest.png")

    LR, thresholds, b1_grid, b2_grid = run_confidence_region(
        pest, beta_parmest, se_parmest, sse_parmest, n, 2, df,
    )
    cov = np.diag(se_parmest ** 2)
    plot_confidence_comparison(
        LR, thresholds, b1_grid, b2_grid, beta_parmest, cov, scipy_fr,
        c_int, c_p, y, FIG_DIR / "nacl_confidence_parmest.png",
    )

    boot_theta = run_bootstrap(pest, beta_parmest, se_parmest)
    plot_bootstrap(boot_theta, beta_parmest, FIG_DIR / "nacl_bootstrap_parmest.png")

    print("=" * 78)
    print("SUMMARY: ParmEst vs. scipy agreement")
    print("=" * 78)
    print(f"  |chi|:   ParmEst={beta_parmest[0]:.6f}  scipy={scipy_fr.beta[0]:.6f}  "
          f"diff={abs(beta_parmest[0]-scipy_fr.beta[0]):.2e}")
    print(f"  delta*:  ParmEst={beta_parmest[1]:.6f}  scipy={scipy_fr.beta[1]:.6f}  "
          f"diff={abs(beta_parmest[1]-scipy_fr.beta[1]):.2e}")
    print(f"  SE(|chi|):  ParmEst={se_parmest[0]:.6f}  scipy={scipy_fr.se[0]:.6f}")
    print(f"  SE(delta*): ParmEst={se_parmest[1]:.6f}  scipy={scipy_fr.se[1]:.6f}")
    print()


if __name__ == "__main__":
    main()
