"""
AR(1)-whitened null-check follow-up to the chi(c) spline first cut
(docs/prompts/chispline_ar1_whitened_nullcheck.md).

Question
--------
The first cut (chispline_cacl2.py) found the chi(c) spline overwhelmingly
beats constant-chi on CaCl2 (i.i.d. F-test), but the SAME machinery also
finds "significant" non-flat chi(c) on NaCl -- a dataset with no adsorption
mechanism, that should stay flat. Suspected cause: within one diafiltration
run, concentration rises ~monotonically with time, so chi(c) is nearly
interchangeable with chi(t), and an unpenalized spline can absorb NaCl's
AR(1) autocorrelation (phi ~ 0.74, sec:naclAR1) instead of chemistry. The
i.i.d. F-test/AIC used in the first cut do not know about that
autocorrelation and are anti-conservative (the same trap sec:naclAR1
diagnosed for the constant-chi regions).

Test
----
Redo the constant-vs-1-knot-spline comparison on AR(1)-WHITENED residuals,
for both NaCl (null check) and CaCl2:
  1. Fit constant-chi i.i.d.; estimate phi_hat from its time-ordered
     residuals via Cochrane-Orcutt (converged, "common" phi -- Step 8's
     Method A convention: GLS/whitening with actual n-p df, NOT combined
     with an n_eff deflation -- no double-counting).
  2. Whiten BOTH the constant and 1-knot-spline residual functions with
     that SAME common phi_hat, so the F-test isolates the mean-model change
     under one shared noise model (fair comparison).
  3. Nested F-test on the whitened SSE, actual n-p df.
  4. Sensitivity check: also let the spline re-estimate its OWN phi via
     Cochrane-Orcutt (mirrors Step 8b's Method A/B framing, but here it's
     common-phi vs per-model-phi, not whitening-vs-n_eff) -- reported
     alongside, not used in the primary F-test.

Verdict logic:
  - NaCl drift evaporates (whitened F not significant, chi(c) flattens)
    AND CaCl2 survives (still strongly favors the spline)
      => GREEN LIGHT: CaCl2 drift is real; autocorrelation was the NaCl
         confound.
  - NaCl improvement persists even whitened
      => confound is NOT just AR(1); needs pooling across differently-paced
         runs to decorrelate time and concentration.

Reuses (does not modify): analysis/step8_autocorrelation/ar1_correction.py
(lag1_autocorr, whiten -- both dataset-agnostic utilities already reused
by Step 9) and analysis/chi_spline/chispline_cacl2.py (f_chispline,
f_nacl_chispline, build_knots, eval_basis, chi_band, nested_f_test,
aic_bic, the CaCl2 dataset/constant-fit helpers). New Cochrane-Orcutt /
whitened-fit glue code lives here because Step 8's own cochrane_orcutt()
is hardwired to the 2-parameter NaCl deltaC_m model and Step 9's
generic version is hardwired to a module-level 2-parameter FIT_BOUNDS --
neither fits the spline's variable-length (n_basis+1)-parameter theta.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/chi_spline/chispline_ar1_whitened.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "cacl2"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step8_autocorrelation"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "chi_spline"))
import plot_style  # noqa: E402
import regress_cacl2 as cacl2  # noqa: E402
import regress_nacl as step2  # noqa: E402
import ar1_correction as s8  # noqa: E402
import chispline_cacl2 as s12  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

ALPHA = 0.05
RNG_SEED = 20260707


# ---------------------------------------------------------------------------
# Generic whitened fit / Cochrane-Orcutt, parameterized by (theta-shape,
# bounds) so it works for both the 2-param constant models and the
# (n_basis+1)-param spline models -- Step 8/9's versions are each hardwired
# to a fixed 2-param shape and can't be reused as-is (see module docstring).
# ---------------------------------------------------------------------------
def whitened_ls_fit(resid_fn, theta0, phi, bounds):
    def resid_w(theta):
        return s8.whiten(resid_fn(theta), phi)

    res = least_squares(resid_w, theta0, method="trf", bounds=bounds)
    sse = float(np.sum(res.fun ** 2))
    return res.x, sse, len(res.fun), len(res.x)


def cochrane_orcutt_fixed_bounds(resid_fn, theta0, bounds, max_iter=50, tol=1e-10):
    theta = np.asarray(theta0, dtype=float)
    r = resid_fn(theta)
    phi = s8.lag1_autocorr(r)
    n_iter = 0
    for n_iter in range(1, max_iter + 1):
        theta_new, sse, n, p = whitened_ls_fit(resid_fn, theta, phi, bounds)
        r_new = resid_fn(theta_new)
        phi_new = s8.lag1_autocorr(r_new)
        d_theta = np.max(np.abs(theta_new - theta) / np.maximum(np.abs(theta), 1e-12))
        d_phi = abs(phi_new - phi)
        theta, phi = theta_new, phi_new
        if d_theta < tol and d_phi < tol:
            break
    sse_final = float(np.sum(s8.whiten(resid_fn(theta), phi) ** 2))
    return theta, phi, n_iter, sse_final


# ---------------------------------------------------------------------------
# Per-dataset pipeline
# ---------------------------------------------------------------------------
def analyze_dataset(name, c_int, c_p, y, resid_fn_const, beta0_const,
                     bounds_const, chi0_flat, k_or_delta0, n_interior=1):
    n = len(y)

    # --- i.i.d. baselines (unwhitened), reusing existing fit code exactly
    #     as the first cut computed them ---
    res_const_iid = least_squares(resid_fn_const, beta0_const, method="trf",
                                   bounds=bounds_const)
    sse_const_iid = float(np.sum(res_const_iid.fun ** 2))
    p_const = len(res_const_iid.x)

    knots, n_basis, c_min, c_max = s12.build_knots(c_int, n_interior)
    s12.check_partition_of_unity(knots, n_basis, c_min, c_max)

    if name == "CaCl2":
        def resid_fn_spline(theta):
            return s12.f_chispline(theta, c_int, c_p, knots, n_basis, c_min, c_max) - y
    else:
        def resid_fn_spline(theta):
            return s12.f_nacl_chispline(theta, c_int, c_p, knots, n_basis, c_min, c_max) - y

    beta0_spline = np.concatenate([np.full(n_basis, chi0_flat), [k_or_delta0]])
    bounds_spline = (np.zeros(n_basis + 1), np.full(n_basis + 1, np.inf))

    res_spline_iid = least_squares(resid_fn_spline, beta0_spline, method="trf",
                                    bounds=bounds_spline)
    sse_spline_iid = float(np.sum(res_spline_iid.fun ** 2))
    p_spline = len(res_spline_iid.x)

    ftest_iid = s12.nested_f_test(sse_const_iid, p_const, sse_spline_iid, p_spline, n)
    aic_const_iid, _ = s12.aic_bic(sse_const_iid, n, p_const)
    aic_spline_iid, _ = s12.aic_bic(sse_spline_iid, n, p_spline)

    # --- AR(1): estimate phi_hat from the constant-chi model (Cochrane-Orcutt,
    #     Method A convention: whitening + actual n-p df, no n_eff double-count) ---
    phi_ols_initial = s8.lag1_autocorr(resid_fn_const(res_const_iid.x))
    theta_const_w, phi_common, n_iter_const, sse_const_w = cochrane_orcutt_fixed_bounds(
        resid_fn_const, res_const_iid.x, bounds_const,
    )

    # --- Whiten the SPLINE with the SAME (common) phi -- fair comparison,
    #     one shared noise model, single refit (not iterated) at that phi ---
    theta_spline_w_common, sse_spline_w_common, _, _ = whitened_ls_fit(
        resid_fn_spline, res_spline_iid.x, phi_common, bounds_spline,
    )

    ftest_whitened = s12.nested_f_test(sse_const_w, p_const, sse_spline_w_common,
                                        p_spline, n)
    aic_const_w, _ = s12.aic_bic(sse_const_w, n, p_const)
    aic_spline_w, _ = s12.aic_bic(sse_spline_w_common, n, p_spline)

    # --- Sensitivity: let the spline re-estimate its OWN phi (Cochrane-Orcutt),
    #     informational only -- not used in the primary (common-phi) F-test ---
    theta_spline_w_own, phi_spline_own, n_iter_spline, sse_spline_w_own = \
        cochrane_orcutt_fixed_bounds(resid_fn_spline, res_spline_iid.x, bounds_spline)

    betas_common = theta_spline_w_common[:-1]
    spread_common = float(np.max(betas_common) - np.min(betas_common))
    rel_spread_common = spread_common / chi0_flat if chi0_flat else np.nan

    betas_own = theta_spline_w_own[:-1]
    spread_own = float(np.max(betas_own) - np.min(betas_own))
    rel_spread_own = spread_own / chi0_flat if chi0_flat else np.nan

    print("=" * 100)
    print(f"{name}: i.i.d. vs. AR(1)-whitened, constant vs. {n_interior}-knot spline")
    print("=" * 100)
    print(f"n = {n}")
    print(f"i.i.d.:    constant SSE = {sse_const_iid:.4f} (p={p_const}), "
          f"spline SSE = {sse_spline_iid:.4f} (p={p_spline}), "
          f"F({ftest_iid['df1']},{ftest_iid['df2']}) = {ftest_iid['F']:.3f}, "
          f"p = {ftest_iid['p_value']:.5f}, dAIC = {aic_spline_iid - aic_const_iid:+.2f} "
          f"-> {'SIGNIFICANT' if ftest_iid['significant'] else 'not significant'}")
    print(f"phi_hat: OLS lag-1 (pre-Cochrane-Orcutt) = {phi_ols_initial:.4f}; "
          f"Cochrane-Orcutt converged ({n_iter_const} iters) common phi = {phi_common:.4f}")
    print(f"Whitened (common phi={phi_common:.4f}): constant SSE = {sse_const_w:.4f}, "
          f"spline SSE = {sse_spline_w_common:.4f}, "
          f"F({ftest_whitened['df1']},{ftest_whitened['df2']}) = {ftest_whitened['F']:.3f}, "
          f"p = {ftest_whitened['p_value']:.5f}, dAIC = {aic_spline_w - aic_const_w:+.2f} "
          f"-> {'SIGNIFICANT' if ftest_whitened['significant'] else 'not significant'}")
    print(f"  chi(c) control-point spread (common phi) = {rel_spread_common:.4f} of chi0 "
          f"({'FLAT' if rel_spread_common < 0.05 else 'not flat'})")
    print(f"Sensitivity -- spline's OWN Cochrane-Orcutt phi (informational, "
          f"{n_iter_spline} iters): phi_spline_own = {phi_spline_own:.4f} "
          f"(vs. common phi = {phi_common:.4f}, delta = {phi_spline_own - phi_common:+.4f})")
    print(f"  chi(c) control-point spread (own phi) = {rel_spread_own:.4f} of chi0 "
          f"({'FLAT' if rel_spread_own < 0.05 else 'not flat'})")
    print()

    return dict(
        name=name, n=n, knots=knots, n_basis=n_basis, c_min=c_min, c_max=c_max,
        sse_const_iid=sse_const_iid, sse_spline_iid=sse_spline_iid,
        ftest_iid=ftest_iid, aic_const_iid=aic_const_iid, aic_spline_iid=aic_spline_iid,
        phi_ols_initial=phi_ols_initial,
        phi_common=phi_common, sse_const_w=sse_const_w,
        sse_spline_w_common=sse_spline_w_common, ftest_whitened=ftest_whitened,
        aic_const_w=aic_const_w, aic_spline_w=aic_spline_w,
        rel_spread_common=rel_spread_common, theta_spline_w_common=theta_spline_w_common,
        phi_spline_own=phi_spline_own, rel_spread_own=rel_spread_own,
        theta_spline_iid=res_spline_iid.x, chi0_flat=chi0_flat,
        beta_const_iid=res_const_iid.x, is_cacl2=(name == "CaCl2"),
    )


# ---------------------------------------------------------------------------
# Data loading (identical to the first cut / Step 2, time-ordered)
# ---------------------------------------------------------------------------
def load_nacl():
    df_full = step2.load_sheet()  # "(2)" sheet, our-CP (APPLY_CP=True)
    df_full = df_full.loc[(df_full.index >= step2.MATLAB_FIRST_ROW)
                           & (df_full.index <= step2.MATLAB_LAST_ROW)].dropna()
    df_full = step2.build_deltaC_m(df_full)
    return (df_full["c_int"].to_numpy(), df_full["c_p"].to_numpy(),
            df_full["deltaC_m"].to_numpy())


def load_cacl2():
    df = cacl2.load_dataset(s12.CACL2_DATASET)
    return df["c_int"].to_numpy(), df["c_p"].to_numpy(), df["deltaC_m"].to_numpy()


# ---------------------------------------------------------------------------
# Figure: whitened chi(c) for NaCl (should flatten) and CaCl2 (should survive)
# ---------------------------------------------------------------------------
def make_figure(result_nacl, result_cacl2, outfile):
    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 2, figsize=(w, w * 0.52))

    for ax, res in zip(axes, [result_nacl, result_cacl2]):
        c_grid = np.linspace(res["c_min"], res["c_max"], 300)
        B = s12.eval_basis(res["knots"], res["n_basis"], res["c_min"], res["c_max"], c_grid)

        chi_iid = B @ res["theta_spline_iid"][:-1]
        chi_whitened = B @ res["theta_spline_w_common"][:-1]

        ax.axhline(res["chi0_flat"], color="k", ls=":", lw=1.2, label="constant $\\chi$")
        ax.plot(c_grid, chi_iid, "--", color=plot_style.OKABE_ITO["blue"], lw=2,
                label="i.i.d. 1-knot spline")
        ax.plot(c_grid, chi_whitened, "-", color=plot_style.OKABE_ITO["vermillion"],
                lw=2, label="AR(1)-whitened 1-knot spline")
        ax.set_xlabel(r"$\boldsymbol{c}$ [mM]")
        ax.set_ylabel(r"$\boldsymbol{\chi(c)}$ [mM]")
        verdict = "flattens" if res["rel_spread_common"] < 0.05 else "still not flat"
        ax.set_title(f"{res['name']} ($\\hat\\phi$={res['phi_common']:.2f}, "
                     f"whitened spread {verdict})")
        ax.set_box_aspect(1)

    axes[0].legend(loc="upper center", bbox_to_anchor=(1.05, -0.38), ncol=3)

    fig.suptitle("AR(1)-whitened $\\boldsymbol{\\chi(c)}$ null-check: NaCl (should flatten) vs. "
                 "CaCl2 (should survive)")
    fig.tight_layout(rect=(0, 0.05, 1, 0.92))
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("#" * 100)
    print("AR(1)-whitened chi(c) null-check (NaCl) / survival-check (CaCl2)")
    print("#" * 100)
    print()

    # --- NaCl ---
    c_int_n, c_p_n, y_n = load_nacl()

    def resid_fn_const_nacl(beta):
        return step2.f_corrected(beta, c_int_n, c_p_n) - y_n

    # (chi,delta enter step2.f_corrected squared, so bounding at >=0 does not
    # change the optimum; kept >=0 here only so the warm-start convention
    # matches CaCl2's bounded fit exactly.)
    res_const_iid = least_squares(resid_fn_const_nacl, step2.BETA0, method="trf",
                                   bounds=([0.0, 0.0], [np.inf, np.inf]))
    chi0_n, delta0_n = res_const_iid.x
    result_nacl = analyze_dataset(
        "NaCl", c_int_n, c_p_n, y_n, resid_fn_const_nacl, step2.BETA0,
        bounds_const=([0.0, 0.0], [np.inf, np.inf]),
        chi0_flat=chi0_n, k_or_delta0=delta0_n,
    )

    # --- CaCl2 ---
    c_int_c, c_p_c, y_c = load_cacl2()

    def resid_fn_const_cacl2(beta):
        return cacl2.f_cacl2(beta, c_int_c, c_p_c) - y_c

    res_const_iid_c = least_squares(resid_fn_const_cacl2, cacl2.BETA0, method="trf",
                                     bounds=cacl2.FIT_BOUNDS)
    chi0_c, K0_c = res_const_iid_c.x
    result_cacl2 = analyze_dataset(
        "CaCl2", c_int_c, c_p_c, y_c, resid_fn_const_cacl2, cacl2.BETA0,
        bounds_const=cacl2.FIT_BOUNDS, chi0_flat=chi0_c, k_or_delta0=K0_c,
    )

    print("=" * 100)
    print("SUMMARY TABLE")
    print("=" * 100)
    header = (f"{'dataset':<8}{'regime':<10}{'F':>10}{'df1,df2':>12}{'p-value':>12}"
              f"{'dAIC':>10}{'chi(c) spread':>16}{'verdict':>16}")
    print(header)
    for res in [result_nacl, result_cacl2]:
        ft = res["ftest_iid"]
        df_str_iid = f"{ft['df1']},{ft['df2']}"
        print(f"{res['name']:<8}{'i.i.d.':<10}{ft['F']:>10.2f}"
              f"{df_str_iid:>12}{ft['p_value']:>12.5f}"
              f"{res['aic_spline_iid']-res['aic_const_iid']:>10.2f}"
              f"{'--':>16}"
              f"{'SIGNIFICANT' if ft['significant'] else 'not sig.':>16}")
        ftw = res["ftest_whitened"]
        df_str_w = f"{ftw['df1']},{ftw['df2']}"
        print(f"{res['name']:<8}{'whitened':<10}{ftw['F']:>10.2f}"
              f"{df_str_w:>12}{ftw['p_value']:>12.5f}"
              f"{res['aic_spline_w']-res['aic_const_w']:>10.2f}"
              f"{res['rel_spread_common']:>16.4f}"
              f"{'SIGNIFICANT' if ftw['significant'] else 'not sig.':>16}")
    print()

    nacl_evaporates = not result_nacl["ftest_whitened"]["significant"] \
        or result_nacl["rel_spread_common"] < 0.10
    cacl2_survives = result_cacl2["ftest_whitened"]["significant"] \
        and result_cacl2["rel_spread_common"] >= 0.10

    print("-" * 100)
    print("VERDICT")
    print("-" * 100)
    if nacl_evaporates and cacl2_survives:
        print("GREEN LIGHT: the NaCl null-check drift evaporates (or is much reduced) "
              "under AR(1) whitening, and CaCl2's chi(c) survives -- the CaCl2 drift "
              "looks like real concentration structure, not an autocorrelation artifact.")
    elif not nacl_evaporates:
        print("CONFOUND PERSISTS: NaCl's spurious drift does NOT evaporate under AR(1) "
              "whitening -- the time-concentration confound is not just AR(1) and must "
              "be broken by pooling across differently-paced runs.")
    else:
        print("MIXED/INCONCLUSIVE: see the numbers above -- neither branch cleanly holds.")
    print()

    make_figure(result_nacl, result_cacl2, FIG_DIR / "cacl2_chispline_whitened.png")

    return dict(nacl=result_nacl, cacl2=result_cacl2)


if __name__ == "__main__":
    main()
