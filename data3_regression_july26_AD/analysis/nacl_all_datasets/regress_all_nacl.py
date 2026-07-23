"""
Regression across all available processed NaCl (1:1) datasets, to assess
reproducibility of the Donnan membrane parameters |chi| (fixed charge, mM)
and delta* (partition factor) across (a) different membrane cuts (MC2 vs
MC5) and (b) a repeat run on the same cut (MC5 NaCl vs NaCl (2)).

Reuses Step 2 (analysis/step2_nacl/regress_nacl.py: build_deltaC_m, f_matlab,
f_corrected, fit_model, load_sheet) and Step 3
(analysis/step3_nacl_uncertainty/scipy_uncertainty.py: f_test_threshold,
sse_surface, profile_likelihood, profile_ci_from_curve) unmodified -- no
model/fitting logic is duplicated here.

Primary model: f_corrected (the physically correct 1/2-form), fit over the
FULL valid data range per dataset (there is no equivalent of Bill's ad hoc
57-748 MATLAB window for MC2/MC5, so full-range is the primary, comparable
fit across datasets). A "startup-trimmed" fit (drop the first
STARTUP_TRIM_FRAC of rows, the low-concentration startup transient) is
reported as a sensitivity check. The uncorrected f_matlab estimate is also
reported (secondary column) for continuity with Bill's original values.

Run with (from the repo root, in the documented conda environment):
    conda activate data3-regression
    python analysis/nacl_all_datasets/regress_all_nacl.py
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
from scipy.stats import norm

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step3_nacl_uncertainty"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import scipy_uncertainty as s3  # noqa: E402  (reuse F-test / profile utilities)

# CP Phase 2 (docs/prompts/cp_phase2_regressions.md): apply concentration
# polarization as the new standing convention. Default primary run is
# CP=on; APPLY_CP=False reproduces the pre-Phase-2 (bulk == interfacial)
# numbers byte-for-byte (regression-tested).
APPLY_CP = True

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

RNG_SEED = 20260705
N_MULTISTART = 100
STARTUP_TRIM_FRAC = 0.07  # drop the first ~7% of (sorted-by-row) data as the
                           # low-concentration startup transient, mirroring
                           # (in spirit, not row count) Bill's MATLAB window
                           # which skipped rows 2-56 of 2-748 (~7.4%).

CHI_LO, CHI_HI = 1.0, 200.0
DELTA_LO, DELTA_HI = 1e-3, 10.0
ALPHA = 0.05  # 95% CIs/regions throughout

# Reference value independently used in the group's back-of-envelope (BoE)
# calculations (see docs/reports/main.tex, Sec. "Potential modeling issue").
CHI_REFERENCE_MM = 44.0

PREPROC_DIR = REPO_ROOT / "analysis" / "preprocessing" / "processed"

DATASETS = [
    # Excel-sourced (already-processed sheets in data/BoE Analysis.xlsx).
    # `already_cp`: whether column B already reflects the AUTHORS' OWN CP
    # correction (do not re-apply ours -- double-applying would be wrong) vs.
    # is their non-CP bulk value (our thin-film CP is applied on top, see
    # cp_authors_bulk.py). Classification is by the c_int/kappa_ret
    # coefficient-of-variation diagnostic (CP Phase 2 refinement audit):
    # MC2 CV=0.096, MC5(2) CV=0.088 -- both in the non-CP cluster, so BOTH
    # get already_cp=False (MC2 was incorrectly marked True in the initial
    # Phase-2 pass, silently skipping CP for it -- fixed here); MC5-main
    # CV=0.255, clearly the odd one out -- genuinely the authors' own CP.
    # `role`: "primary" datasets all carry OUR CP uniformly (MC2, MC5(2), and
    # the 6 raw-only CSVs below); MC5-main is the lone AUTHORS'-CP dataset,
    # kept as an explicitly-labeled "comparison" (not a second primary MC5
    # row) so the table doesn't read as an unexplained MC5 != MC5(2).
    {"key": "MC2", "source": "excel", "sheet": "NF270_MC2 05.07.24_NaCl",
     "label": "MC2 (05.07.24)", "already_cp": False, "role": "primary"},
    {"key": "MC5", "source": "excel", "sheet": "NF270_MC5 05.27.26_NaCl",
     "label": "MC5 (05.27.26, authors' CP) [comparison]", "already_cp": True,
     "role": "comparison"},
    {"key": "MC5(2)", "source": "excel", "sheet": "NF270_MC5 05.27.26_NaCl (2)",
     "label": "MC5 (05.27.26)", "already_cp": False, "role": "primary"},
    # CSV-sourced: raw-only runs processed by
    # analysis/preprocessing/{preprocess_nacl,process_raw_runs}.py. Non-CP
    # (bulk == interfacial) CSVs listed here; CP Phase 2 swaps to the
    # Phase-1 "*_CP.csv" twins (vial-ICP calibration + CP) when APPLY_CP.
    {"key": "MC3-0709", "source": "csv", "path": PREPROC_DIR / "MC3_07.09.24_NaCl.csv",
     "path_cp": PREPROC_DIR / "MC3_07.09.24_NaCl_CP.csv", "label": "MC3 (07.09.24)"},
    {"key": "MC3-0722S", "source": "csv", "path": PREPROC_DIR / "MC3_07.22.24_SNaCl.csv",
     "path_cp": PREPROC_DIR / "MC3_07.22.24_SNaCl_CP.csv", "label": "MC3 (07.22.24, S)"},
    {"key": "MC4-0711S", "source": "csv", "path": PREPROC_DIR / "MC4_07.11.24_SNaCl.csv",
     "path_cp": PREPROC_DIR / "MC4_07.11.24_SNaCl_CP.csv", "label": "MC4 (07.11.24, S)"},
    {"key": "MC5-0723", "source": "csv", "path": PREPROC_DIR / "MC5_07.23.24_NaCl.csv",
     "path_cp": PREPROC_DIR / "MC5_07.23.24_NaCl_CP.csv", "label": "MC5 (07.23.24)"},
    {"key": "MC5-0723S", "source": "csv", "path": PREPROC_DIR / "MC5_07.23.24_SNaCl.csv",
     "path_cp": PREPROC_DIR / "MC5_07.23.24_SNaCl_CP.csv", "label": "MC5 (07.23.24, S)"},
    {"key": "MC5-0723S2", "source": "csv", "path": PREPROC_DIR / "MC5_07.23.24_S2NaCl.csv",
     "path_cp": PREPROC_DIR / "MC5_07.23.24_S2NaCl_CP.csv", "label": "MC5 (07.23.24, S2)"},
]
DATASET_COLORS = plot_style.get_dataset_colors([d["label"] for d in DATASETS])


# ---------------------------------------------------------------------------
# Data loading (reuses step2.load_sheet / build_deltaC_m for Excel-sourced
# datasets; CSV-sourced datasets from analysis/preprocessing/ already give
# c_int/c_p directly and need deltaC_m built from c_p, Jw with the same
# constants Step 2 uses, since they don't share the raw c_int/c_p/J_w
# spreadsheet-column schema build_deltaC_m expects)
# ---------------------------------------------------------------------------
def load_dataset(ds: dict):
    if ds["source"] == "excel":
        df_full = step2.load_sheet(sheet_name=ds["sheet"], apply_cp=APPLY_CP,
                                    already_cp=ds.get("already_cp", False))
        valid = df_full.notna().all(axis=1)
        df_valid = df_full[valid].copy()
        df_valid = step2.build_deltaC_m(df_valid)
    else:
        csv_path = ds["path_cp"] if (APPLY_CP and "path_cp" in ds) else ds["path"]
        df_raw = pd.read_csv(csv_path, comment="#")
        df_full = pd.DataFrame({
            "c_int": df_raw["c_int_mM"],
            "c_p": df_raw["c_p_mM"],
            "deltaC_m": df_raw["c_p_mM"] * df_raw["Jw_m3_m2_s"]
                        * step2.L_MEMBRANE / step2.D_NACL_M,
        })
        # Physically invalid rows (c_p >= c_int, e.g. very start of the run
        # where the retentate/permeate calibration noise dominates) can't
        # feed the f_corrected model (negative under the sqrt argument is
        # fine, but c_p > c_int reverses the sign expected of solute
        # rejection) -- drop them, matching how the Excel-sourced sheets
        # never contain such rows.
        valid = df_full.notna().all(axis=1) & (df_full["c_int"] > df_full["c_p"])
        df_valid = df_full[valid].copy()
    n_trim = int(round(len(df_valid) * STARTUP_TRIM_FRAC))
    df_trimmed = df_valid.iloc[n_trim:].copy()
    return df_valid, df_trimmed, n_trim


# ---------------------------------------------------------------------------
# Multi-start reliability (method='lm', per Step 3's finding that lm/trf are
# reliable while dogbox can get trapped at |chi|->0 -- see step3 README)
# ---------------------------------------------------------------------------
def multistart_reliability(c_int, c_p, y, beta_global, sse_global, rng,
                            n_starts=N_MULTISTART):
    log_chi = rng.uniform(np.log(CHI_LO), np.log(CHI_HI), n_starts)
    log_delta = rng.uniform(np.log(DELTA_LO), np.log(DELTA_HI), n_starts)
    chi0s, delta0s = np.exp(log_chi), np.exp(log_delta)

    n_global = 0
    for chi0, delta0 in zip(chi0s, delta0s):
        def resid(beta):
            return step2.f_corrected(beta, c_int, c_p) - y

        try:
            res = least_squares(resid, np.array([chi0, delta0]), method="lm")
            sse_hat = float(np.sum(res.fun ** 2))
            if res.success and sse_hat <= sse_global * (1 + 1e-6) + 1e-8:
                n_global += 1
        except Exception:
            pass
    return 100.0 * n_global / n_starts


# ---------------------------------------------------------------------------
# Per-dataset analysis
# ---------------------------------------------------------------------------
def analyze_dataset(ds, rng):
    label = ds["label"]
    df_full, df_trim, n_trim = load_dataset(ds)

    c_int_full = df_full["c_int"].to_numpy()
    c_p_full = df_full["c_p"].to_numpy()
    y_full = df_full["deltaC_m"].to_numpy()

    c_int_trim = df_trim["c_int"].to_numpy()
    c_p_trim = df_trim["c_p"].to_numpy()
    y_trim = df_trim["deltaC_m"].to_numpy()

    # Primary: corrected model, full valid range
    fr_corrected_full = step2.fit_model(
        step2.f_corrected, c_int_full, c_p_full, y_full, beta0=step2.BETA0,
        label=f"{label}: corrected, full range",
    )
    # Secondary: uncorrected model (Bill's form), full valid range, for
    # continuity with his originally-reported numbers.
    fr_uncorrected_full = step2.fit_model(
        step2.f_matlab, c_int_full, c_p_full, y_full, beta0=step2.BETA0,
        label=f"{label}: uncorrected, full range",
    )
    # Sensitivity: corrected model, startup-trimmed range
    fr_corrected_trim = step2.fit_model(
        step2.f_corrected, c_int_trim, c_p_trim, y_trim, beta0=step2.BETA0,
        label=f"{label}: corrected, startup-trimmed",
    )

    # Sanity check: the raw-only runs are noisier (Jw smoothing, fitted
    # calibration, no CP treatment) and some show a near-zero apparent
    # rejection (c_int approx c_p over much of the run) that leaves the
    # 2-parameter Donnan model very poorly conditioned -- least_squares
    # can wander off to a nonphysical optimum (|chi| >> O(100) mM). Detect
    # that rather than let it silently poison the multi-start / profile-CI
    # steps (which assume a locally-quadratic SSE surface near a physical
    # optimum) or crash them outright.
    beta_hat_full = fr_corrected_full.beta
    ill_conditioned = (not np.all(np.isfinite(beta_hat_full))
                        or beta_hat_full[0] > 1e4 or beta_hat_full[1] > 1e4)

    if ill_conditioned:
        pct_global = float("nan")
        profile_cis = [(float("nan"), float("nan")), (float("nan"), float("nan"))]
    else:
        pct_global = multistart_reliability(
            c_int_full, c_p_full, y_full, fr_corrected_full.beta,
            fr_corrected_full.sse, rng,
        )

        # 95% profile-likelihood CIs (primary model, full range) -- reuses
        # step3's profile_likelihood/profile_ci_from_curve unmodified.
        beta_hat, se = fr_corrected_full.beta, fr_corrected_full.se
        n, p, sse_min = fr_corrected_full.n, fr_corrected_full.p, fr_corrected_full.sse
        profile_cis = []
        for idx in range(2):
            grid, prof_sse, thresh = s3.profile_likelihood(
                beta_hat, se, sse_min, n, p, c_int_full, c_p_full, y_full, idx,
                alpha=ALPHA,
            )
            lo, hi = s3.profile_ci_from_curve(grid, prof_sse, thresh)
            profile_cis.append((lo, hi))

    return dict(
        ds=ds,
        df_full=df_full,
        n_trim=n_trim,
        fr_corrected_full=fr_corrected_full,
        fr_uncorrected_full=fr_uncorrected_full,
        fr_corrected_trim=fr_corrected_trim,
        ill_conditioned=ill_conditioned,
        pct_global=pct_global,
        profile_cis=profile_cis,
    )


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------
def print_summary_table(results):
    z = norm.ppf(1 - ALPHA / 2)
    print("=" * 100)
    print("SUMMARY: NaCl (1:1) Donnan regression across all processed datasets "
          "(primary = f_corrected, full range)")
    print("=" * 100)
    header = (f"{'Dataset':<18}{'n':>6}{'|chi| (mM)':>14}{'95% CI':>22}"
              f"{'delta*':>10}{'95% CI':>22}{'SSE':>12}{'MS % global':>13}")
    print(header)
    for r in results:
        fr = r["fr_corrected_full"]
        if r["ill_conditioned"]:
            print(f"{r['ds']['label']:<18}{fr.n:>6}  ILL-CONDITIONED FIT "
                  f"(|chi|={fr.beta[0]:.3g} mM, delta*={fr.beta[1]:.3g} -- "
                  f"nonphysical; SSE={fr.sse:.3f}); skipped multi-start/profile CI")
            continue
        lo1, hi1 = r["profile_cis"][0]
        lo2, hi2 = r["profile_cis"][1]
        print(f"{r['ds']['label']:<18}{fr.n:>6}{fr.beta[0]:>14.4f}"
              f"  ({lo1:6.3f}, {hi1:6.3f})"
              f"{fr.beta[1]:>10.5f}"
              f"  ({lo2:6.4f}, {hi2:6.4f})"
              f"{fr.sse:>12.3f}{r['pct_global']:>12.1f}%")
    print()

    print("-" * 100)
    print("Secondary: uncorrected f_matlab estimates (continuity with Bill), full range")
    print("-" * 100)
    for r in results:
        fr = r["fr_uncorrected_full"]
        print(f"{r['ds']['label']:<18}{fr.n:>6}  |chi|={fr.beta[0]:>9.4f} +/- {fr.se[0]:<8.4f}"
              f"  delta*={fr.beta[1]:>9.6f} +/- {fr.se[1]:<9.6f}  SSE={fr.sse:.3f}")
    print()

    print("-" * 100)
    print(f"Startup-trimmed sensitivity (dropped first {STARTUP_TRIM_FRAC*100:.0f}% "
          "of rows, corrected model)")
    print("-" * 100)
    for r in results:
        fr_full = r["fr_corrected_full"]
        fr_trim = r["fr_corrected_trim"]
        d_chi_pct = 100 * (fr_trim.beta[0] - fr_full.beta[0]) / fr_full.beta[0]
        d_delta_pct = 100 * (fr_trim.beta[1] - fr_full.beta[1]) / fr_full.beta[1]
        print(f"{r['ds']['label']:<18}n_trim={r['n_trim']:>4}  "
              f"full: |chi|={fr_full.beta[0]:.4f}, delta*={fr_full.beta[1]:.5f}  "
              f"trimmed: |chi|={fr_trim.beta[0]:.4f} ({d_chi_pct:+.2f}%), "
              f"delta*={fr_trim.beta[1]:.5f} ({d_delta_pct:+.2f}%)")
    print()

    # Cross-dataset comparison
    print("=" * 100)
    print("CROSS-DATASET COMPARISON")
    print("=" * 100)
    by_key = {r["ds"]["key"]: r for r in results}
    mc5_authors_cp = by_key["MC5"]["fr_corrected_full"]     # comparison role (authors' own CP)
    mc5_primary = by_key["MC5(2)"]["fr_corrected_full"]     # primary role (our CP, on the "(2)" sheet's bulk)
    mc2 = by_key["MC2"]["fr_corrected_full"]                # primary role (our CP)

    def overlap(lo1, hi1, lo2, hi2):
        return not (hi1 < lo2 or hi2 < lo1)

    mc5_chi_overlap = overlap(*by_key["MC5"]["profile_cis"][0], *by_key["MC5(2)"]["profile_cis"][0])
    mc5_delta_overlap = overlap(*by_key["MC5"]["profile_cis"][1], *by_key["MC5(2)"]["profile_cis"][1])
    mc2_mc5_chi_overlap = overlap(*by_key["MC2"]["profile_cis"][0], *by_key["MC5(2)"]["profile_cis"][0])
    mc2_mc5_delta_overlap = overlap(*by_key["MC2"]["profile_cis"][1], *by_key["MC5(2)"]["profile_cis"][1])

    print(f"MC5 (05.27.26): our-CP primary (MC5(2)+our-CP) |chi|={mc5_primary.beta[0]:.3f} vs "
          f"authors'-CP comparison |chi|={mc5_authors_cp.beta[0]:.3f} mM, "
          f"95% CIs overlap = {mc5_chi_overlap}; "
          f"delta* {mc5_primary.beta[1]:.5f} vs {mc5_authors_cp.beta[1]:.5f}, "
          f"95% CIs overlap = {mc5_delta_overlap}")
    if APPLY_CP:
        chi_reldiff = abs(mc5_authors_cp.beta[0] - mc5_primary.beta[0]) / mc5_authors_cp.beta[0] * 100
        delta_reldiff = abs(mc5_authors_cp.beta[1] - mc5_primary.beta[1]) / mc5_authors_cp.beta[1] * 100
        print(f"  NOTE: this is the OUR-CP-vs-AUTHORS'-CP implementation")
        print(f"  discrepancy (|chi| rel. diff {chi_reldiff:.2f}%, delta* rel. diff "
              f"{delta_reldiff:.2f}%) -- a third CP-related sensitivity axis "
              f"alongside CP-vs-non-CP (discussion point 11) and calibration "
              f"source (point 9, sec:condmodels). It arises because our "
              f"Leveque-type thin-film k differs from whatever CP mechanism "
              f"the authors applied (Phase 1 found our-CP reproduces the "
              f"authors' CP MC5 sheet only to ~12% median from raw "
              f"reconstruction -- consistent with this ~13% parameter-level "
              f"gap). A minor open item, not a correctness bug.")
    print(f"MC2 vs MC5 (05.27.26) [both our-CP primary]: |chi| {mc2.beta[0]:.3f} vs "
          f"{mc5_primary.beta[0]:.3f} mM, "
          f"95% CIs overlap = {mc2_mc5_chi_overlap}; "
          f"delta* {mc2.beta[1]:.5f} vs {mc5_primary.beta[1]:.5f}, "
          f"95% CIs overlap = {mc2_mc5_delta_overlap}")
    for r in results:
        if r["ill_conditioned"]:
            continue
        fr = r["fr_corrected_full"]
        lo, hi = r["profile_cis"][0]
        contains_ref = lo <= CHI_REFERENCE_MM <= hi
        print(f"{r['ds']['label']:<18}|chi|={fr.beta[0]:.3f} mM vs BoE reference "
              f"{CHI_REFERENCE_MM:.0f} mM: 95% CI ({lo:.3f}, {hi:.3f}) "
              f"{'CONTAINS' if contains_ref else 'does NOT contain'} the reference")
    print()

    csv_results = [r for r in results if r["ds"]["source"] == "csv"]
    if csv_results:
        print("-" * 100)
        print("New (preprocessed raw-only) runs vs. the MC5 (05.27.26) our-CP primary baseline:")
        print("-" * 100)
        for r in csv_results:
            fr = r["fr_corrected_full"]
            if r["ill_conditioned"]:
                print(f"{r['ds']['label']:<18}ILL-CONDITIONED FIT (|chi|={fr.beta[0]:.3g} mM, "
                      f"delta*={fr.beta[1]:.3g}) -- see README for likely cause "
                      "(apparent near-zero rejection in the reconstructed c_int/c_p)")
                continue
            chi_overlap = overlap(*r["profile_cis"][0], *by_key["MC5(2)"]["profile_cis"][0])
            print(f"{r['ds']['label']:<18}|chi|={fr.beta[0]:>8.3f} mM  "
                  f"delta*={fr.beta[1]:>9.5f}   95% CI overlaps MC5 (05.27.26) = {chi_overlap}")
        print()


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
def plot_all_estimates(results, outfile):
    # Displayed at 0.6\linewidth in main.tex (not the "half"=0.49\linewidth
    # minipage convention used elsewhere) -- save at that same fraction of
    # the full page width so on-page text renders at its designed ~11pt
    # scale (see figure-compliance audit, docs/prompts/figure_compliance_audit.md).
    w = plot_style.fig_width("full") * 0.6
    fig, axes = plt.subplots(2, 1, figsize=(w, w * 6.4 / 4.0), sharex=True)
    labels = [r["ds"]["label"] for r in results]
    x = np.arange(len(labels))

    ok = [not r["ill_conditioned"] for r in results]
    colors = [DATASET_COLORS[label] for label in labels]

    ax = axes[0]
    for xi, r, c, keep in zip(x, results, colors, ok):
        fr = r["fr_corrected_full"]
        if keep:
            lo = fr.beta[0] - r["profile_cis"][0][0]
            hi = r["profile_cis"][0][1] - fr.beta[0]
            ax.errorbar(xi, fr.beta[0], yerr=[[lo], [hi]], fmt="o", color=c,
                        capsize=5, markersize=8, elinewidth=2)
        else:
            ax.plot(xi, min(fr.beta[0], CHI_REFERENCE_MM * 3), marker="x", color=c,
                    markersize=10, markeredgewidth=2)
    ax.axhline(CHI_REFERENCE_MM, color="gray", ls="--", lw=1.5,
               label=f"BoE reference ({CHI_REFERENCE_MM:.0f} mM)")
    ax.set_ylabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax.set_ylim(top=CHI_REFERENCE_MM * 3.2)
    ax.legend(loc="best", fontsize=9)

    ax = axes[1]
    for xi, r, c, keep in zip(x, results, colors, ok):
        fr = r["fr_corrected_full"]
        if keep:
            lo = fr.beta[1] - r["profile_cis"][1][0]
            hi = r["profile_cis"][1][1] - fr.beta[1]
            ax.errorbar(xi, fr.beta[1], yerr=[[lo], [hi]], fmt="o", color=c,
                        capsize=5, markersize=8, elinewidth=2)
        else:
            ax.plot(xi, ax.get_ylim()[1] if ax.get_ylim()[1] > 0 else 1, marker="x",
                    color=c, markersize=10, markeredgewidth=2)
    ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=20, ha="right")

    for ax in axes:
        ax.label_outer()

    fig.suptitle("NaCl point estimates $\\boldsymbol{\\pm}$ 95% profile-likelihood CI\n"
                 "(corrected model, full range)")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


WELL_IDENTIFIED_REL_CI_WIDTH = 0.15  # (hi-lo)/point < 15% on |chi| => "well-identified"


def is_well_identified(r):
    if r["ill_conditioned"]:
        return False
    lo, hi = r["profile_cis"][0]
    point = r["fr_corrected_full"].beta[0]
    return (hi - lo) / point < WELL_IDENTIFIED_REL_CI_WIDTH


def _draw_regions(ax, results, linestyles=("solid", "dashed", "dotted", "dashdot")):
    """Draw each dataset's 95% nonlinear (F-test) confidence region as a
    filled contour line (not just a point marker) plus its optimum, in its
    fixed DATASET_COLORS color. Axis limits are left to the caller (autoscale
    from what was actually drawn) so this can be reused for a full view or a
    zoomed-in view without duplicating the region-computation logic."""
    for i, r in enumerate(results):
        fr = r["fr_corrected_full"]
        label = r["ds"]["label"]
        color = DATASET_COLORS[label]
        se = fr.se
        b1_grid = np.linspace(max(fr.beta[0] - 12 * se[0], 1e-3), fr.beta[0] + 12 * se[0], 150)
        b2_grid = np.linspace(max(fr.beta[1] - 12 * se[1], 1e-4), fr.beta[1] + 12 * se[1], 150)
        c_int = r["df_full"]["c_int"].to_numpy()
        c_p = r["df_full"]["c_p"].to_numpy()
        y = r["df_full"]["deltaC_m"].to_numpy()
        SSE = s3.sse_surface(b1_grid, b2_grid, c_int, c_p, y)
        thresh = s3.f_test_threshold(fr.sse, fr.n, fr.p, ALPHA)
        ax.contour(b1_grid, b2_grid, SSE.T, levels=[thresh], colors=[color], linewidths=2.5,
                   linestyles=[linestyles[i % len(linestyles)]])
        ax.plot(fr.beta[0], fr.beta[1], marker="*", color=color, markersize=12,
                markeredgecolor="k", markeredgewidth=0.5, linestyle="none",
                label=label)


def plot_all_confidence_regions(results, outfile):
    """Two-row figure: top row = one PER-DATASET zoom inset per well-identified
    run (each individually zoomed tight enough that its 95% F-test region
    renders as an actual visible ellipse, not a point), bottom row = the FULL
    view across all primary (non-ill-conditioned) runs for context.

    This was a direct fix for a PI-flagged issue: a single *shared* zoomed
    axis still has to span the full 18-62 mM range needed to place all the
    well-identified points, and each individual region (~1-2 mM wide) is
    tiny relative to that span, so it *still* renders as a dot even after
    zooming the shared axis. Only a per-dataset inset -- zoomed to each run's
    own +/-8 SE window -- makes the region visible as a region. This is
    NOT a numerical problem: the F-test region computation
    (s3.sse_surface / s3.f_test_threshold) is identical in both rows; the
    well-identified runs' regions are small because they are, in fact, well
    identified (narrow 95% CIs, e.g. MC5 05.27.26's |chi| CI is 48.6-50.5 mM)."""
    ok = [r for r in results if not r["ill_conditioned"]]
    well = [r for r in ok if is_well_identified(r)]
    poor = [r for r in ok if not is_well_identified(r)]

    w = plot_style.fig_width("full")
    n_well = max(len(well), 1)
    fig = plt.figure(figsize=(w, w * 0.75))
    gs = fig.add_gridspec(2, n_well, height_ratios=[1, 1.3], hspace=0.55, wspace=0.5)

    for j, r in enumerate(well):
        ax = fig.add_subplot(gs[0, j])
        _draw_regions(ax, [r])
        fr = r["fr_corrected_full"]
        se = fr.se
        ax.set_xlim(fr.beta[0] - 8 * se[0], fr.beta[0] + 8 * se[0])
        ax.set_ylim(fr.beta[1] - 8 * se[1], fr.beta[1] + 8 * se[1])
        ax.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]", fontsize=8)
        if j == 0:
            ax.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]", fontsize=8)
        ax.set_title(r["ds"]["label"], fontsize=8)
        ax.tick_params(labelsize=7)
        ax.xaxis.set_major_locator(plt.MaxNLocator(3))
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))

    ax_full = fig.add_subplot(gs[1, :])
    _draw_regions(ax_full, ok)
    ax_full.set_xlabel(r"$\boldsymbol{|\chi|}$ [mM]")
    ax_full.set_ylabel(r"$\boldsymbol{\delta^*}$ [–]")
    ax_full.set_title(f"Full view: all {len(ok)} runs", fontsize=9)
    ax_full.legend(loc="upper center", bbox_to_anchor=(0.5, -0.28), ncol=4,
                   fontsize=7, borderaxespad=0)

    fig.suptitle(f"95% nonlinear confidence regions across NaCl datasets\n"
                 f"(top: per-dataset zoom on the {n_well} well-identified runs, "
                 f"each panel independently scaled)", fontsize=10)
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")
    print(f"  well-identified ({len(well)}): {[r['ds']['label'] for r in well]}")
    print(f"  poorly-identified ({len(poor)}): {[r['ds']['label'] for r in poor]}")


def plot_all_fits(results, outfile):
    w = plot_style.fig_width("full")
    n = len(results)
    ncols = 3
    nrows = int(np.ceil(n / ncols))
    # Each panel roughly square (w/ncols wide) rather than the earlier
    # 0.4 factor, which squashed 3 rows into a short, hard-to-read strip
    # (PI-flagged: "each panel is wide and vertically squished").
    fig, axes = plt.subplots(nrows, ncols, figsize=(w, w * nrows / ncols * 0.95),
                              sharey=True)
    axes = np.atleast_1d(axes).flatten()
    for ax in axes[n:]:
        ax.axis("off")
    for ax, r in zip(axes, results):
        fr = r["fr_corrected_full"]
        label = r["ds"]["label"]
        color = DATASET_COLORS[label]
        df = r["df_full"]
        c_int = df["c_int"].to_numpy()
        c_p = df["c_p"].to_numpy()
        y = df["deltaC_m"].to_numpy()
        order = np.argsort(c_int)
        y_fit = step2.f_corrected(fr.beta, c_int[order], c_p[order])

        ax.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4,
                label=r"data")
        if not r["ill_conditioned"]:
            ax.plot(c_int[order], y_fit, "--", color="k", linewidth=2,
                    label=r"$f_{\mathrm{corr}}$ fit")
        else:
            ax.text(0.5, 0.9, "fit did not converge", transform=ax.transAxes,
                    ha="center", fontsize=8, color="tab:red")
        ax.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]")
        ax.set_title(label, fontsize=10)
        ax.label_outer()
    fig.supylabel(r"$\Delta \boldsymbol{c_m}$ [mM]")

    # Every panel shares the same two legend entries (marker style for
    # "data", line style for the fit) -- dataset color is already keyed to
    # the panel title, so one shared legend replaces 3 identical ones.
    shared_handles = [
        Line2D([], [], marker="o", linestyle="none", color="gray", alpha=0.6,
               label="data"),
        Line2D([], [], color="k", linestyle="--", linewidth=2,
               label=r"$f_{\mathrm{corr}}$ fit"),
    ]
    fig.legend(handles=shared_handles, loc="lower center",
               bbox_to_anchor=(0.5, -0.12), ncol=2, fontsize=9, borderaxespad=0)

    fig.suptitle("NaCl corrected-model fit vs. data, by dataset")
    fig.tight_layout()
    plot_style.save_fig(fig, outfile)
    plt.close(fig)
    print(f"Saved {outfile}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print("=" * 100)
    print("NaCl (1:1) Donnan regression across all processed datasets")
    print("=" * 100)

    print(f"APPLY_CP = {APPLY_CP}")
    rng = np.random.default_rng(RNG_SEED)
    results = []
    for ds in DATASETS:
        if ds["source"] == "excel":
            src = ds["sheet"]
        else:
            src = (ds["path_cp"] if (APPLY_CP and "path_cp" in ds) else ds["path"]).name
        print(f"--- {ds['label']} ({src}) ---")
        t0 = time.time()
        r = analyze_dataset(ds, rng)
        print(f"  done in {time.time()-t0:.1f} s")
        results.append(r)

    print()
    print_summary_table(results)

    # CP Phase 2 refinement: the primary headline figures show only the
    # uniformly-our-CP "primary" datasets (MC2, MC5(2)->relabeled "MC5
    # (05.27.26)", and the 6 raw-only vial-ICP+CP runs). MC5-main
    # (authors'-own CP) is a labeled "comparison" entry, not a second
    # primary MC5 row -- print_summary_table above still sees the full
    # `results` (it needs by_key["MC5"] for that comparison).
    plot_results = [r for r in results if r["ds"].get("role", "primary") == "primary"]

    plot_all_estimates(plot_results, FIG_DIR / "nacl_all_estimates.png")
    plot_all_confidence_regions(plot_results, FIG_DIR / "nacl_all_confidence_regions.png")
    plot_all_fits(plot_results, FIG_DIR / "nacl_all_fits.png")

    return results


if __name__ == "__main__":
    main()
