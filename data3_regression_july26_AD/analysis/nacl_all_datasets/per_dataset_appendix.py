"""
Per-dataset NaCl appendix figures: for every one of the 9 NaCl runs, a single
combined panel figure showing (as applicable):
  - Preprocessing outputs (6 reconstructed raw-only runs only): per-experiment
    conductivity calibration fit, Jw(t), and reconstructed c_int(t)/c_p(t).
  - Fit + residuals (all 9 runs): the corrected-model fit vs. deltaC_m, a
    residual-vs-c_int panel, and a text panel with (|chi|, delta*), CI, SSE, n.

Reuses analysis/nacl_all_datasets/regress_all_nacl.py's DATASETS/analyze_dataset
(no fitting logic duplicated) and analysis/preprocessing/preprocess_nacl.py for
the raw calibration data of the 6 reconstructed runs.

Run with: conda activate data3-regression && python analysis/nacl_all_datasets/per_dataset_appendix.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "analysis"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step2_nacl"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "step3_nacl_uncertainty"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "nacl_all_datasets"))
sys.path.insert(0, str(REPO_ROOT / "analysis" / "preprocessing"))
import plot_style  # noqa: E402
import regress_nacl as step2  # noqa: E402
import regress_all_nacl as ran  # noqa: E402
import preprocess_nacl as pp  # noqa: E402

plot_style.apply_style()

DATA_DIR = REPO_ROOT / "data"
FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

# Raw source (file, sheet) for each CSV-sourced (reconstructed) dataset --
# matches analysis/preprocessing/process_raw_runs.py's RUNS list.
RAW_SOURCE = {
    "MC3-0709": ("NF270_MC3.xlsx", "07.09.24_NaCl"),
    "MC3-0722S": ("NF270_MC3.xlsx", "07.22.24_SNaCl"),
    "MC4-0711S": ("NF270_MC4.xlsx", "07.11.24_SNaCl"),
    "MC5-0723": ("NF270_MC5.xlsx", "07.23.24_NaCl"),
    "MC5-0723S": ("NF270_MC5.xlsx", "07.23.24_SNaCl"),
    "MC5-0723S2": ("NF270_MC5.xlsx", "07.23.24_S2NaCl"),
}

# Output figure name per dataset key (used in main.tex).
FIG_NAME = {
    "MC2": "nacl_appendix_mc2",
    "MC5": "nacl_appendix_mc5_0527",
    "MC5(2)": "nacl_appendix_mc5_0527_2",
    "MC3-0709": "nacl_appendix_mc3_0709",
    "MC3-0722S": "nacl_appendix_mc3_0722s",
    "MC4-0711S": "nacl_appendix_mc4_0711s",
    "MC5-0723": "nacl_appendix_mc5_0723",
    "MC5-0723S": "nacl_appendix_mc5_0723s",
    "MC5-0723S2": "nacl_appendix_mc5_0723s2",
}


def _load_csv_timeseries(ds):
    """Read time_s/Jw/c_int/c_p straight from the processed CSV for
    plotting only (regress_all_nacl.load_dataset intentionally does not
    carry these columns through, to keep its NaN-filtering/regression
    columns exactly as validated -- so this is a separate, read-only path
    that does not touch or duplicate that function's behavior). Applies the
    identical validity filter (finite, c_int > c_p) so lengths/alignment
    match what was actually fit."""
    import pandas as pd
    df_raw = pd.read_csv(ds["path"], comment="#")
    valid = (
        df_raw[["c_int_mM", "c_p_mM", "Jw_m3_m2_s", "time_s"]].notna().all(axis=1)
        & (df_raw["c_int_mM"] > df_raw["c_p_mM"])
    )
    return df_raw[valid].copy()


def make_panel_reconstructed(r, ds):
    """2x3 panel: calibration fit, Jw(t), c_int/c_p(t), fit vs data,
    residual vs c_int, info text."""
    key = ds["key"]
    raw_file, raw_sheet = RAW_SOURCE[key]
    raw = pp.load_raw_sheet(DATA_DIR / raw_file, raw_sheet)
    calib = pp.fit_calibration(raw.vials)

    df = r["df_full"]
    fr = r["fr_corrected_full"]
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    ts_df = _load_csv_timeseries(ds)
    time_s = ts_df["time_s"].to_numpy()

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(2, 3, figsize=(w, w * 0.62))

    # (A) calibration fit
    ax = axes[0, 0]
    if calib is not None:
        vdf = raw.vials.dropna(subset=["cond", "icp_mgL", "sample_vol_mL", "acid_vol_mL"])
        dilution = (vdf["sample_vol_mL"] + vdf["acid_vol_mL"]) / vdf["sample_vol_mL"]
        c_mM = (vdf["icp_mgL"].to_numpy() * dilution.to_numpy()) / pp.ICP_MOLAR_MASS
        cond = vdf["cond"].to_numpy()
        ax.scatter(c_mM, cond, color=plot_style.OKABE_ITO["blue"], s=25, zorder=3)
        c_grid = np.linspace(0, c_mM.max() * 1.05, 50)
        ax.plot(c_grid, calib.intercept + calib.slope * c_grid, "k--", lw=1.5)
        ax.set_title(f"(a) calibration ($\\boldsymbol{{R^2}}$={calib.r2:.4f}, n={calib.n_points})", fontsize=8)
    else:
        ax.text(0.5, 0.5, "calibration unavailable\n(insufficient vial data)",
                ha="center", va="center", transform=ax.transAxes, fontsize=8)
        ax.set_title("(a) calibration", fontsize=8)
    ax.set_xlabel(r"$\boldsymbol{c}$ [mM]", fontsize=8)
    ax.set_ylabel(r"cond. [$\boldsymbol{\mu}$S/cm]", fontsize=8)
    ax.tick_params(labelsize=7)

    # (B) Jw(t)
    ax = axes[0, 1]
    ax.plot(time_s, ts_df["Jw_m3_m2_s"].to_numpy(), color=plot_style.OKABE_ITO["vermillion"], lw=1)
    ax.set_xlabel(r"$\boldsymbol{t}$ [s]", fontsize=8)
    ax.set_ylabel(r"$\boldsymbol{J_w}$ [m/s]", fontsize=8)
    ax.set_title(r"(b) water flux $\boldsymbol{J_w(t)}$", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

    # (C) c_int(t), c_p(t)
    ax = axes[0, 2]
    ax.plot(time_s, ts_df["c_int_mM"].to_numpy(), color=plot_style.OKABE_ITO["blue"], lw=1.2, label=r"$c_{int}$")
    ax.plot(time_s, ts_df["c_p_mM"].to_numpy(), color=plot_style.OKABE_ITO["orange"], lw=1.2, label=r"$c_p$")
    ax.legend(fontsize=7, loc="best")
    ax.set_xlabel(r"$\boldsymbol{t}$ [s]", fontsize=8)
    ax.set_ylabel(r"$\boldsymbol{c}$ [mM]", fontsize=8)
    ax.set_title(r"(c) reconstructed $\boldsymbol{c_{int}(t)}, \boldsymbol{c_p(t)}$", fontsize=8)
    ax.tick_params(labelsize=7)

    _panel_fit_residual_info(axes[1, 0], axes[1, 1], axes[1, 2], r, ds, c_int, c_p, y)

    fig.suptitle(ds["label"], fontsize=11)
    fig.tight_layout()
    return fig


def make_panel_processed(r, ds):
    """1x3 panel (no calibration/Jw/c(t) to show -- arrived pre-processed):
    fit vs data, residual vs c_int, info text."""
    df = r["df_full"]
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 3, figsize=(w, w * 0.32))
    _panel_fit_residual_info(axes[0], axes[1], axes[2], r, ds, c_int, c_p, y)
    fig.suptitle(f"{ds['label']} (arrived pre-processed -- no calibration to show)",
                 fontsize=11)
    fig.tight_layout()
    return fig


def _panel_fit_residual_info(ax_fit, ax_res, ax_info, r, ds, c_int, c_p, y):
    fr = r["fr_corrected_full"]
    label = ds["label"]
    color = ran.DATASET_COLORS[label]

    order = np.argsort(c_int)
    if not r["ill_conditioned"]:
        y_fit = step2.f_corrected(fr.beta, c_int[order], c_p[order])
        y_fit_all = step2.f_corrected(fr.beta, c_int, c_p)
        resid = y - y_fit_all
    else:
        y_fit = None
        resid = None

    ax_fit.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4, label="data")
    if y_fit is not None:
        ax_fit.plot(c_int[order], y_fit, "k--", lw=1.5, label=r"$f_{\mathrm{corr}}$ fit")
    else:
        ax_fit.text(0.5, 0.85, "fit did not converge", transform=ax_fit.transAxes,
                     ha="center", fontsize=7, color="tab:red")
    ax_fit.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]", fontsize=8)
    ax_fit.set_ylabel(r"$\Delta \boldsymbol{c_m}$ [mM]", fontsize=8)
    ax_fit.set_title("fit vs. data", fontsize=8)
    ax_fit.legend(fontsize=6, loc="best")
    ax_fit.tick_params(labelsize=7)

    if resid is not None:
        ax_res.axhline(0, color="k", lw=1, ls=":")
        ax_res.plot(c_int, resid, "o", color=color, markersize=3, alpha=0.5)
    else:
        ax_res.text(0.5, 0.5, "n/a", ha="center", va="center",
                     transform=ax_res.transAxes, fontsize=8)
    ax_res.set_xlabel(r"$\boldsymbol{c_{int}}$ [mM]", fontsize=8)
    ax_res.set_ylabel(r"residual [mM]", fontsize=8)
    ax_res.set_title(r"residual vs. $\boldsymbol{c_{int}}$", fontsize=8)
    ax_res.tick_params(labelsize=7)

    ax_info.axis("off")
    if r["ill_conditioned"]:
        info = (f"n = {fr.n}\n\nILL-CONDITIONED\n"
                f"(diverges to nonphysical\noptimum; see Sec. 2.4)\n\nSSE = {fr.sse:.3g}")
    else:
        lo1, hi1 = r["profile_cis"][0]
        lo2, hi2 = r["profile_cis"][1]
        info = (f"$n$ = {fr.n}\n\n"
                f"$|\\chi|$ = {fr.beta[0]:.3g} mM\n"
                f"95% CI ({lo1:.3g}, {hi1:.3g})\n\n"
                f"$\\delta^*$ = {fr.beta[1]:.4g}\n"
                f"95% CI ({lo2:.3g}, {hi2:.3g})\n\n"
                f"SSE = {fr.sse:.4g}\n"
                f"multi-start: {r['pct_global']:.0f}% global")
    ax_info.text(0.05, 0.95, info, transform=ax_info.transAxes, fontsize=8,
                  va="top", ha="left", family="monospace")


def main():
    rng = np.random.default_rng(ran.RNG_SEED)
    for ds in ran.DATASETS:
        print(f"--- {ds['label']} ---")
        r = ran.analyze_dataset(ds, rng)
        if ds["source"] == "csv":
            fig = make_panel_reconstructed(r, ds)
        else:
            fig = make_panel_processed(r, ds)
        outfile = FIG_DIR / f"{FIG_NAME[ds['key']]}.png"
        plot_style.save_fig(fig, outfile)
        plt.close(fig)
        print(f"  Saved {outfile}")


if __name__ == "__main__":
    main()
