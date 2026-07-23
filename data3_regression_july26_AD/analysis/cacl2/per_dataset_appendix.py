"""
Per-dataset CaCl2 appendix figures: for each of the 3 CaCl2 runs, a single
combined panel figure showing the fit-vs-data, residual-vs-c_int, and
effective-|chi|-vs-c_int diagnostics plus a parameter/CI text summary, so a
collaborator can dig into each experiment individually (mirrors
analysis/nacl_all_datasets/per_dataset_appendix.py's role for NaCl).

Reuses analysis/cacl2/regress_cacl2.py's DATASETS/analyze_dataset (no
fitting/diagnostic logic duplicated).

Run with: conda activate data3-regression && python analysis/cacl2/per_dataset_appendix.py
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
sys.path.insert(0, str(REPO_ROOT / "analysis" / "cacl2"))
import plot_style  # noqa: E402
import regress_cacl2 as m  # noqa: E402

plot_style.apply_style()

FIG_DIR = REPO_ROOT / "docs" / "reports" / "figures"
FIG_DIR.mkdir(parents=True, exist_ok=True)

FIG_NAME = {
    "MC3-0711S": "cacl2_appendix_mc3_0711s",
    "MC5-0527": "cacl2_appendix_mc5_0527",
    "MC5-0627": "cacl2_appendix_mc5_0627",
}


def make_panel(r):
    ds = r["ds"]
    fr = r["fr"]
    label = ds["label"]
    color = m.DATASET_COLORS[label]
    df = r["df_full"]
    c_int = df["c_int"].to_numpy()
    c_p = df["c_p"].to_numpy()
    y = df["deltaC_m"].to_numpy()
    order = np.argsort(c_int)
    y_fit = m.f_cacl2(fr.beta, c_int[order], c_p[order])

    w = plot_style.fig_width("full")
    fig, axes = plt.subplots(1, 4, figsize=(w, w * 0.30))

    ax = axes[0]
    ax.plot(c_int, y, "o", color=color, markersize=3, alpha=0.4, label="data")
    ax.plot(c_int[order], y_fit, "--", color="k", linewidth=1.5, label=r"$f_{2:1}$ fit")
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
        ax.errorbar(centers[rel], locals_[rel], yerr=ses[rel], fmt="o-",
                    color=color, markersize=3, capsize=2, elinewidth=1)
    if np.any(~rel):
        ax.plot(centers[~rel], locals_[~rel], "x", color=color, markersize=4,
                alpha=0.35)
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

    fig.suptitle(label, fontsize=11)
    fig.tight_layout(w_pad=2.0)
    return fig


def main():
    rng = np.random.default_rng(m.RNG_SEED)
    for ds in m.DATASETS:
        print(f"--- {ds['label']} ---")
        r = m.analyze_dataset(ds, rng)
        fig = make_panel(r)
        outfile = FIG_DIR / f"{FIG_NAME[ds['key']]}.png"
        plot_style.save_fig(fig, outfile)
        plt.close(fig)
        print(f"  Saved {outfile}")


if __name__ == "__main__":
    main()
