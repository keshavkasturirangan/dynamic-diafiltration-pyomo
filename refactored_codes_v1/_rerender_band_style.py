#!/usr/bin/env python3
"""Re-render the per-band / per-vial overlays in the canonical figure_s5
contour-LINE style (matching DATA3_collaborator_meeting_2026-06-11.pptx), by
routing the contour through the library's own `_plot_heatmap_frame` and overlaying
the band/window theta markers.  Reads already-saved JSON + the existing
matlab_ports_5ch_full sigma x Lp CSVs — NO refits / recompute.

Run: python3 _rerender_band_style.py
"""
import sys, json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import refactored_ucb_library as lib  # noqa: E402

ART = HERE.parent / "UnifiedFramework" / "DATA3" / "results" / "paper_artifacts" / "nf270"
OUT = ART / "band_value_test"
ZCOL = "Obj_retentate_concentration"


def _contour_df(rid):
    p = ART / "matlab_ports_5ch_full" / rid / "contourdata-x_sigma-y_Lp.csv"
    return pd.read_csv(p) if p.exists() else None


def rerender_band_overlay(r):
    rid = r["rid"]
    df = _contour_df(rid)
    fig, ax = plt.subplots(figsize=(5.6, 5.2))
    if df is not None:
        lib._plot_heatmap_frame(df, "sigma", "Lp", ZCOL, ax=ax, show_title=False)
    # band theta markers (figure_s5 marker idiom: bold edge, distinct shapes)
    pts = [
        ("full sheet", r["full"], "white", "o"),
        (f"low-cF {r['cf_range']['low'][0]:.0f}-{r['cf_range']['low'][1]:.0f}mM", r["low"], "#ff5555", "v"),
        (f"high-cF {r['cf_range']['high'][0]:.0f}-{r['cf_range']['high'][1]:.0f}mM", r["high"], "#4da6ff", "s"),
    ]
    for label, row, color, mk in pts:
        if row.get("sigma") is not None and row.get("Lp") is not None:
            ax.plot([row["sigma"]], [row["Lp"]], mk, markersize=13, markeredgecolor="k",
                    markerfacecolor=color, markeredgewidth=1.4, clip_on=False, zorder=6,
                    label=f"{label}: σ={row['sigma']:.2f} Lp={row['Lp']:.2f} B={row['B']:.2f}")
    ax.set_title(f"{r['role']} · {rid}\nper-band θ on σ×Lp contour (retentate-conc objective)", fontsize=10)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), fontsize=7, frameon=False)
    fig.tight_layout()
    out = OUT / f"band_overlay-{rid}.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


def rerender_per_vial(r):
    if "skip" in r:
        return
    rid = r["rid"]
    df = _contour_df(rid)
    pts = r["points"]
    cfs = np.array([p["cf_mean"] for p in pts])
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.7))
    # (a) sigma & B vs window cF
    ax[0].plot(cfs, [p["sigma"] for p in pts], "o-", color="#d62728", label="σ")
    ax[0].set_xlabel("window mean terminal cF (mM)", fontweight="bold")
    ax[0].set_ylabel("σ", color="#d62728", fontweight="bold"); ax[0].set_ylim(-0.05, 1.05)
    ax[0].tick_params(direction="in")
    axb = ax[0].twinx(); axb.plot(cfs, [p["B"] for p in pts], "s--", color="#1f77b4")
    axb.set_ylabel("B (μm/s)", color="#1f77b4", fontweight="bold")
    ax[0].set_title(f"{r['role']} {rid}\nper-window θ  (corr cF·σ={r['corr_cf_sigma']:.2f}, cF·B={r['corr_cf_B']:.2f})", fontsize=9)
    # (b) Lp vs cF
    ax[1].plot(cfs, [p["Lp"] for p in pts], "^-", color="#2ca02c")
    ax[1].set_xlabel("window mean terminal cF (mM)", fontweight="bold")
    ax[1].set_ylabel("Lp", fontweight="bold"); ax[1].set_title("Lp vs concentration window")
    ax[1].tick_params(direction="in")
    # (c) window optima on the canonical contour, colored by cF
    if df is not None:
        lib._plot_heatmap_frame(df, "sigma", "Lp", ZCOL, ax=ax[2], show_title=False)
    sc = ax[2].scatter([p["sigma"] for p in pts], [p["Lp"] for p in pts], c=cfs, cmap="plasma",
                       s=130, edgecolors="k", linewidths=1.2, zorder=6)
    fig.colorbar(sc, ax=ax[2], label="window cF (mM)")
    ax[2].set_title("window optima on σ×Lp contour\n(color = concentration)", fontsize=10)
    fig.tight_layout()
    out = OUT / f"per_vial_trend-{rid}.png"
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {out}")


if __name__ == "__main__":
    band = json.load(open(OUT / "band_value_test_summary.json"))
    print("[restyle] band overlays:")
    for r in band:
        if "error" not in r:
            rerender_band_overlay(r)
    pv_path = OUT / "per_vial_trend_summary.json"
    if pv_path.exists():
        print("[restyle] per-vial trends:")
        for r in json.load(open(pv_path)):
            rerender_per_vial(r)
    print("[restyle] DONE")
