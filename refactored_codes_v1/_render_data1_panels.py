#!/usr/bin/env python3
"""Re-render the DATA1 three-slice identifiability set in the consistent v7 line-contour
style, per representative salt, from the existing contourdata CSVs (pure plotting — no solves).

Story per salt (the user's framing):
  Row 1  σ×Lp  (B pinned)   ┐  Lp-identifiability evidence: a tight horizontal Lp band
  Row 2  B×Lp  (σ pinned)   ┘  (closed in Lp, flat in σ/B) ⇒ Lp is well-determined per salt.
  Row 3  σ×B   (Lp FIXED)   →  the headline: residual σ–B correlation once Lp is pinned.
Columns = the three measured channels: mass | permeate-conc | retentate-conc.

CSV Obj_* columns are ALREADY log10(WSSE) (from _nf270_contour_grid_dataframe), so plot directly.
"""
import sys, json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
NF270 = HERE.parent / "UnifiedFramework/DATA3/results/paper_artifacts/nf270"
SRC = NF270 / "contour_panels"                       # 3-channel CSVs, complete for all 4 reps
OUT = NF270 / "bform_study" / "Lp_identifiability_panels"

REPS = [("MC3.07.22.24_SNaCl", "NaCl — diluting"), ("MC2.05.07.24_NaCl", "NaCl — concentrating"),
        ("MC2.05.07.24_CaCl2", "CaCl₂"), ("MC2.05.21.24_LaCl3", "LaCl₃")]
PLANES = [("x_sigma-y_Lp", "sigma", "Lp", "σ×Lₚ  (B pinned)"),
          ("x_B-y_Lp",     "B",     "Lp", "B×Lₚ  (σ pinned)"),
          ("x_sigma-y_B",  "sigma", "B",  "σ×B   (Lₚ FIXED at fit)")]
CHAN = [("Obj_mass", "Mass"), ("Obj_concentration", "Permeate-conc"), ("Obj_retentate_concentration", "Retentate-conc")]
LAB = {"Lp": "L$_p$ [L m$^{-2}$h$^{-1}$bar$^{-1}$]", "sigma": "$\\sigma$ [-]", "B": "B [$\\mu$m s$^{-1}$]"}


def _grid(df, xn, yn, zn):
    xs = np.sort(df[xn].unique()); ys = np.sort(df[yn].unique())
    Z = np.full((len(ys), len(xs)), np.nan)
    xi = {v: i for i, v in enumerate(xs)}; yi = {v: i for i, v in enumerate(ys)}
    for _, r in df.iterrows():
        z = r[zn]
        if np.isfinite(z):
            Z[yi[r[yn]], xi[r[xn]]] = z
    return xs, ys, Z


def _panel(ax, xs, ys, Z, title, xlab, ylab):
    ax.set_facecolor("white"); ax.set_title(title, fontsize=8.5)
    ax.set_xlabel(xlab, fontsize=8); ax.set_ylabel(ylab, fontsize=8)
    if Z is None or np.all(np.isnan(Z)):
        ax.text(0.5, 0.5, "(no data)", transform=ax.transAxes, ha="center", color="0.5"); return
    X, Y = np.meshgrid(xs, ys)
    try:
        cs = ax.contour(X, Y, Z, levels=10, cmap="turbo", linewidths=1.3)   # Z already log10
        ax.clabel(cs, inline=True, fontsize=5.5, fmt="%.2f")
    except Exception:
        pass
    ax.grid(alpha=0.2, lw=0.4)
    jm, im = np.unravel_index(np.nanargmin(Z), Z.shape)
    ax.plot(xs[im], ys[jm], "^", ms=10, color="red", mec="white", mew=1.0, zorder=6)
    ax.text(0.5, -0.30, f"min ({xs[im]:.2f}, {ys[jm]:.2f}) log₁₀={Z[jm, im]:.2f}",
            transform=ax.transAxes, ha="center", fontsize=6, color="red")


def render(rid, label):
    sd = SRC / rid
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 11.5))
    for ri, (tag, xn, yn, ptitle) in enumerate(PLANES):
        csv = sd / f"contourdata-{tag}.csv"
        if not csv.exists():
            for ci in range(3):
                axes[ri, ci].text(0.5, 0.5, f"(missing {tag})", transform=axes[ri, ci].transAxes, ha="center", color="0.6")
            continue
        df = pd.read_csv(csv)
        for ci, (zn, cname) in enumerate(CHAN):
            xs, ys, Z = _grid(df, xn, yn, zn)
            _panel(axes[ri, ci], xs, ys, Z, f"{ptitle}\n{cname} residual² (log₁₀)", LAB[xn], LAB[yn])
    fig.suptitle(f"{label}   ({rid})   ·   DATA1-style identifiability slices (DATA3 experiment): "
                 f"Lₚ pinned by mass (rows 1–2) ⇒ focus on σ–B at fixed Lₚ (row 3)", fontsize=12, y=0.997)
    fig.tight_layout(rect=[0, 0, 1, 0.985])
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{rid}_identifiability.png"
    fig.savefig(out, dpi=140); plt.close(fig)
    print(f"[render] {label:22s} -> {out}")
    return out


if __name__ == "__main__":
    for rid, label in REPS:
        render(rid, label)
    print(f"[render] done -> {OUT}")
